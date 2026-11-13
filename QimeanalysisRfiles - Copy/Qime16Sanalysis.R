# 1. LOAD LIBRARIES
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("phyloseq", "biomformat")) # Run if needed
library(phyloseq)
library(tidyverse)
library(scales)

# 2. SET WORKING DIRECTORY (Adjust path if necessary)
setwd("C:/Users/iaadeoye/Downloads/QimeanalysisRfiles")

# 3. IMPORT DATA
# Load OTU table (handles the tricky QIIME 2 header)
otu_df <- read.table("feature-table.tsv", 
                     sep = "\t", 
                     header = TRUE, 
                     skip = 1, 
                     comment.char = "", 
                     check.names = FALSE)
rownames(otu_df) <- otu_df[,1]
otu_mat <- as.matrix(otu_df[,-1])

# Correct for R's tendency to change dashes to dots
colnames(otu_mat) <- gsub("\\.", "-", colnames(otu_mat))

# Load Taxonomy and split into ranks
tax_raw <- read_tsv("taxonomy.tsv")
tax_mat <- tax_raw %>%
  rename(ASV = 1) %>%
  # Fix: Remove the space after the semicolon
  separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
           sep = ";", fill = "right") %>%
  # Fix: Clean prefixes (d__, p__, etc.) more safely
  mutate(across(everything(), ~str_replace_all(., "^[a-z]__", ""))) %>%
  mutate(across(everything(), ~na_if(., ""))) %>%
  mutate(across(everything(), ~replace_na(., "Unknown"))) %>%
  column_to_rownames("ASV") %>%
  select(-any_of("Confidence")) %>% # Remove Confidence column if it exists
  as.matrix()

# Load Metadata and Tree
metadata <- read_tsv("metadata.tsv") %>% rename(SampleID = `sample-id`)
tree     <- read_tree("tree.nwk")

# 4. CREATE PHYLOSEQ OBJECT
# Find samples that exist in both files
common_samples <- intersect(colnames(otu_mat), metadata$SampleID)

ps <- phyloseq(
  otu_table(otu_mat[, common_samples], taxa_are_rows = TRUE),
  tax_table(tax_mat),
  sample_data(metadata %>% filter(SampleID %in% common_samples) %>% column_to_rownames("SampleID")),
  phy_tree(tree)
)
#save.image("step_1-4_filter_trim.RData")
#load("step_1-4_filter_trim.RData")

################################################################################
###############################################################################
###############################################################################
# 5. FILTER: Bulkwater, Second-Run, Weeks 1-7 Genus Level
target_weeks <- paste("Week", 1:7)

ps_filtered <- subset_samples(
  ps,
  SampleCondition == "Second-Run" &
    SampleType %in% c("Bulkwater-Control", "Bulkwater-Spiked") &
    SampleWeek %in% target_weeks
)

# 6. RELATIVE ABUNDANCE AT GENUS LEVEL
ps_rel <- transform_sample_counts(ps_filtered, function(x) x / sum(x))
glom <- tax_glom(ps_rel, taxrank = "Genus", NArm = FALSE)

# 7. PREPARE DATA FRAME
plot_df_all <- psmelt(glom) %>%
  mutate(
    # 1. 
    SampleType = recode(SampleType, 
                        "Bulkwater-Control" = "Bulk Water-Control", 
                        "Bulkwater-Spiked"  = "Bulk Water-Spiked"),
    
    SampleWeek = factor(SampleWeek, levels = target_weeks),
    
    # 2. 
    SampleType = factor(SampleType, levels = c("Bulk Water-Control", "Bulk Water-Spiked")),
    
    Genus = ifelse(is.na(Genus) | Genus == "", "Unknown", as.character(Genus))
  )

# 8. SELECT TOP 15 (EXCLUDING UNKNOWN & INCERTAE_SEDIS)
top15 <- plot_df_all %>%
  filter(!Genus %in% c("Unknown", "Incertae_Sedis")) %>%
  group_by(Genus) %>%
  summarise(Total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  slice_head(n = 15) %>%
  pull(Genus)

# 9. GROUP UNKNOWN, INCERTAE_SEDIS AND NON-TOP GENERA AS OTHER
plot_df <- plot_df_all %>%
  mutate(Genus = ifelse(Genus %in% top15, Genus, "Other")) %>%
  group_by(SampleWeek, SampleType, Genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# 10. ORDER GENERA, KEEP OTHER LAST
genus_order <- plot_df %>%
  group_by(Genus) %>%
  summarise(Total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  pull(Genus)

genus_order <- c(setdiff(genus_order, "Other"), "Other")
plot_df$Genus <- factor(plot_df$Genus, levels = genus_order)

# 11. PALETTE
pal <- c(
  "#4682B4", "#CD5C5C", "#00FF00", "#FF00FF",
  "#00FFFF", "#FFFF00", "#1869FF", "#800000",
  "#008000", "#800080", "#008080", "#808000",
  "#808080", "#F5F5DC", "#FFA500", "#CCCCCC"
)

# 12. PLOT
p <- ggplot(plot_df, aes(x = SampleType, y = Abundance, fill = Genus)) +
  geom_col(position = "stack", color = NA) +
  facet_wrap(~SampleWeek, nrow = 1) +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = pal) +
  labs(
    x = "Sample Type",
    y = "Relative Abundance (%)",
    #title = "Top 15 Genera: Bulkwater Weekly Samples"
  ) +
  theme_classic(base_size = 38) +
  theme(
    axis.title.x = element_text(face = "bold", size = 38, margin = margin(t = 20)),
    axis.title.y = element_text(face = "bold", size = 38, margin = margin(r = 20)),
    axis.text.x = element_text(size = 38, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 38),
    legend.title = element_text(face = "bold", size = 34),
    legend.text = element_text(size = 38),
    legend.key.size = unit(1.5, "cm"),
    plot.title = element_text(size = 38, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 34, face = "bold"),
    strip.background = element_rect(fill = "grey95", color = "black", linewidth = 1.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.8),
    panel.spacing.x = unit(0.15, "cm")
  )
print(p)

# 13. EXPORT CSV
genus_table <- plot_df %>%
  arrange(SampleWeek, SampleType, desc(Abundance))

write.csv(
  genus_table,
  "Taxonomy_Top15_Genus_Bulkwater_SecondRun_ByWeek.csv",
  row.names = FALSE
)

# 14. SAVE OUTPUTS
ggsave("Taxonomy_Top15_Genus_Bulkwater_SecondRun_ByWeek.png", p, width = 22, height = 14, dpi = 300)
ggsave("Taxonomy_Top15_Genus_Bulkwater_SecondRun_ByWeek.pdf", p, width = 22, height = 14)
#######################################################################################
#####################################################################################
#################################################################################
# 5. FILTER: Bulkwater, Weeks 1-7 by FAMILY Level
target_weeks <- paste("Week", 1:7)
ps_filtered <- subset_samples(
  ps, 
  SampleCondition == "Second-Run" & 
    SampleType %in% c("Bulkwater-Control", "Bulkwater-Spiked") & 
    SampleWeek %in% target_weeks
)

# 6. RELATIVE ABUNDANCE AT FAMILY LEVEL
ps_rel <- transform_sample_counts(ps_filtered, function(x) x / sum(x))
glom <- tax_glom(ps_rel, taxrank = "Family", NArm = FALSE)

# 7. PREPARE DATA FRAME 
plot_df_all <- psmelt(glom) %>% 
  mutate( 
    # 1.  
    SampleType = recode(SampleType, "Bulkwater-Control" = "Bulk Water-Control", "Bulkwater-Spiked" = "Bulk Water-Spiked"), 
    SampleWeek = factor(SampleWeek, levels = target_weeks), 
    # 2.  
    SampleType = factor(SampleType, levels = c("Bulk Water-Control", "Bulk Water-Spiked")), 
    
    
    Family = ifelse(is.na(Family) | Family == "", "Unknown", as.character(Family)) 
  ) 

# 8. SELECT TOP 15 (EXCLUDING UNKNOWN & INCERTAE_SEDIS) 
top15 <- plot_df_all %>% 
  filter(!Family %in% c("Unknown", "Incertae_Sedis")) %>% 
  group_by(Family) %>% 
  summarise(Total = sum(Abundance), .groups = "drop") %>% 
  arrange(desc(Total)) %>% 
  slice_head(n = 15) %>% 
  pull(Family)

# 9. GROUP UNKNOWN, INCERTAE_SEDIS AND NON-TOP FAMILIES AS OTHER
plot_df <- plot_df_all %>%
  mutate(Family = ifelse(Family %in% top15, Family, "Other")) %>%
  group_by(SampleWeek, SampleType, Family) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# 10. ORDER FAMILIES, KEEP OTHER LAST
family_order <- plot_df %>%
  group_by(Family) %>%
  summarise(Total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  pull(Family)

family_order <- c(setdiff(family_order, "Other"), "Other")
plot_df$Family <- factor(plot_df$Family, levels = family_order)

# 11. PALETTE
pal <- c(
  "#008A00", "#1869FF", "#00FF00", "#FF00FF", "#00FFFF", "#FFFF00", 
  "#000080", "#800000", "#008000", "#800080", "#008080", "#7600A1", 
  "#808080", "#F5F5DC", "#D21820", "#CCCCCC"
)

# 12. PLOT
p <- ggplot(plot_df, aes(x = SampleType, y = Abundance, fill = Family)) +
  geom_col(position = "stack", color = NA) +
  facet_wrap(~SampleWeek, nrow = 1) +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = pal) +
  labs(
    x = "Sample Type",
    y = "Relative Abundance (%)",
    fill = "Family"
  ) +
  theme_classic(base_size = 30) +
  theme(
    axis.title.x = element_text(face = "bold", size = 32, margin = margin(t = 20)),
    axis.title.y = element_text(face = "bold", size = 32, margin = margin(r = 20)),
    axis.text.x = element_text(size = 28, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 28),
    legend.title = element_text(face = "bold", size = 32),
    legend.text = element_text(size = 28),
    legend.key.size = unit(1.5, "cm"),
    strip.text = element_text(size = 30, face = "bold"),
    strip.background = element_rect(fill = "grey95", color = "black", linewidth = 1.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.8),
    panel.spacing.x = unit(0.15, "cm")
  )

print(p)

# 13. EXPORT CSV
family_table <- plot_df %>% 
  arrange(SampleWeek, SampleType, desc(Abundance))

write.csv(
  family_table, 
  "Taxonomy_Top15_Family_Bulkwater_SecondRun_ByWeek.csv", 
  row.names = FALSE
)

# 14. SAVE OUTPUTS
ggsave("Taxonomy_Top15_Family_Bulkwater_SecondRun_ByWeek.png", p, width = 28, height = 14, dpi = 600)
ggsave("Taxonomy_Top15_Family_Bulkwater_SecondRun_ByWeek.pdf", p, width = 28, height = 14)

################################################################################
###############################################################################
###############################################################################
# 5. FILTER: Bulkwater, Second-Run, Weeks 8-10 at Genus Level
target_weeks <- paste("Week", 8:10)

ps_filtered <- subset_samples(
  ps,
  SampleCondition == "Second-Run" &
    SampleType %in% c("Bulkwater-Control", "Bulkwater-Spiked") &
    SampleWeek %in% target_weeks
)

# 6. RELATIVE ABUNDANCE AT GENUS LEVEL
ps_rel <- transform_sample_counts(ps_filtered, function(x) x / sum(x))
glom <- tax_glom(ps_rel, taxrank = "Genus", NArm = FALSE)

# 7. PREPARE DATA FRAME
plot_df_all <- psmelt(glom) %>%
  mutate(
    # 1. 
    SampleType = recode(SampleType, 
                        "Bulkwater-Control" = "Bulk Water-Control", 
                        "Bulkwater-Spiked"  = "Bulk Water-Spiked"),
    
    SampleWeek = factor(SampleWeek, levels = target_weeks),
    
    # 2. 
    SampleType = factor(SampleType, levels = c("Bulk Water-Control", "Bulk Water-Spiked")),
    
    Genus = ifelse(is.na(Genus) | Genus == "", "Unknown", as.character(Genus))
  )

# 8. SELECT TOP 15 (EXCLUDING UNKNOWN & INCERTAE_SEDIS)
top15 <- plot_df_all %>%
  filter(!Genus %in% c("Unknown", "Incertae_Sedis")) %>%
  group_by(Genus) %>%
  summarise(Total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  slice_head(n = 15) %>%
  pull(Genus)

# 9. GROUP UNKNOWN, INCERTAE_SEDIS AND NON-TOP GENERA AS OTHER
plot_df <- plot_df_all %>%
  mutate(Genus = ifelse(Genus %in% top15, Genus, "Other")) %>%
  group_by(SampleWeek, SampleType, Genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# 10. ORDER GENERA, KEEP OTHER LAST
genus_order <- plot_df %>%
  group_by(Genus) %>%
  summarise(Total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  pull(Genus)

genus_order <- c(setdiff(genus_order, "Other"), "Other")
plot_df$Genus <- factor(plot_df$Genus, levels = genus_order)

# 11. PALETTE
pal <- c(
  "#008A00", "#00FF00", "#938290", "#CD5C5C",
  "#4682B4", "#088800", "#000080", "#800000",
  "#F5F5DC", "#808000", "#463524", "#845323",
  "#808080", "#231213", "#008080", "#CCCCCC"
)

# 12. PLOT
p <- ggplot(plot_df, aes(x = SampleType, y = Abundance, fill = Genus)) +
  geom_col(position = "stack", color = NA) +
  facet_wrap(~SampleWeek, nrow = 1) +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = pal) +
  labs(
    x = "Sample Type",
    y = "Relative Abundance (%)",
    #title = "Top 15 Genera: Bulkwater after coupon swap"
  ) +
  theme_classic(base_size = 38) +
  theme(
    axis.title.x = element_text(face = "bold", size = 38, margin = margin(t = 20)),
    axis.title.y = element_text(face = "bold", size = 38, margin = margin(r = 20)),
    axis.text.x = element_text(size = 38, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 38),
    legend.title = element_text(face = "bold", size = 38),
    legend.text = element_text(size = 38),
    legend.key.size = unit(1.5, "cm"),
    plot.title = element_text(size = 38, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 34, face = "bold"),
    strip.background = element_rect(fill = "grey95", color = "black", linewidth = 1.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.8),
    panel.spacing.x = unit(0.15, "cm")
  )
print(p)

# 13. EXPORT CSV
genus_table <- plot_df %>%
  arrange(SampleWeek, SampleType, desc(Abundance))

write.csv(
  genus_table,
  "Taxonomy_Top15_Genus_Bulkwater_SecondRun_ByWeek.csv",
  row.names = FALSE
)

# 14. SAVE OUTPUTS
ggsave("Taxonomy_Top15_Genus_Bulkwater_SecondRun_Week8-10.png", p, width = 16, height = 14, dpi = 300)
ggsave("Taxonomy_Top15_Genus_Bulkwater_SecondRun_Week8-10.pdf", p, width = 16, height = 14)
################################################################################
###############################################################################
###############################################################################
# 5. FILTER: Biofilm groups, Second-Run Genus Level
target_weeks <- paste("Week", 1:7)

ps_filtered <- subset_samples(
  ps,
  SampleCondition == "Second-Run" &
    SampleType %in% c("Biofilm-Control", "Biofilm-Spiked") &
    SampleWeek %in% target_weeks
)

# 6. GENUS-LEVEL RELATIVE ABUNDANCE
glom <- tax_glom(ps_filtered, taxrank = "Genus", NArm = FALSE)
ps_rel <- transform_sample_counts(glom, function(x) x / sum(x))

# 7. PREPARE DATA
plot_df_all <- psmelt(ps_rel) %>%
  mutate(
    SampleWeek = factor(SampleWeek, levels = target_weeks),
    SampleType = factor(SampleType, levels = c("Biofilm-Control", "Biofilm-Spiked")),
    Genus = ifelse(is.na(Genus) | Genus == "", "Unknown", as.character(Genus))
  )

# 8. SELECT TOP 15, EXCLUDING UNKNOWN AND INCERTAE_SEDIS
top15 <- plot_df_all %>%
  filter(!Genus %in% c("Unknown", "Incertae_Sedis")) %>%
  group_by(Genus) %>%
  summarise(Total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  slice_head(n = 15) %>%
  pull(Genus)

# 9. COLLAPSE REMAINING GENERA INTO OTHER
plot_df <- plot_df_all %>%
  mutate(Genus = ifelse(Genus %in% top15, Genus, "Other")) %>%
  group_by(SampleWeek, SampleType, Genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# 10. ORDER GENERA, KEEP OTHER LAST
genus_order <- plot_df %>%
  group_by(Genus) %>%
  summarise(Total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  pull(Genus)

genus_order <- c(setdiff(genus_order, "Other"), "Other")
plot_df$Genus <- factor(plot_df$Genus, levels = genus_order)

# 11. WEEK 1 EMPTY PANEL WITH LABELS

week1_blank <- data.frame(
  SampleWeek = factor(c("Week 1", "Week 1"), levels = target_weeks),
  SampleType = factor(c("Biofilm-Control", "Biofilm-Spiked"),
                      levels = c("Biofilm-Control", "Biofilm-Spiked")),
  y = 1
)
week1_note <- data.frame(
  SampleWeek = factor("Week 1", levels = target_weeks),
  x = 1.5,
  y = 0.5,
  Label = "No biofilm\nformed"
)

# 12. COLOR PALETTE
pal <- c("#008A00","#1869FF","#D21820","#FFB000","#7600A1","#00BFC4",
         "#FF6F00","#8B4513","#F564E3","#7CAE00","#C77CFF","#00A08A",
         "#E76BF3","#A3A500","#4682B4","#CCCCCC")

# 13. PLOT
p <- ggplot(plot_df, aes(x = SampleType, y = Abundance, fill = Genus)) +
  geom_col(position = "stack", color = NA) +
  geom_blank(data = week1_blank, aes(x = SampleType, y = 1), inherit.aes = FALSE) +
  geom_text(data = week1_note, aes(x = x, y = y, label = Label),
            inherit.aes = FALSE, size = 10, fontface = "bold") +
  facet_wrap(~SampleWeek, nrow = 1, drop = FALSE) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1), expand = c(0, 0)) +
  scale_fill_manual(values = pal) + 
  labs(x = "Sample Type", y = "Relative Abundance (%)") + 
  theme_classic(base_size = 38)+
  theme(
    axis.title.x = element_text(face = "bold", size = 38, margin = margin(t = 20)),
    axis.title.y = element_text(face = "bold", size = 38, margin = margin(r = 20)),
    axis.text.x = element_text(size = 38, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 38),
    legend.title = element_text(face = "bold", size = 36),
    legend.text = element_text(size = 38),
    legend.key.size = unit(1.5, "cm"),
    plot.title = element_text(size = 38, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 34, face = "bold"),
    strip.background = element_rect(fill = "grey95", color = "black", linewidth = 1.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.8),
    panel.spacing.x = unit(0.15, "cm")
  )

print(p)

# 14. EXPORT TABLE
genus_table <- plot_df %>% arrange(SampleWeek, SampleType, desc(Abundance))
write.csv(genus_table, "Taxonomy_Top15_Genus_Biofilm_SecondRun_ByWeek.csv", row.names = FALSE)
# 15. SAVE OUTPUTS
ggsave("Taxonomy_Top15_Genus_Biofilm_SecondRun_ByWeek.png", p, width = 22, height = 14, dpi = 300)
ggsave("Taxonomy_Top15_Genus_Biofilm_SecondRun_ByWeek.pdf", p, width = 22, height = 14)
################################################################################
###############################################################################
################################################################################
# 5. FILTER: Biofilm groups, Second-Run Family Level
target_weeks <- paste("Week", 1:7)
ps_filtered <- subset_samples(
  ps, 
  SampleCondition == "Second-Run" & 
    SampleType %in% c("Biofilm-Control", "Biofilm-Spiked") & 
    SampleWeek %in% target_weeks
)

# 6. FAMILY-LEVEL RELATIVE ABUNDANCE
glom <- tax_glom(ps_filtered, taxrank = "Family", NArm = FALSE)
ps_rel <- transform_sample_counts(glom, function(x) x / sum(x))

# 7. PREPARE DATA
plot_df_all <- psmelt(ps_rel) %>%
  mutate(
    SampleWeek = factor(SampleWeek, levels = target_weeks),
    SampleType = factor(SampleType, levels = c("Biofilm-Control", "Biofilm-Spiked")),
    Family = ifelse(is.na(Family) | Family == "", "Unknown", as.character(Family))
  )

# 8. SELECT TOP 15, EXCLUDING UNKNOWN AND INCERTAE_SEDIS
top15 <- plot_df_all %>%
  filter(!Family %in% c("Unknown", "Incertae_Sedis")) %>%
  group_by(Family) %>%
  summarise(Total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  slice_head(n = 15) %>%
  pull(Family)

# 9. COLLAPSE REMAINING FAMILIES INTO OTHER
plot_df <- plot_df_all %>%
  mutate(Family = ifelse(Family %in% top15, Family, "Other")) %>%
  group_by(SampleWeek, SampleType, Family) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# 10. ORDER FAMILIES, KEEP OTHER LAST
family_order <- plot_df %>%
  group_by(Family) %>%
  summarise(Total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  pull(Family)

family_order <- c(setdiff(family_order, "Other"), "Other")
plot_df$Family <- factor(plot_df$Family, levels = family_order)

# 11. WEEK 1 EMPTY PANEL WITH LABELS
week1_blank <- data.frame(
  SampleWeek = factor(c("Week 1", "Week 1"), levels = target_weeks),
  SampleType = factor(c("Biofilm-Control", "Biofilm-Spiked"), levels = c("Biofilm-Control", "Biofilm-Spiked")),
  y = 1
)

week1_note <- data.frame(
  SampleWeek = factor("Week 1", levels = target_weeks),
  x = 1.5,
  y = 0.5,
  Label = "No biofilm\nformed"
)

# 12. COLOR PALETTE
pal <- c(
  "#008A00", "#1869FF", "#D21820", "#FFB000", "#7600A1", "#00BFC4", 
  "#FF6F00", "#8B4513", "#F564E3", "#7CAE00", "#C77CFF", "#00A08A", 
  "#E76BF3", "#A3A500", "#4682B4", "#CCCCCC"
)

# 13. PLOT
p <- ggplot(plot_df, aes(x = SampleType, y = Abundance, fill = Family)) +
  geom_col(position = "stack", color = NA) +
  geom_blank(data = week1_blank, aes(x = SampleType, y = 1), inherit.aes = FALSE) +
  geom_text(data = week1_note, aes(x = x, y = y, label = Label), inherit.aes = FALSE, size = 10, fontface = "bold") +
  facet_wrap(~SampleWeek, nrow = 1, drop = FALSE) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1), expand = c(0, 0)) +
  scale_fill_manual(values = pal) +
  labs(
    x = "Sample Type", 
    y = "Relative Abundance (%)", 
    fill = "Family"
  ) +
  theme_classic(base_size = 30) +
  theme(
    axis.title.x = element_text(face = "bold", size = 32, margin = margin(t = 20)),
    axis.title.y = element_text(face = "bold", size = 32, margin = margin(r = 20)),
    axis.text.x = element_text(size = 28, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 28),
    legend.title = element_text(face = "bold", size = 32),
    legend.text = element_text(size = 28),
    legend.key.size = unit(1.5, "cm"),
    strip.text = element_text(size = 30, face = "bold"),
    strip.background = element_rect(fill = "grey95", color = "black", linewidth = 1.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.8),
    panel.spacing.x = unit(0.15, "cm")
  )

print(p)

# 14. EXPORT TABLE
family_table <- plot_df %>% 
  arrange(SampleWeek, SampleType, desc(Abundance))

write.csv(
  family_table, 
  "Taxonomy_Top15_Family_Biofilm_SecondRun_ByWeek.csv", 
  row.names = FALSE
)

# 15. SAVE OUTPUTS
ggsave("Taxonomy_Top15_Family_Biofilm_SecondRun_ByWeek.png", p, width = 28, height = 14, dpi = 600)
ggsave("Taxonomy_Top15_Family_Biofilm_SecondRun_ByWeek.pdf", p, width = 28, height = 14)
################################################################################
###############################################################################
###############################################################################
# 5. FILTER: Biofilm, Second-Run, Weeks 8–10
target_weeks <- paste("Week", 8:10)

ps_filtered <- subset_samples(
  ps,
  SampleCondition == "Second-Run" &
    SampleType %in% c("Biofilm-Control", "Biofilm-Spiked") &
    SampleWeek %in% target_weeks
)

# 6. GENUS-LEVEL RELATIVE ABUNDANCE
glom <- tax_glom(ps_filtered, taxrank = "Genus", NArm = FALSE)
ps_rel <- transform_sample_counts(glom, function(x) x / sum(x))

# 7. PREPARE DATA
plot_df_all <- psmelt(ps_rel) %>%
  mutate(
    SampleWeek = factor(SampleWeek, levels = target_weeks),
    SampleType = factor(SampleType, levels = c("Biofilm-Control", "Biofilm-Spiked")),
    Genus = ifelse(is.na(Genus) | Genus == "", "Unknown", as.character(Genus))
  )

# 8. SELECT TOP 15 GENERA, EXCLUDING UNKNOWN AND INCERTAE_SEDIS
top15 <- plot_df_all %>%
  filter(!Genus %in% c("Unknown", "Incertae_Sedis")) %>%
  group_by(Genus) %>%
  summarise(Total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  slice_head(n = 15) %>%
  pull(Genus)

# 9. COLLAPSE ALL REMAINING GENERA INTO OTHER
plot_df <- plot_df_all %>%
  mutate(Genus = ifelse(Genus %in% top15, Genus, "Other")) %>%
  group_by(SampleWeek, SampleType, Genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# 10. ORDER GENERA AND KEEP OTHER LAST
genus_order <- plot_df %>%
  group_by(Genus) %>%
  summarise(Total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  pull(Genus)

plot_df$Genus <- factor(plot_df$Genus, levels = c(setdiff(genus_order, "Other"), "Other"))

# 11. COLOR PALETTE
pal <- c("#008A00","#1869FF","#F564E3","#808000","#FFB000","#7600A1",
         "#FF6F00","#00BFC4","#FFB000","#CD5C5C","#C77CFF","#00A08A",
         "#7CAE00","#8B4513","#4682B4","#CCCCCC")

# 12. PLOT
p <- ggplot(plot_df, aes(SampleType, Abundance, fill = Genus)) +
  geom_col() +
  facet_wrap(~SampleWeek, nrow = 1, drop = FALSE) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1), expand = c(0, 0)) +
  scale_fill_manual(values = pal) + 
  labs(x = "Sample Type", y = "Relative Abundance (%)") + 
  theme_classic(base_size = 38)+
  theme(
    axis.title.x = element_text(face = "bold", size = 38, margin = margin(t = 20)),
    axis.title.y = element_text(face = "bold", size = 38, margin = margin(r = 20)),
    axis.text.x = element_text(size = 38, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 38),
    legend.title = element_text(face = "bold", size = 32),
    legend.text = element_text(size = 38),
    legend.key.size = unit(1.5, "cm"),
    plot.title = element_text(size = 38, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 34, face = "bold"),
    strip.background = element_rect(fill = "grey95", color = "black", linewidth = 1.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.8),
    panel.spacing.x = unit(0.15, "cm")
  )

print(p)

# 13. EXPORT DATA AND FIGURES
genus_table <- plot_df %>% arrange(SampleWeek, SampleType, desc(Abundance))

write.csv(genus_table, "Taxonomy_Top15_Genus_Biofilm_SecondRun_Week8-10.csv", row.names = FALSE)
ggsave("Taxonomy_Top15_Genus_Biofilm_SecondRun_Week8-10.png", p, width = 18, height = 14, dpi = 600)
ggsave("Taxonomy_Top15_Genus_Biofilm_SecondRun_Week8-10.pdf", p, width = 18, height = 14)
################################################################################
###############################################################################
###############################################################################
# 5. FILTER: Influent, Second-Run, Weeks 1-10
target_weeks <- paste("Week", 1:10)

ps_filtered <- subset_samples(
  ps,
  SampleCondition == "Second-Run" &
    SampleType == "Influent" &
    SampleWeek %in% target_weeks
)

# 6. GENUS-LEVEL RELATIVE ABUNDANCE
glom <- tax_glom(ps_filtered, taxrank = "Genus", NArm = FALSE)
ps_rel <- transform_sample_counts(glom, function(x) x / sum(x))

# 7. PREPARE DATA
plot_df_all <- psmelt(ps_rel) %>%
  mutate(
    SampleWeek = factor(SampleWeek, levels = target_weeks),
    SampleType = factor(SampleType, levels = "Influent"),
    Genus = ifelse(is.na(Genus) | Genus == "", "Unknown", as.character(Genus))
  )

# 8. SELECT TOP 15, EXCLUDING UNKNOWN AND INCERTAE_SEDIS
top15 <- plot_df_all %>%
  filter(!Genus %in% c("Unknown", "Incertae_Sedis")) %>%
  group_by(Genus) %>%
  summarise(Total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  slice_head(n = 15) %>%
  pull(Genus)

# 9. COLLAPSE REMAINING GENERA INTO OTHER
plot_df <- plot_df_all %>%
  mutate(Genus = ifelse(Genus %in% top15, Genus, "Other")) %>%
  group_by(SampleWeek, SampleType, Genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# 10. ORDER GENERA, KEEP OTHER LAST
genus_order <- plot_df %>%
  group_by(Genus) %>%
  summarise(Total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  pull(Genus)

genus_order <- c(setdiff(genus_order, "Other"), "Other")
plot_df$Genus <- factor(plot_df$Genus, levels = genus_order)

# 11. COLOR PALETTE
pal <- c("#567467", "#785743", "#980", "#CD5C5C", "#00FFFF",
         "#FFFF00", "#000080", "#800000", "#008000", "#800080",
         "#E76BF3","#008080","#F5F5FF","#485746", "#FFA500","#CCCCCC")

# 12. PLOT
p <- ggplot(plot_df, aes(x = SampleWeek, y = Abundance, fill = Genus)) +
  geom_col(position = "stack", color = NA) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1), expand = c(0, 0)) +
  scale_fill_manual(values = pal) +
  labs(x = "Sample Week", y = "Relative Abundance (%)") +
  theme_classic(base_size = 32) +
  theme(
    axis.title.x = element_text(face = "bold", size = 32, margin = margin(t = 20)),
    axis.title.y = element_text(face = "bold", size = 32, margin = margin(r = 20)),
    axis.text.x = element_text(size = 28, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 28),
    legend.title = element_text(face = "bold", size = 32),
    legend.text = element_text(size = 28),
    legend.key.size = unit(1.5, "cm"),
    plot.title = element_text(size = 36, face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.8)
  )

print(p)

# 13. EXPORT TABLE
genus_table <- plot_df %>% arrange(SampleWeek, desc(Abundance))
write.csv(genus_table, "Taxonomy_Top15_Genus_Influent_SecondRun_ByWeek.csv", row.names = FALSE)

# 14. SAVE OUTPUTS
ggsave("Taxonomy_Top15_Genus_Influent_SecondRun_ByWeek.png", p, width = 20, height = 14, dpi = 600)
ggsave("Taxonomy_Top15_Genus_Influent_SecondRun_ByWeek.pdf", p, width = 20, height = 14)
################################################################################
###############################################################################
#field biofilm samples
library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(scales)

# 1. Filter Week 0 field biofilm samples
ps_filtered <- subset_samples(ps, SampleCondition == "Second-Run" &
                                SampleType == "Field-Biofilm" & SampleWeek == "Week0")
ps_filtered <- prune_taxa(taxa_sums(ps_filtered) > 0, ps_filtered)
stopifnot(nsamples(ps_filtered) == 6)

# 2. Relative abundance and genus aggregation
ps_rel <- transform_sample_counts(ps_filtered, function(x) x / sum(x))
ps_genus <- tax_glom(ps_rel, taxrank = "Genus", NArm = FALSE)

# 3. Prepare data and extract replicate numbers
plot_df_all <- psmelt(ps_genus) %>%
  mutate(Replicate = factor(paste("Rep", str_extract(Sample, "\\d+$")), levels = paste("Rep", 1:6)),
         SampleWeek = "Week 0",
         Genus = ifelse(is.na(Genus) | Genus == "", "Unknown", as.character(Genus)))

# 4. Identify top 15 genera
top15 <- plot_df_all %>%
  filter(!Genus %in% c("Unknown", "Incertae_Sedis")) %>%
  group_by(Genus) %>% summarise(Total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(Total)) %>% slice_head(n = 15) %>% pull(Genus)

# 5. Group remaining genera as Other
plot_df <- plot_df_all %>%
  mutate(Genus = ifelse(Genus %in% top15, Genus, "Other")) %>%
  group_by(SampleWeek, Replicate, Genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

genus_order <- c(top15, "Other")
plot_df$Genus <- factor(plot_df$Genus, levels = genus_order)

# 6. Color palette
pal <- c("#FFA500", "#485746", "#567467", "#F5F5FF", "#785743","#CD5C5C","#000080",
         "#800000", "#FFFF00", "#008000", "#800080","#008080",
         "#E76BF3","#980","#00FFFF","#CCCCCC")

names(pal) <- genus_order

# 7. Plot
p <- ggplot(plot_df, aes(x = Replicate, y = Abundance, fill = Genus)) +
  geom_col(width = 0.85) +
  facet_grid(~SampleWeek, labeller = labeller(SampleWeek = c("Week 0" = "Sample"))) +
  scale_y_continuous(labels = percent_format(), breaks = seq(0, 1, 0.25),
                     limits = c(0, 1), expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values = pal, breaks = genus_order) +
  labs(x = "Field Biofilm Replicates", y = "Relative Abundance (%)", fill = "Genus") +
  theme_classic(base_size = 30) +
  theme(
    axis.title = element_text(size = 30, face = "bold"),
    axis.text = element_text(size = 30, color = "black"),
    strip.text = element_text(size = 30, face = "bold"),
    strip.background = element_rect(fill = "grey95", color = "black"),
    legend.position = "right",
    legend.title = element_text(size = 30, face = "bold"),
    legend.text = element_text(size = 30),
    legend.key.size = grid::unit(0.6, "cm"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.margin = margin(12, 15, 12, 12)
  )

print(p)

# 8. Verify each replicate sums to 100%
print(plot_df %>% group_by(Replicate) %>%
        summarise(Total = sum(Abundance), .groups = "drop"))

# 9. Export results
write.csv(plot_df, "Taxonomy_Top15_Genus_FieldBiofilm_Week0.csv", row.names = FALSE)
ggsave("Taxonomy_Top15_Genus_FieldBiofilm_Week0.png", p, width = 17, height = 9, dpi = 300, bg = "white")
ggsave("Taxonomy_Top15_Genus_FieldBiofilm_Week0.pdf", p, width = 17, height = 9, bg = "white")

################################################################################
###############################################################################
###############################################################################
#family_tree_top30.png
library(phyloseq)
library(tidyverse)
library(ape)
library(ggplot2)

load("step_1-4_filter_trim.RData")

ps <- prune_taxa(taxa_sums(ps) > 0, ps)

taxmat <- as(tax_table(ps), "matrix")
taxmat[, "Family"] <- ifelse(is.na(taxmat[, "Family"]) | taxmat[, "Family"] == "", "Unknown", taxmat[, "Family"])
tax_table(ps) <- tax_table(taxmat)

ps_family <- tax_glom(ps, taxrank = "Family", NArm = FALSE)

top_families <- names(sort(taxa_sums(ps_family), decreasing = TRUE))[1:30]
ps_family_top <- prune_taxa(top_families, ps_family)

p_family_tree <- plot_tree(
  ps_family_top,
  color = "Family",
  label.tips = "Family",
  size = "abundance",
  ladderize = TRUE
) +
  theme_bw() +
  ggtitle("Phylogenetic Tree of Top 30 Families")

print(p_family_tree)

ggsave("family_tree_top30.png", p_family_tree, width = 14, height = 10, dpi = 300)

save(ps_family, ps_family_top, file = "family_tree_objects.RData")
################################################################################
###############################################################################
###############################################################################
# ============================================================
# TARGET SPIKED ORGANISMS GENERA LINE PLOT
# Biofilm-Control vs Biofilm-Spiked, Weeks 1 to 10 without greybackround for wks 8-10
# ============================================================
OUT_DIR <- "C:/Users/iaadeoye/Downloads/QimeanalysisRfiles"

ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
ps_genus_rel <- transform_sample_counts(ps_genus, function(x) x / sum(x))
genus_df <- psmelt(ps_genus_rel)

# Define the master ordered chronological sequence explicitly
weeks_order <- c("Week 1", "Week 2", "Week 3", "Week 4", "Week 5", "Week 6", "Week 7", "Week 8", "Week 9", "Week 10")

plot_df <- genus_df %>%
  filter(SampleType %in% c("Biofilm-Control", "Biofilm-Spiked")) %>%
  filter(SampleWeek %in% weeks_order) %>%
  mutate(
    # CRITICAL FIX: Convert to factor IMMEDIATELY here to override alphabetical sorting traps
    SampleWeek = factor(SampleWeek, levels = weeks_order, ordered = TRUE),
    TargetGenus = case_when(
      Genus %in% c("Escherichia/Shigella", "Escherichia-Shigella", "Escherichia") ~ "Escherichia",
      Genus == "Enterococcus" ~ "Enterococcus",
      Genus == "Pseudomonas" ~ "Pseudomonas",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(TargetGenus)) %>%
  group_by(SampleType, SampleWeek, TargetGenus) %>%
  summarise(RelAbund = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  # Fill structurally missing slots using our strict weeks_order definition
  complete(
    SampleType = factor(c("Biofilm-Control", "Biofilm-Spiked"), levels = c("Biofilm-Control", "Biofilm-Spiked")),
    SampleWeek = factor(weeks_order, levels = weeks_order, ordered = TRUE),
    TargetGenus = factor(c("Escherichia", "Enterococcus", "Pseudomonas"), levels = c("Escherichia", "Enterococcus", "Pseudomonas")),
    fill = list(RelAbund = NA)
  )

p_target_genera <- ggplot(
  plot_df, 
  aes(x = SampleWeek, y = RelAbund, group = SampleType, color = SampleType, shape = SampleType, linetype = SampleType)
) +
  geom_line(linewidth = 1.2, na.rm = TRUE) +
  geom_point(size = 3.8, na.rm = TRUE) +
  
  # Text remains bound to Week 1 at the baseline of each free facet
  annotate(
    geom = "text",
    x = "Week 1", 
    y = -Inf, 
    label = "No visible biofilm formed",
    angle = 90,
    hjust = -0.05,
    vjust = 0.5,
    color = "black",
    size = 6.5
  ) +
  
  facet_wrap(~ TargetGenus, nrow = 1, scales = "free_y") +
  scale_color_manual(values = c("Biofilm-Control" = "#333333", "Biofilm-Spiked" = "#1F78B4")) +
  scale_shape_manual(values = c("Biofilm-Control" = 16, "Biofilm-Spiked" = 17)) +
  scale_linetype_manual(values = c("Biofilm-Control" = "dashed", "Biofilm-Spiked" = "solid")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01, scale = 100)) +
  labs(
    title = "Spiked organisms tracking in Biofilm-Control and Biofilm-Spiked samples",
    x = "Sample Week",
    y = "Relative abundance (%)",
    color = NULL,
    shape = NULL,
    linetype = NULL
  ) +
  theme_bw(base_size = 24) +
  theme(
    plot.title = element_text(size = 28, face = "plain", hjust = 0.5),
    axis.title.x = element_text(size = 26, face = "bold", margin = margin(t = 20)),
    axis.title.y = element_text(size = 26, face = "bold", margin = margin(r = 20)),
    axis.text.x = element_text(size = 22, angle = 35, hjust = 1),
    axis.text.y = element_text(size = 20),
    strip.text = element_text(size = 24, face = "bold"),
    strip.background = element_rect(fill = "grey95", color = "black"),
    legend.position = "right",
    legend.text = element_text(size = 20),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

print(p_target_genera)

ggsave(file.path(OUT_DIR, "Target_Genera_LinePlot_RelAbund.png"), p_target_genera, width = 19, height = 7, dpi = 600)
ggsave(file.path(OUT_DIR, "Target_Genera_LinePlot_RelAbund.pdf"), p_target_genera, width = 19, height = 7)
write_csv(plot_df, file.path(OUT_DIR, "Target_Genera_LinePlot_Data.csv"))
################################################################################
################################################################################
################################################################################
###############################################################################
# ============================================================ 
# TARGET SPIKED ORGANISMS GENERA LINE PLOT 
# Biofilm-Control vs Biofilm-Spiked, Weeks 1 to 10 with Grey Background for wks 8-10
# ============================================================ 
OUT_DIR <- "C:/Users/iaadeoye/Downloads/QimeanalysisRfiles" 

ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE) 
ps_genus_rel <- transform_sample_counts(ps_genus, function(x) x / sum(x)) 
genus_df <- psmelt(ps_genus_rel) 

# Define the master ordered chronological sequence explicitly 
weeks_order <- c("Week 1", "Week 2", "Week 3", "Week 4", "Week 5", "Week 6", "Week 7", "Week 8", "Week 9", "Week 10") 

plot_df <- genus_df %>% 
  filter(SampleType %in% c("Biofilm-Control", "Biofilm-Spiked")) %>% 
  filter(SampleWeek %in% weeks_order) %>% 
  mutate( 
    # This is to Convert to factor IMMEDIATELY here to override alphabetical sorting traps 
    SampleWeek = factor(SampleWeek, levels = weeks_order, ordered = TRUE), 
    TargetGenus = case_when( 
      Genus %in% c("Escherichia/Shigella", "Escherichia-Shigella", "Escherichia") ~ "Escherichia", 
      Genus == "Enterococcus" ~ "Enterococcus", 
      Genus == "Pseudomonas" ~ "Pseudomonas", 
      TRUE ~ NA_character_ 
    ) 
  ) %>% 
  filter(!is.na(TargetGenus)) %>% 
  group_by(SampleType, SampleWeek, TargetGenus) %>% 
  summarise(RelAbund = sum(Abundance, na.rm = TRUE), .groups = "drop") %>% 
  # Fill structurally missing slots using our strict weeks_order definition 
  complete( 
    SampleType = factor(c("Biofilm-Control", "Biofilm-Spiked"), levels = c("Biofilm-Control", "Biofilm-Spiked")), 
    SampleWeek = factor(weeks_order, levels = weeks_order, ordered = TRUE), 
    TargetGenus = factor(c("Escherichia", "Enterococcus", "Pseudomonas"), levels = c("Escherichia", "Enterococcus", "Pseudomonas")), 
    fill = list(RelAbund = NA) 
  ) 

p_target_genera <- ggplot( 
  plot_df, 
  aes(x = SampleWeek, y = RelAbund, group = SampleType, color = SampleType, shape = SampleType, linetype = SampleType) 
) + 
  # Background Shading layer for Weeks 8 to 10
  # xmin = 7.5 drops the vertical edge exactly halfway between Week 7 and Week 8
  annotate(
    "rect",
    xmin = 7.5, xmax = Inf, ymin = -Inf, ymax = Inf,
    fill = "grey92", alpha = 1
  ) +
  geom_line(linewidth = 1.2, na.rm = TRUE) + 
  geom_point(size = 3.8, na.rm = TRUE) + 
  
  # To make text remain bound to Week 1 at the baseline of each free facet 
  annotate( 
    geom = "text", 
    x = "Week 1", 
    y = -Inf, 
    label = "No visible biofilm formed", 
    angle = 90, 
    hjust = -0.05, 
    vjust = 0.5, 
    color = "black", 
    size = 10 
  ) + 
  facet_wrap(~ TargetGenus, nrow = 1, scales = "free_y") + 
  scale_color_manual(values = c("Biofilm-Control" = "#333333", "Biofilm-Spiked" = "#1F78B4")) + 
  scale_shape_manual(values = c("Biofilm-Control" = 16, "Biofilm-Spiked" = 17)) + # Keeps original 16/17 shapes
  scale_linetype_manual(values = c("Biofilm-Control" = "dashed", "Biofilm-Spiked" = "solid")) + 
  # Strips out "Week " to leave only numbers on the x-axis
  scale_x_discrete(labels = function(x) gsub("Week ", "", x)) + 
  
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01, scale = 100)) +
  # ========================================================================= 
  labs( 
    #title = "Spiked organisms tracking in Biofilm-Control and Biofilm-Spiked samples", 
    x = "Sample Week", 
    y = "Relative abundance (%)", 
    color = NULL, 
    shape = NULL, 
    linetype = NULL 
  ) + 
  theme_bw(base_size = 32) + 
  theme( 
    plot.title = element_text(size = 32, face = "plain", hjust = 0.5), 
    axis.title.x = element_text(size = 32, face = "bold", margin = margin(t = 20)), 
    axis.title.y = element_text(size = 32, face = "bold", margin = margin(r = 20)), 
    axis.text.x = element_text(size = 32, angle = 0, hjust = 0.5, vjust = 0.5, color = "black"), 
    axis.text.y = element_text(size = 32, color = "black"), 
    strip.text = element_text(size = 32, face = "bold"), 
    strip.background = element_rect(fill = "grey95", color = "black"), 
    legend.position = "right", 
    legend.text = element_text(size = 32), 
    legend.title = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8) 
  ) 

print(p_target_genera) 

ggsave(file.path(OUT_DIR, "Target_Genera_LinePlot_RelAbundgreybackground.png"), p_target_genera, width = 20.5, height = 8, dpi = 600) 
ggsave(file.path(OUT_DIR, "Target_Genera_LinePlot_RelAbundgreybackground.pdf"), p_target_genera, width = 20.5, height = 8) 
write_csv(plot_df, file.path(OUT_DIR, "Target_Genera_LinePlot_Datagreybackround.csv"))
################################################################################
################################################################################
################################################################################
# # ============================================================ 
# # TARGET SPIKED ORGANISMS GENERA LINE PLOT 
# # Bulkwater-Control vs Bulkwater-Spiked, Weeks 1 to 10
# # ============================================================ 
OUT_DIR <- "C:/Users/iaadeoye/Downloads/QimeanalysisRfiles"

ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE) 
ps_genus_rel <- transform_sample_counts(ps_genus, function(x) x / sum(x)) 
genus_df <- psmelt(ps_genus_rel) 

# Define explicit 1 to 10 sequence limits
weeks_order <- c("Week 1", "Week 2", "Week 3", "Week 4", "Week 5", "Week 6", "Week 7", "Week 8", "Week 9", "Week 10")

plot_df <- genus_df %>% 
  filter(SampleType %in% c("Bulkwater-Control", "Bulkwater-Spiked")) %>% 
  filter(SampleWeek %in% weeks_order) %>% 
  mutate( 
    # CRITICAL FIX: Convert SampleWeek to a factor IMMEDIATELY here to block alphabetical sorting traps
    SampleWeek = factor(SampleWeek, levels = weeks_order, ordered = TRUE),
    
    TargetGenus = case_when( 
      Genus %in% c("Escherichia/Shigella", "Escherichia-Shigella", "Escherichia") ~ "Escherichia", 
      Genus == "Enterococcus" ~ "Enterococcus", 
      Genus == "Pseudomonas" ~ "Pseudomonas", 
      TRUE ~ NA_character_ 
    ) 
  ) %>% 
  filter(!is.na(TargetGenus)) %>% 
  group_by(SampleType, SampleWeek, TargetGenus) %>% 
  summarise(RelAbund = sum(Abundance, na.rm = TRUE), .groups = "drop") %>% 
  complete( 
    SampleType = factor(c("Bulkwater-Control", "Bulkwater-Spiked"), levels = c("Bulkwater-Control", "Bulkwater-Spiked")), 
    SampleWeek = factor(weeks_order, levels = weeks_order, ordered = TRUE), 
    TargetGenus = factor(c("Escherichia", "Enterococcus", "Pseudomonas"), levels = c("Escherichia", "Enterococcus", "Pseudomonas")), 
    fill = list(RelAbund = 0) 
  ) 

p_target_genera <- ggplot( 
  plot_df, 
  aes(x = SampleWeek, y = RelAbund, group = SampleType, color = SampleType, shape = SampleType, linetype = SampleType) 
) + 
  # Shading box layer for Weeks 8 to 10
  annotate(
    "rect",
    xmin = 7.5, xmax = Inf, ymin = -Inf, ymax = Inf,
    fill = "grey92", alpha = 1
  ) +
  geom_line(linewidth = 1.2) + 
  geom_point(size = 3.8) + 
  facet_wrap(~ TargetGenus, nrow = 1, scales = "free_y") + 
  # =========================================================================
# =========================================================================
# UPDATED PLOT LAYERS (LEGEND LABELS & X-AXIS CLEANUP)
# =========================================================================
scale_color_manual(
  values = c("Bulkwater-Control" = "#333333", "Bulkwater-Spiked" = "#1F78B4"),
  labels = c("Bulkwater-Control" = "Bulk Water-Control", "Bulkwater-Spiked" = "Bulk Water-Spiked")
) +
  scale_shape_manual(
    values = c("Bulkwater-Control" = 16, "Bulkwater-Spiked" = 17),
    labels = c("Bulkwater-Control" = "Bulk Water-Control", "Bulkwater-Spiked" = "Bulk Water-Spiked")
  ) +
  scale_linetype_manual(
    values = c("Bulkwater-Control" = "dashed", "Bulkwater-Spiked" = "solid"),
    labels = c("Bulkwater-Control" = "Bulk Water-Control", "Bulkwater-Spiked" = "Bulk Water-Spiked")
  ) +
  
  # Strips out "Week " to leave only numbers on the x-axis
  scale_x_discrete(labels = function(x) gsub("Week ", "", x)) + 
  
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01, scale = 100)) +
  # =========================================================================
  labs( 
    #title = "Spiked organisms tracking in Bulkwater-Control and Bulkwater-Spiked samples", 
    x = "Sample Week", 
    y = "Relative abundance (%)", 
    color = NULL, 
    shape = NULL, 
    linetype = NULL 
  ) + 
  theme_bw(base_size = 34) + 
  theme( 
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 32, face = "plain", hjust = 0.5), 
    axis.title.x = element_text(size = 32, face = "bold", margin = margin(t = 20)), 
    axis.title.y = element_text(size = 32, face = "bold", margin = margin(r = 20)), 
    axis.text.x = element_text(size = 32, angle = 0, hjust = 0.5, vjust = 0.5, color = "black"), 
    axis.text.y = element_text(size = 32, color = "black"), 
    strip.text = element_text(size = 32, face = "bold"), 
    strip.background = element_rect(fill = "grey95", color = "black"), 
    legend.position = "right", 
    legend.text = element_text(size = 32), 
    legend.title = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8) 
  ) 

print(p_target_genera) 

ggsave(file.path(OUT_DIR, "Target_Genera_LinePlot_RelAbund_Bulkwater.png"), p_target_genera, width = 20.5, height = 8, dpi = 600) 
ggsave(file.path(OUT_DIR, "Target_Genera_LinePlot_RelAbund_Bulkwater.pdf"), p_target_genera, width = 20.5, height = 8) 
write_csv(plot_df, file.path(OUT_DIR, "Target_Genera_LinePlot_Data_Bulkwater.csv"))

################################################################################
###############################################################################
###############################################################################
# =============================================================================
# 7. PUBLICATION-QUALITY NMDS PLOT (UPDATED)
# =============================================================================
library(ggplot2)
library(dplyr)

# Define shapes exactly matching the image (using updated publication labels)
sample_shapes <- c(
  "Biofilm-Control"    = 21, # Circle
  "Biofilm-Spiked"     = 24, # Up-triangle
  "Bulk Water-Control" = 22, # Square
  "Bulk Water-Spiked"  = 23, # Diamond
  "Biofilm-Field"      = 25, # Down-triangle
  "Influent"           = 8   # Asterisk (*)
)

week_cols <- c("#3B0F70","#623F8D","#365C8D","#277F8E","#1FA187","#4AC16D", "#A0DA39","#FDE725","#FFB000","#F87600","#D7191C")


nmds_df <- nmds_df %>%
  mutate(
    SampleType = case_match(as.character(SampleType),
                            "Bulkwater-Control" ~ "Bulk Water-Control",
                            "Bulkwater-Spiked"  ~ "Bulk Water-Spiked",
                            "Field-Biofilm"     ~ "Biofilm-Field",
                            .default            = as.character(SampleType)
    ),
    SampleType = factor(SampleType, levels = names(sample_shapes))
  )

p_nmds <- ggplot(nmds_df, aes(NMDS1, NMDS2)) + 
  # 1. Plot all points (Color and Fill are both mapped to Week so asterisks change color by week)
  geom_point(aes(shape = SampleType, color = Week, fill = Week), size = 12, stroke = 1.7, alpha = 0.95) + 
  
  # 2. Text labels (including Influent)
  geom_text(
    aes(label = Week), 
    color = "black", 
    size = 8, 
    fontface = "bold"
  ) + 
  
  # 3. Apply configurations and styling
  scale_shape_manual(name = "Sample Type", values = sample_shapes) + 
  scale_fill_gradientn(name = "Week", colors = colorRampPalette(week_cols)(11), limits = range(nmds_df$Week, na.rm = TRUE), breaks = sort(unique(nmds_df$Week))) + 
  scale_color_gradientn(name = "Week", colors = colorRampPalette(week_cols)(11), limits = range(nmds_df$Week, na.rm = TRUE), breaks = sort(unique(nmds_df$Week))) + 
  labs(
    #title = "NMDS of Microbial Communities", 
    subtitle = paste0("Bray–Curtis dissimilarity · Stress = ", sprintf("%.3f", ord_nmds$stress)), 
    x = "NMDS1", 
    y = "NMDS2"
  ) + 
  coord_equal() + 
  
  # 4. Clean up legend appearances to match the reference image
  guides(
    shape = guide_legend(order = 1, override.aes = list(size = 9, fill = "black", color = "black")), 
    # UPDATED: Changed to guide_colorbar for a smooth gradient, and added reverse = TRUE
    fill = guide_colorbar(
      order = 2, 
      barheight = unit(10, "cm"), # Adjusts how tall the continuous bar is
      barwidth = unit(1.1, "cm")  # Matches your theme's legend width
    ), 
    color = "none" # Hides the duplicate color legend entry
  ) + 
  theme_classic(base_size = 36) + 
  theme(
    plot.title = element_text(face = "bold", size = 26, hjust = 0.5), 
    plot.subtitle = element_text(size = 28, hjust = 0.5, margin = margin(b = 14)), 
    axis.title = element_text(face = "bold", size = 34), 
    axis.text = element_text(size = 34, color = "black"), 
    legend.title = element_text(face = "bold", size = 34), 
    legend.text = element_text(size = 34), 
    legend.position = "right", 
    legend.key.size = unit(1.1, "cm"), 
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2), 
    plot.margin = margin(15, 20, 15, 15)
  )
print(p_nmds)
# Save the final outputs
ggsave("NMDS_BrayCurtis_SecondRun_Final.png", p_nmds, width = 16, height = 11, dpi = 600, bg = "white")
ggsave("NMDS_BrayCurtis_SecondRun_Final.pdf", p_nmds, width = 16, height = 11)
#save.image("step_nmds.RData")
#load("step_NMDS.RData")
################################################################################
###############################################################################
##############################################################################
# =============================================================================
# GLOBAL PERMANOVA: DETAILED TIME VS. TREATMENT BREAKDOWN
# =============================================================================
library(vegan)
library(phyloseq)

# 1. Standardize matrix structure (using your filtered ps_norm)
global_otu <- as(otu_table(ps_norm), "matrix")
if (taxa_are_rows(ps_norm)) {
  global_otu <- t(global_otu)
}

# 2. Build Bray-Curtis distance and format metadata
global_bc <- vegdist(global_otu, method = "bray")
global_metadata <- data.frame(sample_data(ps_norm), check.names = FALSE)

# 3. Ensure variables are treated cleanly as categorical factor blocks
global_metadata$SampleType <- as.factor(global_metadata$SampleType)
global_metadata$SampleWeek <- as.factor(global_metadata$SampleWeek)

# 4. Run the split multi-factor PERMANOVA (using sequential term margins)
set.seed(123)
global_permanova_split <- adonis2(global_bc ~ SampleWeek + SampleType, 
                                  data = global_metadata, 
                                  permutations = 999, 
                                  by = "terms") # <--- Explicitly forces individual rows

# 5. Output the individual factor breakdown table
cat("\n--- DETAILED GLOBAL PERMANOVA BREAKDOWN ---\n")
print(global_permanova_split)
################################################################################
#============================================================================= 
# 2. DATA PREPARATION & ALPHA DIVERSITY CALCULATION 
# ============================================================================= 
# Calculate alpha diversity measures directly from phyloseq 
alpha_raw <- estimate_richness(ps, measures = c("Observed", "Shannon")) %>% 
  rownames_to_column("SampleID") 

# Clean, isolate, and format metadata factors 
metadata_clean <- data.frame(sample_data(ps), check.names = FALSE) %>% 
  rownames_to_column("SampleID") 

# Convert dots to dashes so sample identifiers match perfectly 
alpha_raw$SampleID <- gsub("\\.", "-", alpha_raw$SampleID) 

# Merge alpha diversity metrics with cleaned metadata 
alpha_df <- left_join(alpha_raw, metadata_clean, by = "SampleID") %>% 
  mutate( 
    # Safely convert SampleWeek strings ("Week 1") into matching numeric values (1) 
    Week_Num = as.numeric(str_extract(as.character(SampleWeek), "\\d+")), 
    
    # Standardize sample names matching your production labels
    SampleType = case_match(as.character(SampleType), 
                            "Bulkwater-Control" ~ "Bulk Water-Control", 
                            "Bulkwater-Spiked" ~ "Bulk Water-Spiked", 
                            .default = as.character(SampleType))
  ) %>% 
  filter(!is.na(Week_Num)) 

# TARGET GROUPS: Expanded to include all 5 updated publication categories
target_types <- c("Biofilm-Control", "Biofilm-Spiked", "Bulk Water-Control", "Bulk Water-Spiked", "Influent") 

# Filter down to targets and expand week range constraint to 1-10
plot_data <- alpha_df %>% 
  filter( 
    SampleCondition == "Second-Run", # Keeps exactly 1 sample per week 
    SampleType %in% target_types, 
    Week_Num >= 1 & Week_Num <= 10 
  ) 

# ============================================================================= 
# 3. STATISTICAL ANALYSIS: SPEARMAN TEMPORAL TREND WITH BONFERRONI CORRECTION
# ============================================================================= 
# Clean statistical data base (ensuring values exist)
stat_data_base <- plot_data %>% 
  filter(!is.na(Week_Num), !is.na(Observed), !is.na(Shannon)) 

# Create empty vectors to store values for final correction
raw_p_list   <- c()
raw_rho_list <- c()
test_names   <- c()

# Loop 1: Run the tests and collect the raw statistics
for (group in target_types) { 
  
  if (str_detect(group, "Biofilm")) {
    group_subset <- filter(stat_data_base, SampleType == group, Week_Num >= 2 & Week_Num <= 10)
  } else {
    group_subset <- filter(stat_data_base, SampleType == group, Week_Num >= 1 & Week_Num <= 10)
  }
  
  if(nrow(group_subset) >= 3) {
    x_weeks <- as.numeric(group_subset$Week_Num) 
    y_obs   <- as.numeric(group_subset$Observed) 
    y_shan  <- as.numeric(group_subset$Shannon) 
    
    obs_test     <- cor.test(x_weeks, y_obs, method = "spearman", exact = FALSE) 
    shannon_test <- cor.test(x_weeks, y_shan, method = "spearman", exact = FALSE) 
    
    # Store raw p-values and rho estimates
    raw_p_list   <- c(raw_p_list, obs_test$p.value, shannon_test$p.value)
    raw_rho_list <- c(raw_rho_list, obs_test$estimate, shannon_test$estimate)
    
    # Store clean labels for the final table
    test_names   <- c(test_names, paste0(group, " [Observed]"), paste0(group, " [Shannon]"))
    
    # Store estimates directly into the object attributes for downstream reference if needed
    assign(paste0("rho_obs_", group), obs_test$estimate)
    assign(paste0("rho_shan_", group), shannon_test$estimate)
  } 
}

# Loop 2: Apply Bonferroni Adjustment and Print the Summary Table
names(raw_p_list) <- test_names
corrected_p_list  <- p.adjust(raw_p_list, method = "bonferroni")

# Combine Rho estimates, original p-values, and Bonferroni p-values into a single data frame
alpha_trend_summary <- data.frame(
  Rho_estimate = round(raw_rho_list, 4),
  Original_p   = round(raw_p_list, 4),
  Bonferroni_p = round(corrected_p_list, 4)
)

cat("\n=======================================================\n")
cat("       ALPHA DIVERSITY TEMPORAL TREND SUMMARY          ")
cat("\n=======================================================\n")
print(alpha_trend_summary)


# ============================================================================= 
# 4. FACETED GGPLOT LINE GRAPH GENERATION 
# ============================================================================= 
# Pivot long for layout scaling behavior across panels 
long_plot_data <- plot_data %>% 
  select(SampleID, Week_Num, SampleType, Observed, Shannon) %>% 
  pivot_longer(cols = c(Observed, Shannon), names_to = "Metric", values_to = "Diversity") %>% 
  mutate( 
    # Format structural factor names and panel assignments to match your figure headers 
    Metric = ifelse(Metric == "Observed", "Observed ASVs", "Shannon"), 
    Metric = factor(Metric, levels = c("Observed ASVs", "Shannon")), 
    SampleType = factor(SampleType, levels = target_types) 
  ) 

# Set custom aesthetics matching colors, shapes, and lines from your layout reference 
palette_cols <- c(
  "Biofilm-Control"    = "#00BA38", 
  "Biofilm-Spiked"     = "#619CFF",
  "Bulk Water-Control" = "#B79F00", 
  "Bulk Water-Spiked"  = "#F564E3", 
  "Influent"           = "#F8766D"
) 

custom_shapes <- c(
  "Biofilm-Control"    = 21, # Circle
  "Biofilm-Spiked"     = 24, # Up-triangle
  "Bulk Water-Control" = 22, # Square
  "Bulk Water-Spiked"  = 23, # Diamond
  "Influent"           = 8   # Asterisk (*)
) 

custom_lines <- c(
  "Biofilm-Control"    = "longdash", 
  "Biofilm-Spiked"     = "dashed",
  "Bulk Water-Control" = "dotdash", 
  "Bulk Water-Spiked"  = "twodash", 
  "Influent"           = "solid"
) 

p_alpha <- ggplot(long_plot_data, aes(
  x = Week_Num, 
  y = Diversity, 
  color = SampleType, 
  fill = SampleType,    # Added fill mapping so shapes 21-24 are solid-filled
  group = SampleType
)) + 
  
  # Highlights post-coupon change phase background window (Weeks 8 to 10)
  annotate("rect", xmin = 7.5, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey92", alpha = 1) +
  
  # Draw trend lines with a uniform solid line style for all groups
  geom_line(linetype = "solid", linewidth = 1.2) + 
  
  # Render point sets (stroke adds clear borders to the filled shapes)
  geom_point(aes(shape = SampleType), size = 5, stroke = 1.2) + 
  
  # Structural facet splitting with independent free y-axes scales 
  facet_wrap(~Metric, scales = "free_y") + 
  
  # Forces the y-axes for both panels to start cleanly at 0
  expand_limits(y = 0) + 
  
  # Map visual scales 
  scale_color_manual(values = palette_cols) + 
  scale_fill_manual(values = palette_cols) +   # Map fill to the same palette colors
  scale_shape_manual(values = custom_shapes) + 
  scale_x_continuous(breaks = 1:10, labels = paste("Week", 1:10)) + 
  
  # Set typography parameters matching original reference 
  labs( 
    x = "Sample Week", 
    y = "Alpha diversity" 
  ) + 
  theme_bw(base_size = 30) + 
  theme( 
    axis.title = element_text(face = "bold", size = 30), 
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"), 
    axis.text.y = element_text(color = "black"), 
    strip.text = element_text(face = "bold", size = 30), # Panel Headers 
    strip.background = element_rect(fill = "gray90", color = "black"), 
    legend.title = element_blank(), # Clean legend layout without a title header 
    legend.text = element_text(size = 28), 
    legend.position = "right", 
    panel.grid.minor = element_blank(), 
    plot.margin = margin(15, 15, 15, 15) 
  )
  # Render and save visual outputs 
  print(p_alpha) 
  ggsave("Alpha_Diversity_Trends_SecondRun.png", p_alpha, width = 20, height = 8, dpi = 600, bg = "white")

