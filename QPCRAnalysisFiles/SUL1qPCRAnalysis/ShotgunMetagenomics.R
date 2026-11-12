################################################################################
################################################################################
# Load required libraries
#Bioinformatic Analysis of Metagenomics
#Plot ARGs/ SCGs from the CARD results. tsv 
# Load required libraries
library(readxl)
library(dplyr)
library(ggplot2)

# Step 1: Load the Excel data
# Replace with your actual file path
card_data <- read_excel("CARD_Results.xlsx")
# Convert TPM_Normalization to numeric (safely)
card_data$TPM_Normalization <- as.numeric(card_data$TPM_Normalization)
# Step 2: Filter out rows with zero TPM
filtered_data <- card_data %>%
  filter(TPM_Normalization > 0)

# Step 3: Summarize TPM per ARG gene across all samples
top_genes <- filtered_data %>%
  group_by(gene) %>%
  summarise(Total_TPM = sum(TPM_Normalization, na.rm = TRUE)) %>%
  arrange(desc(Total_TPM)) %>%
  slice_max(Total_TPM, n = 15)

# Step 4: Plot top 15 ARGs by total TPM
ggplot(top_genes, aes(x = reorder(gene, Total_TPM), y = Total_TPM)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  coord_flip() +
  labs(
    title = "Top 15 ARGs by Total TPM (All Samples)",
    x = "ARG Gene",
    y = "Total TPM"
  ) +
  theme_minimal(base_size = 16)

# Step 5: Summarize TPM per gene per sample
sample_summary <- filtered_data %>%
  group_by(directory, gene) %>%
  summarise(TPM = sum(TPM_Normalization, na.rm = TRUE)) %>%
  ungroup()

# Optional: Limit to top 10 most abundant genes for cleaner plotting
top_gene_list <- top_genes$gene  # Use genes from step 3
sample_filtered <- sample_summary %>%
  filter(gene %in% top_gene_list)

# Step 6: Plot ARG abundance per gene by sample
ggplot(sample_filtered, aes(x = gene, y = TPM, fill = directory)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Top ARGs per Sample (TPM Normalized)",
    x = "ARG Gene",
    y = "TPM"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
# Display the plot
print(ggplot)

# Save the plot as an image
ggsave("ARGs_SCGs_plot.png", ggplot, width = 10, height = 6)

#Statistical Analysis of Metagenomics- Shotgun Samples
# Load required libraries
library(tidyverse)
library(vegan)
library(readxl)
library(ggpubr)
library(reshape2)
library(readxl)
# 1. Read Data
card_data <- read_excel("Card_results_Updated2.xlsx", sheet = "card_results")

# 2. Clean & Prepare Data
card_data$Sample_ID <- card_data$Sample
card_data$Sample_Date <- as.Date(card_data$Sample_Date)
card_data <- card_data %>%
  mutate(
    TPM_Normalization = as.numeric(TPM_Normalization),
    FPKM_Normalization = as.numeric(FPKM_Normalization),
    `16S_Normalization` = as.numeric(`16S_Normalization`),
    RpoB_Normalization = as.numeric(RpoB_Normalization),
    SCGS_Normalization = as.numeric(SCGS_Normalization)
  )

# 3. Summarize ARG abundance by Sample and Resistance Type
normalize_and_summarize <- function(df, norm_col, label) {
  df %>% 
    filter(!is.na(.data[[norm_col]])) %>%
    group_by(Sample_ID, Sample_Day, Sample_Type, Resistance = type) %>% 
    summarise(Normalized_Abundance = sum(as.numeric(.data[[norm_col]]), na.rm = TRUE), .groups = 'drop') %>%
    mutate(Normalization = label)
}

scgs_data <- normalize_and_summarize(card_data, "SCGS_Normalization", "SCGS")

# 4. Boxplot of ARG Abundance by Sample Day (SCGS only)
p1 <- ggplot(scgs_data, aes(x = factor(Sample_Day), y = Normalized_Abundance + 1e-5, fill = Sample_Type)) +
  geom_boxplot() +
  scale_y_log10(
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    #breaks = scales::trans_breaks("log10", function(x) 10^x)
    breaks = 10^seq(-4, 0, by = 1),
  ) +
  labs(
    title = "SCGS-Normalized ARG Abundance by Sample Day",
    y = expression("Normalized Abundance ("*10^-x*")"),
    x = "Sample Day"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),    # remove major grid lines
    panel.grid.minor = element_blank(),    # remove minor grid lines
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
print(p1)
ggsave("SCGS.png", plot = p1, 
       width = 16, height = 10, dpi = 300)
#########################################################################
#SecondCheck of the SCGS-Normalized abundance


# 5. ARGs/SCGs Ratio Calculation (for visualization)
arg_scg_ratio <- card_data %>% 
  group_by(Sample_ID, Sample_Type) %>% 
  summarise(ArgSum_SCGS = sum(SCGS_Normalization, na.rm = TRUE),
            ArgSum_TPM = sum(TPM_Normalization, na.rm = TRUE), .groups = 'drop')

# 6. Kruskal-Wallis Tests
kw_scgs <- kruskal.test(ArgSum_SCGS ~ Sample_Type, data = arg_scg_ratio)
kw_tpm  <- kruskal.test(ArgSum_TPM ~ Sample_Type, data = arg_scg_ratio)

print(kw_scgs)
print(kw_tpm)
##########################################################################
#PCOA Plot
library(tidyverse)
library(vegan)
library(readxl)
library(ape)

# 1) Read and clean
df <- read_excel("Card_results_Updated2.xlsx") %>%
  mutate(
    Sample       = as.character(Sample),
    Sample_Type  = as.factor(Sample_Type),
    SCGS         = as.numeric(SCGS_Normalization),
    gene         = as.character(gene)
  )

# 2) Build samples × genes abundance matrix
abund <- df %>%
  filter(!is.na(SCGS)) %>%
  group_by(Sample, gene) %>%
  summarise(SCGS = sum(SCGS), .groups="drop") %>%
  pivot_wider(names_from = gene, values_from = SCGS, values_fill = 0)

# 3) Compute Bray–Curtis distance
mat <- abund %>% column_to_rownames("Sample")
bc_dist <- vegdist(mat, method = "bray")

# 4) Extract metadata
meta <- tibble(Sample = rownames(mat)) %>%
  left_join(df %>% select(Sample, Sample_Type) %>% distinct(), by = "Sample")

# 5) Run PCoA
pcoa_res <- pcoa(bc_dist)
pcoa_scores <- as.data.frame(pcoa_res$vectors[, 1:2]) %>%
  rownames_to_column("Sample") %>%
  left_join(meta, by = "Sample")

# 6) Filter to six desired Sample_Types
desired_groups <- c("Influent", "Biofilm-Control", "Bulkwater-Control",
                    "Biofilm-Spiked", "Bulkwater-Spiked")

pcoa_scores_filtered <- pcoa_scores %>%
  filter(Sample_Type %in% desired_groups)

pcoa_scores_filtered$Sample_Type <- factor(
  pcoa_scores_filtered$Sample_Type,
  levels = desired_groups
)

# 7) Plot with solid filled shapes
p2 <- ggplot(pcoa_scores_filtered, aes(
  x = Axis.1, y = Axis.2,
  fill = Sample_Type,
  shape = Sample_Type
)) +
  geom_point(size = 6, color = "black") +  # black border for contrast
  scale_shape_manual(values = c(21, 22, 23, 24, 25)) +  # 5 solid shapes
  theme_minimal(base_size = 14) +
  labs(
    title = "PCoA – SCGS (Bray–Curtis)",
    x = paste0("PCoA1 (", round(100 * pcoa_res$values$Relative_eig[1], 1), "%)"),
    y = paste0("PCoA2 (", round(100 * pcoa_res$values$Relative_eig[2], 1), "%)")
  ) +
  theme(
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.line         = element_line(color = "black"),
    axis.ticks        = element_line(color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    axis.title        = element_text(size = 32),
    axis.text         = element_text(size = 32),
    legend.title      = element_text(size = 28),
    legend.text       = element_text(size = 26)
  )

print(p2)
# Save PCoA plot to file
# width and height in inches, dpi in pixels per inch
ggsave("PCoA_SCGS_BrayCurtis.png", plot = p2, width = 9, height = 6, dpi = 300)
ggsave("PCoA_SCGS_BrayCurtis.pdf", plot = p2, width = 9, height = 6, dpi = 300)
# 4b) NMDS
nmds_res <- metaMDS(mat, distance="bray", k=2, trymax=100)
nmds_scores <- as.data.frame(scores(nmds_res)) %>%
  rownames_to_column("Sample") %>%
  left_join(meta, by="Sample")

# plot NMDS
ggplot(nmds_scores, aes(x=NMDS1, y=NMDS2, color=Sample_Type)) +
  geom_point(size=3) +
  theme_minimal() +
  labs(
    title    = "NMDS – SCGS (Bray–Curtis)",
    subtitle = paste("Stress =", round(nmds_res$stress, 4))
  )
##########################################################################
# ==============================
# Load required libraries
# ==============================
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
# 1. Read Data
# ==============================
card_data_raw <- read_excel(
  "Card_results_latest.xlsx",
  sheet = "card_results"
)
# 2. Clean & Prepare Data
# ==============================
card_data <- card_data_raw %>%
  mutate(
    Sample_ID = as.character(Sample),
    gene = as.character(gene),
    SCGS_Normalization = as.numeric(SCGS_Normalization),
    Sample_Type = as.character(Sample_Type),
    Treatment = sub(".*-", "", Sample_Type),
    Matrix_Type = ifelse(
      grepl("^Influent", Sample_Type),
      "Bulkwater",
      sub("-.*", "", Sample_Type)
    )
  ) %>%
  separate(type, into = c("type", "genetype"), sep = "_")

card_data$type <- toupper(card_data$type)
# 3. Import manual resistance classification
# ==============================
FRSheet <- read_excel("CARD - Manual Classification.xlsx")
FRSheet$Resistance_Genes <- toupper(FRSheet$Resistance_Genes)
# 4. Merge with resistance categories
# ==============================
Alldata_MergedResistances <- merge.data.frame(
  card_data,
  FRSheet,
  by.x = "type",
  by.y = "Resistance_Genes",
  sort = TRUE,
  all.x = TRUE
)
# 5. Get top 15 genes by mean abundance
# ==============================
top_genes <- card_data %>%
  group_by(type) %>%
  summarise(
    mean_abundance = mean(SCGS_Normalization, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abundance)) %>%
  slice_head(n = 15) %>%
  pull(type)
# 6. Filter to top 15 genes
# ==============================
filtered_data <- card_data %>%
  filter(type %in% top_genes)
# 7. Aggregate data
# ==============================
aggregated_data <- filtered_data %>%
  group_by(Sample_ID, Sample_Week, type, Treatment, Matrix_Type) %>%
  summarise(
    Total_Abundance = sum(SCGS_Normalization, na.rm = TRUE),
    .groups = "drop"
  )
# 8. Define a deeper and more distinct color palette
# ==============================
color_palette <- c(
  "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
  "#E6AB02", "#A6761D", "#666666", "#1F78B4", "#B2DF8A",
  "#FB9A99", "#FDBF6F", "#CAB2D6", "#FF7F00", "#6A3D9A"
)
# 9. Set factor levels for ARGs and Treatment
# ==============================
gene_levels <- rev(top_genes)  # most abundant at bottom

aggregated_data$type <- factor(
  aggregated_data$type,
  levels = gene_levels
)

aggregated_data$Treatment <- factor(
  aggregated_data$Treatment,
  levels = c("Influent", "Control", "Spiked")
)
# 10. Create the plot
# ==============================
p <- ggplot(
  aggregated_data,
  aes(x = Sample_Week, y = Total_Abundance, fill = type)
) +
  geom_bar(stat = "identity") +
  labs(
    title = "Top 15 ARGs",
    x = "Sample Week",
    y = "Abundance",
    fill = "ARG"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    strip.background = element_rect(fill = "lightgray"),
    strip.placement = "outside",
    strip.text.y.right = element_text(face = "bold", size = 13),
    strip.text.x = element_text(face = "bold", size = 13),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    )
  ) +
  facet_grid(
    rows = vars(Matrix_Type),
    cols = vars(Treatment),
    scales = "free_x",
    space = "free",
    switch = "y"
  ) +
  scale_fill_manual(values = color_palette)
# 11. Display the plot
# ==============================
print(p)
# ==============================
# 12. Save plot
ggsave("Top15_ARGs_Stacked_Bar_Plot.png", plot = p, width = 12, height = 8, dpi = 300)
ggsave("Top15_ARGs_Stacked_Bar_Plot.pdf", plot = p, width = 12, height = 8, dpi = 300)
###############################################################################
###############################################################################
# Load required libraries
#Relative abundance across samples and matrix. ARG abundance list from bottom to top
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

# 1. Read Data
card_data_raw <-- read_excel("Card_results_latest.xlsx", sheet = "card_results")

# 2. Clean & Prepare Data
card_data <- card_data_raw %>%
  mutate(
    Sample_ID = as.character(Sample),
    gene = as.character(gene),
    SCGS_Normalization = as.numeric(SCGS_Normalization),
    Sample_Type = as.character(Sample_Type),
    Treatment = sub(".*-", "", Sample_Type),
    Matrix_Type = ifelse(grepl("^Influent", Sample_Type), "Bulkwater", sub("-.*", "", Sample_Type)),
    Sample_Week = as.character(Sample_Week)
  ) %>%
  separate(type, into = c("type", "genetype"), sep = "_")

card_data$type <- toupper(card_data$type)

# ---- FIX 1: Force Sample_Week ordering (Week1...Week10...) ----
card_data <- card_data %>%
  mutate(
    Week_num = suppressWarnings(as.integer(gsub("[^0-9]", "", Sample_Week)))
  )

week_levels <- card_data %>%
  filter(!is.na(Week_num)) %>%
  distinct(Week_num) %>%
  arrange(Week_num) %>%
  pull(Week_num)

card_data <- card_data %>%
  mutate(
    Sample_Week = factor(Sample_Week, levels = paste0("Week", week_levels))
  )

# 3. Import manual resistance classification (kept as you had it)
FRSheet <- read_excel("CARD - Manual Classification.xlsx")
FRSheet$Resistance_Genes <- toupper(FRSheet$Resistance_Genes)

Alldata_MergedResistances <- merge.data.frame(
  card_data,
  FRSheet,
  by.x = "type",
  by.y = "Resistance_Genes",
  sort = TRUE,
  all.x = TRUE
)

# 4. Get top 15 genes by mean abundance (on full dataset)
top_genes <- card_data %>%
  group_by(type) %>%
  summarise(mean_abundance = mean(SCGS_Normalization, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abundance)) %>%
  slice_head(n = 15) %>%
  pull(type)

# 5. Filter to top 15 genes
filtered_data <- card_data %>%
  filter(type %in% top_genes)

# 6A. Per-sample aggregation first (ARG totals per sample/week)
per_sample <- filtered_data %>%
  group_by(Sample_ID, Sample_Week, type, Treatment, Matrix_Type) %>%
  summarise(Total_Abundance = sum(SCGS_Normalization, na.rm = TRUE), .groups = "drop")

# 6B. Aggregate across Sample_ID to get ONE value per Week × Treatment × Matrix × ARG
# Using MEAN across replicates
aggregated_data <- per_sample %>%
  group_by(Sample_Week, Treatment, Matrix_Type, type) %>%
  summarise(Mean_Abundance = mean(Total_Abundance, na.rm = TRUE), .groups = "drop")

# ---- CRITICAL FIX: make missing genes explicit zeros so every bar stacks to 1 ----
aggregated_data <- aggregated_data %>%
  tidyr::complete(
    Sample_Week,
    Treatment,
    Matrix_Type,
    type = top_genes,
    fill = list(Mean_Abundance = 0)
  )

# ---- FIX 2: Normalize within each Week × Treatment × Matrix so each bar hits 1.0 ----
aggregated_data <- aggregated_data %>%
  group_by(Sample_Week, Treatment, Matrix_Type) %>%
  mutate(
    Total = sum(Mean_Abundance, na.rm = TRUE),
    Rel_Abundance = ifelse(Total > 0, Mean_Abundance / Total, 0)
  ) %>%
  ungroup() %>%
  select(-Total)

# 7. Colors
color_palette <- c(
  "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
  "#E6AB02", "#A6761D", "#666666", "#1F78B4", "#B2DF8A",
  "#FB9A99", "#FDBF6F", "#CAB2D6", "#FF7F00", "#6A3D9A"
)

# 8. Factor levels
gene_levels <- rev(top_genes)
aggregated_data$type <- factor(aggregated_data$type, levels = gene_levels)
aggregated_data$Treatment <- factor(aggregated_data$Treatment, levels = c("Influent", "Control", "Spiked"))

# 9. Optional QC: verify each bar sums to 1 (within rounding)
qc <- aggregated_data %>%
  group_by(Sample_Week, Treatment, Matrix_Type) %>%
  summarise(sum_rel = sum(Rel_Abundance, na.rm = TRUE), .groups = "drop")

print(qc)

# 10. Plot (ONE bar per week per facet, always reaches 1 when Total>0)
p <- ggplot(aggregated_data, aes(x = Sample_Week, y = Rel_Abundance, fill = type)) +
  geom_col(width = 0.9) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(
    title = "Top 15 ARGs (Relative Abundance)",
    x = "Sample Week",
    y = "Relative abundance (0–1)",
    fill = "ARG"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    strip.background = element_rect(fill = "lightgray"),
    strip.placement = "outside",
    strip.text.y.right = element_text(face = "bold", size = 13),
    strip.text.x = element_text(face = "bold", size = 13),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  ) +
  facet_grid(rows = vars(Matrix_Type), cols = vars(Treatment),
             scales = "free_x", space = "free", switch = "y") +
  scale_fill_manual(values = color_palette)

print(p)

ggsave("Top15_ARGs_RelativeAbundance_OneBarPerWeek.png", plot = p, width = 12, height = 6, dpi = 300)
ggsave("Top15_ARGs_RelativeAbundance_OneBarPerWeek.pdf", plot = p, width = 12, height = 6, dpi = 300)
############################################################################################################
###########################################################################################################
# Load required libraries
# Relative abundance across samples and matrix, excluding the influent. ARG abundance list from bottom to top
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

# 1. Read Data
card_data_raw <- read_excel("Card_results_latest.xlsx", sheet = "card_results")

# 2. Clean & Prepare Data
card_data <- card_data_raw %>%
  mutate(
    Sample_ID = as.character(Sample),
    gene = as.character(gene),
    SCGS_Normalization = as.numeric(SCGS_Normalization),
    Sample_Type = as.character(Sample_Type),
    Treatment = sub(".*-", "", Sample_Type),
    Matrix_Type = ifelse(grepl("^Influent", Sample_Type), "Bulk Wwater", sub("-.*", "", Sample_Type)),
    Sample_Week = as.character(Sample_Week)
  ) %>%
  separate(type, into = c("type", "genetype"), sep = "_")

card_data$type <- toupper(card_data$type)

# ---- FIX 1: Force Sample_Week ordering (Week1...Week10...) ----
card_data <- card_data %>%
  mutate(
    Week_num = suppressWarnings(as.integer(gsub("[^0-9]", "", Sample_Week)))
  )

# ---- Keep only Week 1 to Week 7 for analysis ----
card_data <- card_data %>%
  filter(!is.na(Week_num) & Week_num >= 1 & Week_num <= 7)

week_levels <- card_data %>%
  distinct(Week_num) %>%
  arrange(Week_num) %>%
  pull(Week_num)

card_data <- card_data %>%
  mutate(
    Sample_Week = factor(Sample_Week, levels = paste0("Week", week_levels))
  )

# ---- EXCLUDE INFLUENT SAMPLES ----
card_data <- card_data %>%
  filter(Treatment != "Influent")

# 3. Import manual resistance classification (kept as you had it)
FRSheet <- read_excel("CARD - Manual Classification.xlsx")
FRSheet$Resistance_Genes <- toupper(FRSheet$Resistance_Genes)

Alldata_MergedResistances <- merge.data.frame(
  card_data,
  FRSheet,
  by.x = "type",
  by.y = "Resistance_Genes",
  sort = TRUE,
  all.x = TRUE
)

# 4. Get top 15 genes by mean abundance (NOW only using Week 1-7 and excluding Influent)
top_genes <- card_data %>%
  group_by(type) %>%
  summarise(mean_abundance = mean(SCGS_Normalization, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abundance)) %>%
  slice_head(n = 15) %>%
  pull(type)

# 5. Filter to top 15 genes
filtered_data <- card_data %>%
  filter(type %in% top_genes)

# 6A. Per-sample aggregation first (ARG totals per sample/week)
per_sample <- filtered_data %>%
  group_by(Sample_ID, Sample_Week, type, Treatment, Matrix_Type) %>%
  summarise(Total_Abundance = sum(SCGS_Normalization, na.rm = TRUE), .groups = "drop")

# 6B. Aggregate across Sample_ID to get ONE value per Week × Treatment × Matrix × ARG
# Using MEAN across replicates
aggregated_data <- per_sample %>%
  group_by(Sample_Week, Treatment, Matrix_Type, type) %>%
  summarise(Mean_Abundance = mean(Total_Abundance, na.rm = TRUE), .groups = "drop")

# ---- CRITICAL FIX: make missing genes explicit zeros so every bar stacks to 1 ----
aggregated_data <- aggregated_data %>%
  tidyr::complete(
    Sample_Week,
    Treatment,
    Matrix_Type,
    type = top_genes,
    fill = list(Mean_Abundance = 0)
  )

# ---- Normalize within each Week × Treatment × Matrix so each bar hits 1.0 ----
aggregated_data <- aggregated_data %>%
  group_by(Sample_Week, Treatment, Matrix_Type) %>%
  mutate(
    Total = sum(Mean_Abundance, na.rm = TRUE),
    Rel_Abundance = ifelse(Total > 0, Mean_Abundance / Total, 0)
  ) %>%
  ungroup() %>%
  select(-Total)

# 7. Colors
color_palette <- c(
  "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
  "#E6AB02", "#A6761D", "#666666", "#1F78B4", "#B2DF8A",
  "#FB9A99", "#FDBF6F", "#CAB2D6", "#FF7F00", "#6A3D9A"
)

# 8. Factor levels
gene_levels <- rev(top_genes)
aggregated_data$type <- factor(aggregated_data$type, levels = gene_levels)
aggregated_data$Treatment <- factor(aggregated_data$Treatment, levels = c("Control", "Spiked"))

# 9. Optional QC: verify each bar sums to 1 (within rounding)
qc <- aggregated_data %>%
  group_by(Sample_Week, Treatment, Matrix_Type) %>%
  summarise(sum_rel = sum(Rel_Abundance, na.rm = TRUE), .groups = "drop")

print(qc)

# 10. Plot (ONE bar per week per facet, always reaches 1 when Total>0)
p <- ggplot(aggregated_data, aes(x = Sample_Week, y = Rel_Abundance, fill = type)) +
  geom_col(width = 0.9) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(
    title = "Top 15 ARGs (Relative Abundance) – Weeks 1–7 (Control vs Spiked)",
    x = "Sample Week",
    y = "Relative Abundance",
    fill = "ARG"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    strip.background = element_rect(fill = "lightgray"),
    strip.placement = "outside",
    strip.text.y.right = element_text(face = "bold", size = 13),
    strip.text.x = element_text(face = "bold", size = 13),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.spacing.y = unit(0.7, "lines"),
    plot.margin = margin(t = 10, r = 10, b = 15, l = 10)
  ) +
  facet_grid(
    rows = vars(Matrix_Type),
    cols = vars(Treatment),
    scales = "free_x",
    space = "free",
    switch = "y"
  ) +
  scale_fill_manual(values = color_palette)

print(p)

ggsave("Top15_ARGs_RelativeAbundance_Weeks1to7_ControlVsSpiked.png", plot = p, width = 10, height = 6.5, dpi = 300)
ggsave("Top15_ARGs_RelativeAbundance_Weeks1to7_ControlVsSpiked.pdf", plot = p, width = 10, height = 6.5, dpi = 300)
############################################################################################################
###########################################################################################################
# Load required libraries
# Absolute abundance (SCG-normalized) across samples and matrix, INCLUDING influent.
# One bar per Week × Treatment × Matrix, stacked by top 15 ARGs (bottom-to-top by abundance order).
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

# 1. Read Data
card_data_raw <- read_excel("Card_results_latest.xlsx", sheet = "card_results")

# 2. Clean & Prepare Data
card_data <- card_data_raw %>%
  mutate(
    Sample_ID = as.character(Sample),
    gene = as.character(gene),
    SCGS_Normalization = as.numeric(SCGS_Normalization),
    Sample_Type = as.character(Sample_Type),
    Treatment = sub(".*-", "", Sample_Type),
    Matrix_Type = ifelse(grepl("^Influent", Sample_Type), "Bulk Wwater", sub("-.*", "", Sample_Type)),
    Sample_Week = as.character(Sample_Week)
  ) %>%
  separate(type, into = c("type", "genetype"), sep = "_")

card_data$type <- toupper(card_data$type)

# ---- FIX 1: Force Sample_Week ordering (Week1...Week10...) ----
card_data <- card_data %>%
  mutate(
    Week_num = suppressWarnings(as.integer(gsub("[^0-9]", "", Sample_Week)))
  )

# ---- Keep only Week 1 to Week 7 for analysis ----
card_data <- card_data %>%
  filter(!is.na(Week_num) & Week_num >= 1 & Week_num <= 7)

week_levels <- card_data %>%
  distinct(Week_num) %>%
  arrange(Week_num) %>%
  pull(Week_num)

card_data <- card_data %>%
  mutate(
    Sample_Week = factor(Sample_Week, levels = paste0("Week", week_levels))
  )

# 3. Import manual resistance classification (kept as you had it)
FRSheet <- read_excel("CARD - Manual Classification.xlsx")
FRSheet$Resistance_Genes <- toupper(FRSheet$Resistance_Genes)

Alldata_MergedResistances <- merge.data.frame(
  card_data,
  FRSheet,
  by.x = "type",
  by.y = "Resistance_Genes",
  sort = TRUE,
  all.x = TRUE
)

# 4. Get top 15 genes by mean abundance (Weeks 1-7, INCLUDING influent)
top_genes <- card_data %>%
  group_by(type) %>%
  summarise(mean_abundance = mean(SCGS_Normalization, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abundance)) %>%
  slice_head(n = 15) %>%
  pull(type)

# 5. Filter to top 15 genes
filtered_data <- card_data %>%
  filter(type %in% top_genes)

# 6A. Per-sample aggregation first (ARG totals per sample/week)
per_sample <- filtered_data %>%
  group_by(Sample_ID, Sample_Week, type, Treatment, Matrix_Type) %>%
  summarise(Total_Abundance = sum(SCGS_Normalization, na.rm = TRUE), .groups = "drop")

# 6B. Aggregate across Sample_ID to get ONE value per Week × Treatment × Matrix × ARG
# Using MEAN across replicates (SCG-normalized absolute abundance)
aggregated_data <- per_sample %>%
  group_by(Sample_Week, Treatment, Matrix_Type, type) %>%
  summarise(Mean_Abundance = mean(Total_Abundance, na.rm = TRUE), .groups = "drop")

# ---- Make missing genes explicit zeros (so stacking is consistent week-to-week) ----
aggregated_data <- aggregated_data %>%
  tidyr::complete(
    Sample_Week,
    Treatment,
    Matrix_Type,
    type = top_genes,
    fill = list(Mean_Abundance = 0)
  )

# 7. Colors
color_palette <- c(
  "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
  "#E6AB02", "#A6761D", "#666666", "#1F78B4", "#B2DF8A",
  "#FB9A99", "#FDBF6F", "#CAB2D6", "#FF7F00", "#6A3D9A"
)

# 8. Factor levels
gene_levels <- rev(top_genes)
aggregated_data$type <- factor(aggregated_data$type, levels = gene_levels)
aggregated_data$Treatment <- factor(aggregated_data$Treatment, levels = c("Influent", "Control", "Spiked"))

# 9. Optional QC: check totals per bar (absolute SCG-normalized totals across top 15)
qc_abs <- aggregated_data %>%
  group_by(Sample_Week, Treatment, Matrix_Type) %>%
  summarise(total_top15 = sum(Mean_Abundance, na.rm = TRUE), .groups = "drop")

print(qc_abs)

# 10. Plot (Absolute SCG-normalized abundance; bars can be <1, >1, etc.)
p <- ggplot(aggregated_data, aes(x = Sample_Week, y = Mean_Abundance, fill = type)) +
  geom_col(width = 0.9) +
  labs(
    title = "Top 15 ARGs (SCG-normalized Absolute Abundance) – Weeks 1–7",
    x = "Sample Week",
    y = "SCG-normalized abundance",
    fill = "ARG"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    strip.background = element_rect(fill = "lightgray"),
    strip.placement = "outside",
    strip.text.y.right = element_text(face = "bold", size = 13),
    strip.text.x = element_text(face = "bold", size = 13),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.spacing.y = unit(0.7, "lines"),
    plot.margin = margin(t = 10, r = 10, b = 15, l = 10)
  ) +
  facet_grid(
    rows = vars(Matrix_Type),
    cols = vars(Treatment),
    scales = "free_x",
    space = "free",
    switch = "y"
  ) +
  scale_fill_manual(values = color_palette)

print(p)

ggsave("Top15_ARGs_AbsoluteAbundance_Weeks1to7.png", plot = p, width = 10, height = 6.5, dpi = 300)
ggsave("Top15_ARGs_AbsoluteAbundance_Weeks1to7.pdf", plot = p, width = 10, height = 6.5, dpi = 300)
