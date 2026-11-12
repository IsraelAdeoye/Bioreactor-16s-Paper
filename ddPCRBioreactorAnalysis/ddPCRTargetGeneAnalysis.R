###############################################################################
###############################################################################
###############################################################################
# Influent: geometric mean, geometric SD error bars, and LOQ lines
library(ggplot2)
library(dplyr)
library(readxl)

# 1. LOAD AND PREPARE DATA
ddpcr_1 <- read_excel("Israel_YFP_ARG_ddPCR.xlsx")
ddpcr_2 <- read_excel("Israel_YFP_ARG_ddPCR.xlsx", sheet = "Combined_VanATetA")

names(ddpcr_1)[names(ddpcr_1) == "Conc(copies/mL)"] <- "Concentration"
names(ddpcr_2)[names(ddpcr_2) == "Conc(copies/mL)"] <- "Concentration"

ddpcr_all <- bind_rows(
  ddpcr_1 %>% select(Description, SampleWeek, Probe, Concentration),
  ddpcr_2 %>% select(Description, SampleWeek, Probe, Concentration) %>% 
    mutate(Concentration = as.numeric(Concentration) * 10)
)

cols <- c("NDM" = "steelblue", "YFP" = "firebrick", "vanA" = "darkgreen", "tetA" = "purple")
shps <- c("NDM" = 16, "YFP" = 17, "vanA" = 15, "tetA" = 18)

# Define strict chronological order to prevent Week 10 from jumping to the front
weeks_order <- paste("Week", 1:10)

# 2. FILTER AND REFACTOR INFLUENT DATA
influent_data <- ddpcr_all %>% 
  filter(
    Description == "Influent", 
    SampleWeek %in% weeks_order, 
    Probe %in% names(cols)
  ) %>% 
  mutate(
    #Convert to factor here with explicit levels to maintain timeline order
    SampleWeek = factor(SampleWeek, levels = weeks_order, ordered = TRUE),
    Probe = factor(Probe, levels = names(cols)),
    Concentration = as.numeric(Concentration),
    Concentration_plot = ifelse(is.na(Concentration), NA, pmax(Concentration, 0.001))
  )

# 3. CALCULATE GEOMETRIC MEANS AND ERROR BOUNDARIES
summary_df <- influent_data %>% 
  group_by(Probe, SampleWeek) %>% 
  summarise(
    GeoMean = exp(mean(log(Concentration_plot), na.rm = TRUE)),
    SDlog = sd(log(Concentration_plot), na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  mutate(
    Lower = GeoMean / exp(SDlog),
    Upper = GeoMean * exp(SDlog)
  )

# 4. DEFINE LOQ LIMITS DATA FRAME
loq_df <- data.frame(
  Probe = factor(names(cols), levels = names(cols)), 
  LOQ = c(263, 420, 7.65, 79.9)
) %>%
  mutate(LOQ_plot = ifelse(LOQ <= 0, 0.001, LOQ))

# 5. GENERATE TEMPORAL TREND PLOT
p1 <- ggplot(summary_df, aes(x = SampleWeek, y = GeoMean, color = Probe, shape = Probe, group = Probe)) +
  geom_line(aes(linetype = "Influent"), linewidth = 1.6) +
  geom_point(size = 7.0) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.25, linewidth = 1.2) +
  
  # Dynamic LOQ lines mapping
  geom_hline(
    data = loq_df, 
    aes(yintercept = LOQ_plot, color = Probe, linetype = "LOQ"), 
    linewidth = 1.2
  ) +
  
  # Manual Scale Aesthetics mappings
  scale_color_manual(
    name = "Target Gene", 
    values = cols,
    labels = c("NDM" = "NDM", "YFP" = "YFP", "vanA" = expression(italic("vanA")), "tetA" = expression(italic("tetA")))
  ) + 
  scale_shape_manual(
    name = "Target Gene", 
    values = shps,
    labels = c("NDM" = "NDM", "YFP" = "YFP", "vanA" = expression(italic("vanA")), "tetA" = expression(italic("tetA")))
  ) + 
  scale_linetype_manual(
    "Line Type", 
    # CRITICAL FIX: Changed "LOD" key to "LOQ" to match your data aesthetic mapping safely
    values = c("Influent" = "solid", "LOQ" = "dotdash")
  ) +
  scale_y_log10(
    limits = c(1e-5, 1e7),
    breaks = 10^seq(-5, 7, by = 2),
    labels = c(
      expression(10^-5),
      expression(10^-3),
      expression(10^-1),
      expression(10^1),
      expression(10^3),
      expression(10^5),
      expression(10^7)
    )
  ) +
  labs(
    #title = "Temporal trends of target gene concentrations in influent samples",
    x = "Sample Week",
    y = "Concentration (copies/mL)"
  ) +
  # Controls the exact vertical stack order of your legends
  guides(
    linetype = guide_legend(order = 1), # Forces Line Type to the top
    color = guide_legend(order = 2),    # Forces Target Gene below it
    shape = guide_legend(order = 2)     # Keeps shapes synced with colors
  ) +
  theme_classic(base_size = 28) + 
  theme(
    plot.title = element_text(face = "bold", size = 26, hjust = 0.5),
    plot.subtitle = element_text(size = 26, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 26),
    axis.text.x = element_text(size = 26, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 26),
    legend.title = element_text(face = "bold", size = 26),
    legend.text = element_text(size = 26),
    legend.position = "right",
    legend.box = "vertical", # Ensures the two legend groups stack cleanly on top of each other
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7)
  )

print(p1)

# 6. SAVE PLOTS (PNG & PDF)
ggsave("Influent_All_Targets_GeoMean_LOQ_errorbars.png", p1, width = 14, height = 8, dpi = 600)
ggsave("Influent_All_Targets_GeoMean_LOQ_errorbars.pdf", p1, width = 14, height = 8)
###############################################################################
###############################################################################
###############################################################################
#Week 1-7: # Bulkwater-Control/Bulkwater-Spiked: all target genes line plot with geometric mean and LOD legend

library(ggplot2)
library(dplyr)
library(readxl)
library(scales)

ddpcr_1 <- read_excel("Israel_YFP_ARG_ddPCR.xlsx")
names(ddpcr_1)[names(ddpcr_1) == "Conc(copies/mL)"] <- "Concentration"

ddpcr_2 <- read_excel("Israel_YFP_ARG_ddPCR.xlsx", sheet = "Combined_VanATetA")
names(ddpcr_2)[names(ddpcr_2) == "Conc(copies/mL)"] <- "Concentration"

ddpcr_1 <- ddpcr_1 %>% select(Description, SampleWeek, Probe, Concentration)
ddpcr_2 <- ddpcr_2 %>% select(Description, SampleWeek, Probe, Concentration) %>%
  mutate(Concentration = as.numeric(Concentration) * 10)

ddpcr_all <- bind_rows(ddpcr_1, ddpcr_2)

cols <- c("NDM" = "steelblue", "YFP" = "firebrick", "vanA" = "darkgreen", "tetA" = "purple")
shps <- c("NDM" = 16, "YFP" = 17, "vanA" = 15, "tetA" = 18)

bulkwater_data <- ddpcr_all %>%
  filter(
    Description %in% c("Bulkwater-Control", "Bulkwater-Spiked"),
    SampleWeek %in% paste("Week", 1:7),
    Probe %in% c("NDM", "YFP", "vanA", "tetA")
  ) %>%
  mutate(
    Description = factor(Description, levels = c("Bulkwater-Control", "Bulkwater-Spiked")),
    SampleWeek = factor(SampleWeek, levels = paste("Week", 1:7)),
    Probe = factor(Probe, levels = c("NDM", "YFP", "vanA", "tetA")),
    Concentration = as.numeric(Concentration),
    Concentration_plot = ifelse(is.na(Concentration), NA,
                                ifelse(Concentration <= 0, 0.001, Concentration))
  )

summary_df <- bulkwater_data %>%
  group_by(Description, Probe, SampleWeek) %>%
  summarise(
    GeoMean = exp(mean(log(Concentration_plot), na.rm = TRUE)),
    SDlog = sd(log(Concentration_plot), na.rm = TRUE),
    N = sum(!is.na(Concentration_plot)),
    .groups = "drop"
  ) %>%
  mutate(
    Lower = GeoMean / exp(SDlog),
    Upper = GeoMean * exp(SDlog)
  )

lod_df <- data.frame(
  Probe = factor(c("NDM", "YFP", "vanA", "tetA"),
                 levels = c("NDM", "YFP", "vanA", "tetA")),
  LOD = c(
    17.51641,      # NDM
    15.40229885,   # YFP
    3.3353,       # vanA
    3.94        # tetA
  )
) %>%
  mutate(
    LOD_plot = ifelse(LOD <= 0, 0.001, LOD)
  )
p1 <- ggplot(
  summary_df,
  aes(
    x = SampleWeek,
    y = GeoMean,
    color = Probe,
    shape = Probe,
    group = interaction(Description, Probe),
    linetype = Description
  )
) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3.8) +
  geom_errorbar(
    aes(ymin = Lower, ymax = Upper),
    width = 0.15,
    linewidth = 0.8
  ) +
  geom_hline(
    data = lod_df,
    aes(
      yintercept = LOD_plot,
      color = Probe,
      linetype = "LOD"
    ),
    linewidth = 0.5,
    show.legend = TRUE
  ) +
  scale_linetype_manual(
    name = "Line Type",
    values = c(
      "Bulkwater-Control" = "dashed",
      "Bulkwater-Spiked" = "solid",
      "LOD" = "dotdash"
    )
  ) +
  scale_color_manual(name = "Target Gene", values = cols) +
  scale_shape_manual(name = "Target Gene", values = shps) +
  scale_y_log10(
    limits = c(1e-5, 1e7),
    breaks = 10^seq(-5, 7, by = 2),
    labels = c(
      expression(10^-5),
      expression(10^-3),
      expression(10^-1),
      expression(10^1),
      expression(10^3),
      expression(10^5),
      expression(10^7)
    )
  ) +
  labs(
    title = "Temporal trends of target gene concentrations in bulkwater samples",
    x = "Sample Week",
    y = "Concentration (copies/mL)"
  ) +
  scale_x_discrete(expand = expansion(add = c(0.4, 0.4))) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 15, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 14, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 14),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7)
  )

print(p1)

ggsave("Bulkwater_All_Targets_log_lineplot_errorbars.png", p1, width = 10, height = 7, dpi = 600)
#################################################################################
###################################################################################
##################################################################################
# Week 1-10: Bulkwater-Control/Bulkwater-Spiked: all target genes line plot with geometric mean and LOQ legend
library(ggplot2)
library(dplyr)
library(readxl)
library(scales)

# 1. LOAD AND PREPARE DATA
ddpcr_1 <- read_excel("Israel_YFP_ARG_ddPCR.xlsx")
names(ddpcr_1)[names(ddpcr_1) == "Conc(copies/mL)"] <- "Concentration"

ddpcr_2 <- read_excel("Israel_YFP_ARG_ddPCR.xlsx", sheet = "Combined_VanATetA")
names(ddpcr_2)[names(ddpcr_2) == "Conc(copies/mL)"] <- "Concentration"

ddpcr_1 <- ddpcr_1 %>% select(Description, SampleWeek, Probe, Concentration)
ddpcr_2 <- ddpcr_2 %>% select(Description, SampleWeek, Probe, Concentration) %>% 
  mutate(Concentration = as.numeric(Concentration) * 10)

ddpcr_all <- bind_rows(ddpcr_1, ddpcr_2)

cols <- c("NDM" = "steelblue", "YFP" = "firebrick", "vanA" = "darkgreen", "tetA" = "purple")
shps <- c("NDM" = 16, "YFP" = 17, "vanA" = 15, "tetA" = 18)

# 2. FILTER FOR WEEKS 1 TO 10
target_weeks <- paste("Week", 1:10)

bulkwater_data <- ddpcr_all %>% 
  filter(
    Description %in% c("Bulkwater-Control", "Bulkwater-Spiked"),
    SampleWeek %in% target_weeks,
    Probe %in% c("NDM", "YFP", "vanA", "tetA")
  ) %>% 
  mutate(
    Description = factor(Description, levels = c("Bulkwater-Control", "Bulkwater-Spiked")),
    SampleWeek = factor(SampleWeek, levels = target_weeks),
    Probe = factor(Probe, levels = c("NDM", "YFP", "vanA", "tetA")),
    Concentration = as.numeric(Concentration),
    Concentration_plot = ifelse(is.na(Concentration), NA, ifelse(Concentration <= 0, 0.001, Concentration))
  )

# 3. CALCULATE GEOMETRIC MEANS AND STANDARD DEVIATIONS
summary_df <- bulkwater_data %>% 
  group_by(Description, Probe, SampleWeek) %>% 
  summarise(
    GeoMean = exp(mean(log(Concentration_plot), na.rm = TRUE)),
    SDlog = sd(log(Concentration_plot), na.rm = TRUE),
    N = sum(!is.na(Concentration_plot)), 
    .groups = "drop"
  ) %>% 
  mutate(
    Lower = GeoMean / exp(SDlog),
    Upper = GeoMean * exp(SDlog)
  )

# 4. DEFINE NEW LOQ VALUES
loq_df <- data.frame(
  Probe = factor(c("NDM", "YFP", "vanA", "tetA"), levels = c("NDM", "YFP", "vanA", "tetA")),
  LOQ = c(
    145, # NDM
    231, # YFP
    4.22,   # vanA
    44.1   # tetA
  )
) %>% 
  mutate(LOQ_plot = ifelse(LOQ <= 0, 0.001, LOQ))

# 5. GENERATE PLOT
p1 <- ggplot(summary_df, aes(x = SampleWeek, y = GeoMean, color = Probe, shape = Probe, group = interaction(Description, Probe), linetype = Description)) +
  
  # CRITICAL LOG-SCALE FIX: Replaced -Inf/Inf with your exact axis view limits (1e-5 and 1e9)
  # This stops R from encountering log(0) calculation errors so the gray background displays properly.
  annotate(
    "rect",
    xmin = 7.8, xmax = 10.4, ymin = 1e-5, ymax = 1e9,
    fill = "grey92", alpha = 1
  ) +
  
  # Data trends layers
  geom_line(linewidth = 1.6) + 
  geom_point(size = 7.0) + 
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.25, linewidth = 1.2) + 
  
  # Horizontal LOQ lines
  geom_hline(
    data = loq_df, 
    aes(yintercept = LOQ_plot, color = Probe, linetype = "LOQ"), 
    linewidth = 1.2, show.legend = TRUE
  ) + 
  
  # Legend and scale configurations
  scale_linetype_manual(
    name = "Line Type", 
    values = c(
      "Bulkwater-Control" = "dashed", 
      "Bulkwater-Spiked" = "solid", 
      "LOQ" = "dotdash"
    ),
    labels = c(
      "Bulkwater-Control" = "Bulk Water - Control", 
      "Bulkwater-Spiked" = "Bulk Water - Spiked", 
      "LOQ" = "LOQ"
    )
  ) +
  scale_color_manual(
    name = "Target Gene", 
    values = cols,
    labels = c("NDM" = "NDM", "YFP" = "YFP", "vanA" = expression(italic("vanA")), "tetA" = expression(italic("tetA")))
  ) + 
  scale_shape_manual(
    name = "Target Gene", 
    values = shps,
    labels = c("NDM" = "NDM", "YFP" = "YFP", "vanA" = expression(italic("vanA")), "tetA" = expression(italic("tetA")))
  ) +  
  scale_y_log10(
    limits = c(1e-5, 1e9), 
    breaks = 10^seq(-5, 9, by = 2), 
    labels = c(
      expression(10^-5), expression(10^-3), expression(10^-1), 
      expression(10^1), expression(10^3), expression(10^5), 
      expression(10^7), expression(10^9)
    )
  ) + 
  labs(
    #title = "Temporal trends of target gene concentrations in bulkwater samples", 
    x = "Sample Week", 
    y = "Concentration (copies/mL)"
  ) + 
  scale_x_discrete(expand = expansion(add = c(0.4, 0.4))) + 
  
  # Vertical stacked legends order override setup
  guides(
    linetype = guide_legend(order = 1),
    color = guide_legend(order = 2),
    shape = guide_legend(order = 2)
  ) +
  
  theme_classic(base_size = 28) + 
  theme( 
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold", size = 26, hjust = 0.5), 
    plot.subtitle = element_text(size = 26, hjust = 0.5), 
    axis.title = element_text(face = "bold", size = 24), 
    axis.text.x = element_text(size = 26, angle = 30, hjust = 1), 
    axis.text.y = element_text(size = 26), 
    legend.title = element_text(face = "bold", size = 26), 
    legend.text = element_text(size = 26), 
    legend.position = "right", 
    legend.box = "vertical",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7) 
  )

print(p1)

# 6. SAVE OUTPUTS
ggsave("Bulkwater_All_Targets_log_lineplot_errorbars_Weeks1-10.png", p1, width = 14, height = 8, dpi = 600)
ggsave("Bulkwater_All_Targets_log_lineplot_errorbars_Weeks1-10.pdf", p1, width = 14, height = 8)
################
################
###############
#Week 8-10: # Bulkwater-Control/Bulkwater-Spiked: all target genes line plot with geometric mean and LOD legend

library(ggplot2)
library(dplyr)
library(readxl)
library(scales)

ddpcr_1 <- read_excel("Israel_YFP_ARG_ddPCR.xlsx")
names(ddpcr_1)[names(ddpcr_1) == "Conc(copies/mL)"] <- "Concentration"

ddpcr_2 <- read_excel("Israel_YFP_ARG_ddPCR.xlsx", sheet = "Combined_VanATetA")
names(ddpcr_2)[names(ddpcr_2) == "Conc(copies/mL)"] <- "Concentration"

ddpcr_1 <- ddpcr_1 %>% select(Description, SampleWeek, Probe, Concentration)
ddpcr_2 <- ddpcr_2 %>% select(Description, SampleWeek, Probe, Concentration) %>%
  mutate(Concentration = as.numeric(Concentration) * 10)

ddpcr_all <- bind_rows(ddpcr_1, ddpcr_2)

cols <- c("NDM" = "steelblue", "YFP" = "firebrick", "vanA" = "darkgreen", "tetA" = "purple")
shps <- c("NDM" = 16, "YFP" = 17, "vanA" = 15, "tetA" = 18)

bulkwater_data <- ddpcr_all %>%
  filter(
    Description %in% c("Bulkwater-Control", "Bulkwater-Spiked"),
    SampleWeek %in% paste("Week", 8:10),
    Probe %in% c("NDM", "YFP", "vanA", "tetA")
  ) %>%
  mutate(
    Description = factor(Description, levels = c("Bulkwater-Control", "Bulkwater-Spiked")),
    SampleWeek = factor(SampleWeek, levels = paste("Week", 8:10)),
    Probe = factor(Probe, levels = c("NDM", "YFP", "vanA", "tetA")),
    Concentration = as.numeric(Concentration),
    Concentration_plot = ifelse(is.na(Concentration), NA,
                                ifelse(Concentration <= 0, 0.001, Concentration))
  )

summary_df <- bulkwater_data %>%
  group_by(Description, Probe, SampleWeek) %>%
  summarise(
    GeoMean = exp(mean(log(Concentration_plot), na.rm = TRUE)),
    SDlog = sd(log(Concentration_plot), na.rm = TRUE),
    N = sum(!is.na(Concentration_plot)),
    .groups = "drop"
  ) %>%
  mutate(
    Lower = GeoMean / exp(SDlog),
    Upper = GeoMean * exp(SDlog)
  )

lod_df <- data.frame(
  Probe = factor(c("NDM", "YFP", "vanA", "tetA"),
                 levels = c("NDM", "YFP", "vanA", "tetA")),
  LOD = c(
    17.51641,      # NDM
    15.40229885,   # YFP
    3.3353,       # vanA
    3.94        # tetA
  )
) %>%
  mutate(
    LOD_plot = ifelse(LOD <= 0, 0.001, LOD)
  )
p1 <- ggplot(
  summary_df,
  aes(
    x = SampleWeek,
    y = GeoMean,
    color = Probe,
    shape = Probe,
    group = interaction(Description, Probe),
    linetype = Description
  )
) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3.8) +
  geom_errorbar(
    aes(ymin = Lower, ymax = Upper),
    width = 0.15,
    linewidth = 0.8
  ) +
  geom_hline(
    data = lod_df,
    aes(
      yintercept = LOD_plot,
      color = Probe,
      linetype = "LOD"
    ),
    linewidth = 0.5,
    show.legend = TRUE
  ) +
  scale_linetype_manual(
    name = "Line Type",
    values = c(
      "Bulkwater-Control" = "dashed",
      "Bulkwater-Spiked" = "solid",
      "LOD" = "dotdash"
    )
  ) +
  scale_color_manual(name = "Target Gene", values = cols) +
  scale_shape_manual(name = "Target Gene", values = shps) +
  scale_y_log10(
    limits = c(1e-3, 1e9),
    breaks = 10^seq(-3, 9, by = 2),
    labels = c(
      expression(10^-3),
      expression(10^-1),
      expression(10^1),
      expression(10^3),
      expression(10^5),
      expression(10^7),
      expression(10^9)
    )
  ) +
  labs(
    title = "Temporal trends in target gene concentrations after coupon swap",
    x = "Sample Week",
    y = "Concentration (copies/mL)"
  ) +
  scale_x_discrete(expand = expansion(add = c(0.15, 0.15))) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 15, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 14, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 14),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7)
  )

print(p1)

ggsave("Bulkwater_WK8-10_All_Targets_log_lineplot_errorbars.png", p1, width = 8, height = 7, dpi = 600)

###############################################################################
###############################################################################
###############################################################################
#Reactor Bulkwater- Added for reshaping data to calculate C/C0
library(ggplot2)
library(dplyr)
library(readxl)
library(scales)
library(tidyr) # Added for reshaping data to calculate C/C0

# 1. LOAD AND PREPARE DATA
ddpcr_1 <- read_excel("Israel_YFP_ARG_ddPCR.xlsx")
names(ddpcr_1)[names(ddpcr_1) == "Conc(copies/mL)"] <- "Concentration"

ddpcr_2 <- read_excel("Israel_YFP_ARG_ddPCR.xlsx", sheet = "Combined_VanATetA")
names(ddpcr_2)[names(ddpcr_2) == "Conc(copies/mL)"] <- "Concentration"

ddpcr_1 <- ddpcr_1 %>% select(Description, SampleWeek, Probe, Concentration)
ddpcr_2 <- ddpcr_2 %>% select(Description, SampleWeek, Probe, Concentration) %>% 
  mutate(Concentration = as.numeric(Concentration) * 10)

ddpcr_all <- bind_rows(ddpcr_1, ddpcr_2)

# Aesthetics Setup
cols <- c("NDM" = "steelblue", "YFP" = "firebrick", "vanA" = "darkgreen", "tetA" = "purple")
shps <- c("NDM" = 16, "YFP" = 17, "vanA" = 15, "tetA" = 18)

# 2. FILTER FOR WEEKS 1 TO 10
target_weeks <- paste("Week", 1:10)

bulkwater_data <- ddpcr_all %>%
  filter(
    Description %in% c("Bulkwater-Control", "Bulkwater-Spiked"),
    SampleWeek %in% target_weeks,
    Probe %in% c("NDM", "YFP", "vanA", "tetA")
  ) %>%
  mutate(
    Description = factor(Description, levels = c("Bulkwater-Control", "Bulkwater-Spiked")),
    SampleWeek = factor(SampleWeek, levels = target_weeks),
    Probe = factor(Probe, levels = c("NDM", "YFP", "vanA", "tetA")),
    Concentration = as.numeric(Concentration),
    Concentration_plot = ifelse(is.na(Concentration), NA, ifelse(Concentration <= 0, 0.001, Concentration))
  )

# 3. CALCULATE GEOMETRIC MEANS AND LOG-VARIANCES FOR ERROR PROPAGATION
summary_df <- bulkwater_data %>%
  group_by(Description, Probe, SampleWeek) %>%
  summarise(
    GeoMean = exp(mean(log(Concentration_plot), na.rm = TRUE)),
    VarLog = var(log(Concentration_plot), na.rm = TRUE),
    .groups = "drop"
  )

# 4. RESHAPE TO CALCULATE C/C0 RATIOS AND PROPAGATED ERROR BARS
ratio_df <- summary_df %>%
  pivot_wider(
    names_from = Description, 
    values_from = c(GeoMean, VarLog)
  ) %>%
  rename(
    C_0 = `GeoMean_Bulkwater-Control`, 
    C   = `GeoMean_Bulkwater-Spiked`
  ) %>%
  mutate(
    # C/C0 Ratio calculation
    C_C0 = C / C_0,
    
    # Propagated standard error for the ratio of two log-normal distributions
    SDlog_ratio = sqrt(replace_na(`VarLog_Bulkwater-Spiked`, 0) + replace_na(`VarLog_Bulkwater-Control`, 0)),
    Lower = C_C0 / exp(SDlog_ratio),
    Upper = C_C0 * exp(SDlog_ratio)
  )

# 5. DEFINE LOQ VALUES AND TRANSLATE TO DIMENSIONLESS LOQ RATIO (LOQ / C0)
loq_values <- data.frame(
  Probe = factor(c("NDM", "YFP", "vanA", "tetA"), levels = c("NDM", "YFP", "vanA", "tetA")),
  LOQ = c(145, 231, 4.22, 44.1)
)

plot_df <- ratio_df %>%
  left_join(loq_values, by = "Probe") %>%
  mutate(
    # The LOQ threshold on a C/C0 plot is dynamic: LOQ / C0
    LOQ_ratio = LOQ / C_0
  )

# 6. GENERIC PLOT GENERATION
p1 <- ggplot(plot_df, aes(x = SampleWeek, y = C_C0, color = Probe, shape = Probe, group = Probe)) +
  
  # SHADING LAYER: Background window for post coupon exchange phase (Weeks 8 to 10)
  annotate("rect", xmin = 7.8, xmax = 10.4, ymin = 1e-6, ymax = 1e8, fill = "grey92", alpha = 1) +
  # C/C0 Trends
  geom_line(linewidth = 1.1) +
  geom_point(size = 3.8) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.15, linewidth = 0.8) +
  
  # Dynamic LOQ threshold curves (LOQ / C0) per gene
  geom_line(aes(y = LOQ_ratio, linetype = "LOQ Threshold"), linewidth = 0.6) +
  
  # Manual Scalings
  scale_color_manual(
    name = "Target Gene", 
    values = cols,
    labels = c("NDM" = "NDM", "YFP" = "YFP", "vanA" = expression(italic("vanA")), "tetA" = expression(italic("tetA")))
  ) + 
  scale_shape_manual(
    name = "Target Gene", 
    values = shps,
    labels = c("NDM" = "NDM", "YFP" = "YFP", "vanA" = expression(italic("vanA")), "tetA" = expression(italic("tetA")))
  ) + 
  scale_linetype_manual(name = "Reference", values = c("LOQ Threshold" = "dotdash")) +
  
  # Log scale for ratios (Centered around 10^0 = 1 where C == C0)
  scale_y_log10(
    limits = c(1e-6, 1e8),
    breaks = 10^seq(-6, 8, by = 2),
    labels = c(
      expression(10^-6),
      expression(10^-4),
      expression(10^-2),
      expression(10^0),
      expression(10^2),
      expression(10^4),
      expression(10^6),
      expression(10^8)
    )
  ) +
  labs(
    x = "Sample Week",
    y = expression(C[reactor] == C[Spiked] / C[Control])
  ) +
  scale_x_discrete(expand = expansion(add = c(0.4, 0.4))) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5, face = "italic"),
    axis.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 14, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 14),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7)
  )

print(p1)
# =============================================================================
# 7. EXPORT DATA POINTS TO CSV
# =============================================================================
# Select and rename columns cleanly to mirror your variable mapping
creactor_table <- plot_df %>%
  select(
    Probe, 
    SampleWeek, 
    C_0_Control = C_0, 
    C_Spiked = C, 
    C_reactor = C_C0, 
    Lower_Error = Lower, 
    Upper_Error = Upper,
    LOQ_ratio
  ) %>%
  arrange(Probe, SampleWeek)

# Write out to a clean CSV file matching your formatting style
write.csv(creactor_table, "Reactor_Creactor_Data_Points_ByWeek.csv", row.names = FALSE)
# 8. SAVE OUTPUTS
ggsave("Bulkwater_All_Targets_C_C0_Ratio_Plot.png", p1, width = 11, height = 7, dpi = 600)
ggsave("Bulkwater_All_Targets_C_C0_Ratio_Plot.pdf", p1, width = 11, height = 7)
################################################################################
###############################################################################
##############################################################################
#WEEK 1-7: # Biofilm-Control/Biofilm-Spiked: all target genes line plot with geometric mean and LOD legend

library(ggplot2)
library(dplyr)
library(readxl)
library(scales)

ddpcr_1 <- read_excel("Israel_YFP_ARG_ddPCR.xlsx")
names(ddpcr_1)[names(ddpcr_1) == "Conc(copies/mL)"] <- "Concentration"

ddpcr_2 <- read_excel("Israel_YFP_ARG_ddPCR.xlsx", sheet = "Combined_VanATetA")
names(ddpcr_2)[names(ddpcr_2) == "Conc(copies/mL)"] <- "Concentration"

ddpcr_1 <- ddpcr_1 %>% select(Description, SampleWeek, Probe, Concentration)
ddpcr_2 <- ddpcr_2 %>% select(Description, SampleWeek, Probe, Concentration) %>%
  mutate(Concentration = as.numeric(Concentration) * 10)

ddpcr_all <- bind_rows(ddpcr_1, ddpcr_2)

cols <- c("NDM" = "steelblue", "YFP" = "firebrick", "vanA" = "darkgreen", "tetA" = "purple")
shps <- c("NDM" = 16, "YFP" = 17, "vanA" = 15, "tetA" = 18)

biofilm_data <- ddpcr_all %>%
  filter(
    Description %in% c("Biofilm-Control", "Biofilm-Spiked"),
    SampleWeek %in% paste("Week", 1:7),
    Probe %in% c("NDM", "YFP", "vanA", "tetA")
  ) %>%
  mutate(
    Description = factor(Description, levels = c("Biofilm-Control", "Biofilm-Spiked")),
    SampleWeek = factor(SampleWeek, levels = paste("Week", 1:7)),
    Probe = factor(Probe, levels = c("NDM", "YFP", "vanA", "tetA")),
    Concentration = as.numeric(Concentration),
    Concentration_plot = ifelse(is.na(Concentration), NA,
                                ifelse(Concentration <= 0, 0.001, Concentration))
  )

summary_df <- biofilm_data %>%
  group_by(Description, Probe, SampleWeek) %>%
  summarise(
    GeoMean = exp(mean(log(Concentration_plot), na.rm = TRUE)),
    SDlog = sd(log(Concentration_plot), na.rm = TRUE),
    N = sum(!is.na(Concentration_plot)),
    .groups = "drop"
  ) %>%
  mutate(
    Lower = GeoMean / exp(SDlog),
    Upper = GeoMean * exp(SDlog)
  )

loq_df <- data.frame(
  Probe = factor(c("NDM", "YFP", "vanA", "tetA"),
                 levels = c("NDM", "YFP", "vanA", "tetA")),
  LOQ = c(
    1010,    # NDM
    1610,         # YFP
    29.4,   # vanA
    307    # tetA
  )
) %>%
  mutate(
    LOQ_plot = ifelse(LOQ <= 0, 0.001, LOQ)
  )

p1 <- ggplot(
  summary_df,
  aes(
    x = SampleWeek,
    y = GeoMean,
    color = Probe,
    shape = Probe,
    group = interaction(Description, Probe),
    linetype = Description
  )
) +
  geom_line(linewidth = 1.6) +
  geom_point(size = 7) +
  geom_errorbar(
    aes(ymin = Lower, ymax = Upper),
    width = 0.25,
    linewidth = 1.2
  ) +
  geom_hline(
    data = loq_df,
    aes(
      yintercept = LOQ_plot,
      color = Probe,
      linetype = "LOQ"
    ),
    linewidth = 1.2,
    show.legend = TRUE
  ) +
  annotate(
    "text",
    x = 1.0,
    y = 10,
    label = "No visible biofilm formed",
    angle = 90,
    fontface = "italic",
    size = 9
  ) +
  scale_linetype_manual(
    name = "Line Type",
    values = c(
      "Biofilm-Control" = "dashed",
      "Biofilm-Spiked" = "solid",
      "LOQ" = "dotdash"
    )
  ) +
  scale_color_manual(
    name = "Target Gene", 
    values = cols,
    labels = c("NDM" = "NDM", "YFP" = "YFP", "vanA" = expression(italic("vanA")), "tetA" = expression(italic("tetA")))
  ) + 
  scale_shape_manual(
    name = "Target Gene", 
    values = shps,
    labels = c("NDM" = "NDM", "YFP" = "YFP", "vanA" = expression(italic("vanA")), "tetA" = expression(italic("tetA")))
  ) + 
  guides(
    linetype = guide_legend(order = 1),
    color = guide_legend(order = 2),
    shape = guide_legend(order = 2)
  )+
  scale_y_log10(
    limits = c(1e-5, 1e7),
    breaks = 10^seq(-5, 7, by = 2),
    labels = c(
      expression(10^-5),
      expression(10^-3),
      expression(10^-1),
      expression(10^1),
      expression(10^3),
      expression(10^5),
      expression(10^7)
    )
  ) +
  labs(
    #title = "Temporal trends of target gene concentrations in biofilm samples",
    x = "Sample Week",
    y = expression(bold("Concentration (copies/cm"^2*")"))
  ) +
  scale_x_discrete(expand = expansion(add = c(0.4, 0.4))) +
  theme_classic(base_size = 28) +
  theme(
    plot.title = element_text(face = "bold", size = 26, hjust = 0.5),
    plot.subtitle = element_text(size = 26, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 26),
    axis.text.x = element_text(size = 26, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 26),
    legend.title = element_text(face = "bold", size = 26),
    legend.text = element_text(size = 26),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7)
  )

print(p1)

ggsave("Biofilm_All_Targets_log_lineplot_errorbars.png", p1, width = 14, height = 8, dpi = 600)
ggsave("Biofilm_All_Targets_log_lineplot_errorbars.pdf", p1, width = 14, height = 8)
###########
##########
##########
################################################################################
###############################################################################
################################################################################
# WEEK 8-10: # Biofilm-Control/Biofilm-Spiked: all target genes line plot with geometric mean and LOQ legend
library(ggplot2)
library(dplyr)
library(readxl)
library(scales)

ddpcr_1 <- read_excel("Israel_YFP_ARG_ddPCR.xlsx")
names(ddpcr_1)[names(ddpcr_1) == "Conc(copies/mL)"] <- "Concentration"

ddpcr_2 <- read_excel("Israel_YFP_ARG_ddPCR.xlsx", sheet = "Combined_VanATetA")
names(ddpcr_2)[names(ddpcr_2) == "Conc(copies/mL)"] <- "Concentration"

ddpcr_1 <- ddpcr_1 %>% select(Description, SampleWeek, Probe, Concentration)
ddpcr_2 <- ddpcr_2 %>% select(Description, SampleWeek, Probe, Concentration) %>% mutate(Concentration = as.numeric(Concentration) * 10)

ddpcr_all <- bind_rows(ddpcr_1, ddpcr_2)

cols <- c("NDM" = "steelblue", "YFP" = "firebrick", "vanA" = "darkgreen", "tetA" = "purple")
shps <- c("NDM" = 16, "YFP" = 17, "vanA" = 15, "tetA" = 18)

biofilm_data <- ddpcr_all %>% 
  filter(
    Description %in% c("Biofilm-Control", "Biofilm-Spiked"),
    SampleWeek %in% paste("Week", 8:10),
    Probe %in% c("NDM", "YFP", "vanA", "tetA")
  ) %>% 
  mutate(
    # Map the old names to your new panel names
    Description = case_match(Description,
                             "Biofilm-Control" ~ "Detachment/Sloughing",
                             "Biofilm-Spiked"  ~ "Colonization"),
    # FIXED: Levels must precisely match the string values assigned above
    Description = factor(Description, levels = c("Detachment/Sloughing", "Colonization")),
    SampleWeek = factor(SampleWeek, levels = paste("Week", 8:10)),
    Probe = factor(Probe, levels = c("NDM", "YFP", "vanA", "tetA")),
    Concentration = as.numeric(Concentration),
    Concentration_plot = ifelse(is.na(Concentration), NA, ifelse(Concentration <= 0, 0.001, Concentration))
  )

summary_df <- biofilm_data %>% 
  group_by(Description, Probe, SampleWeek) %>% 
  summarise(
    GeoMean = exp(mean(log(Concentration_plot), na.rm = TRUE)),
    SDlog = sd(log(Concentration_plot), na.rm = TRUE),
    N = sum(!is.na(Concentration_plot)),
    .groups = "drop"
  ) %>% 
  mutate(
    Lower = GeoMean / exp(SDlog),
    Upper = GeoMean * exp(SDlog)
  )

# Updated from LOD to LOQ thresholds
loq_df <- data.frame(
  Probe = factor(c("NDM", "YFP", "vanA", "tetA"), levels = c("NDM", "YFP", "vanA", "tetA")),
  LOQ = c(
    1010, # NDM
    1610,      # YFP
    29.4,     # vanA
    307     # tetA
  )
) %>% 
  mutate(LOQ_plot = ifelse(LOQ <= 0, 0.001, LOQ))

p1 <- ggplot(
  summary_df, 
  aes(
    x = SampleWeek, 
    y = GeoMean, 
    color = Probe, 
    shape = Probe, 
    group = interaction(Description, Probe), 
    linetype = Description
  )
) +
  geom_line(linewidth = 1.6) +
  geom_point(size = 7) +
  geom_errorbar(
    aes(ymin = Lower, ymax = Upper), width = 0.25, linewidth = 1.2
  ) +
  geom_hline(
    data = loq_df, 
    aes(
      yintercept = LOQ_plot, 
      color = Probe, 
      linetype = "LOQ"
    ), 
    linewidth = 1.2, 
    show.legend = TRUE
  ) +
  # Facet wrap splits data into individual side-by-side panels
  facet_wrap(~Description) + 
  # Legend names to reflect updated factor names to properly draw lines
  scale_linetype_manual(
    name = "Line Type", 
    values = c(
      "Detachment/Sloughing" = "dashed", 
      "Colonization"         = "solid", 
      "LOQ"                  = "dotdash"
    )
  ) +
  scale_color_manual(
    name = "Target Gene", 
    values = cols,
    labels = c("NDM" = "NDM", "YFP" = "YFP", "vanA" = expression(italic("vanA")), "tetA" = expression(italic("tetA")))
  ) + 
  scale_shape_manual(
    name = "Target Gene", 
    values = shps,
    labels = c("NDM" = "NDM", "YFP" = "YFP", "vanA" = expression(italic("vanA")), "tetA" = expression(italic("tetA")))
  ) + 
  guides(
    linetype = guide_legend(order = 1),
    color = guide_legend(order = 2),
    shape = guide_legend(order = 2)
  ) + 
  scale_y_log10(
    limits = c(1e-1, 1e7),
    breaks = 10^seq(-1, 7, by = 2),
    labels = c(
      expression(10^-1),
      expression(10^1),
      expression(10^3),
      expression(10^5),
      expression(10^7)
    )
  ) + 
  labs(
    #title = "Temporal trends in target gene concentrations after coupon swap", 
    x = "Sample Week", 
    y = expression(bold("Concentration (copies/cm"^2*")"))
  ) + 
  scale_x_discrete(expand = expansion(add = c(0.15, 0.15))) + 
  theme_classic(base_size = 28) + 
  theme(
    plot.title = element_text(face = "bold", size = 26, hjust = 0.5),
    plot.subtitle = element_text(size = 26, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 26),
    axis.text.x = element_text(size = 26, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 26),
    legend.title = element_text(face = "bold", size = 26),
    legend.text = element_text(size = 26),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    # Styling for the facet headers to match your attached style
    strip.background = element_rect(color = "black", fill = "gray90"),
    strip.text = element_text(face = "bold", size = 20)
  )

print(p1)
ggsave("Biofilm_WEEK8-10_All_Targets_log_lineplot_errorbars.png", p1, width = 14, height = 8, dpi = 600)
ggsave("Biofilm_WEEK8-10_All_Targets_log_lineplot_errorbars.pdf", p1, width = 14, height = 8)
###############################################################################################
#################################################################################################
##############################################################################################
# ==============================================================================
# Verify tetA Log10 Mean and SD for all groups (UPDATED NAMES)
# ==============================================================================
tetA_check_all <- ddpcr_all %>% 
  filter( 
    # Match your exact sample names here
    Description %in% c("Influent", "Bulkwater-Control", "Bulkwater-Spiked", "Biofilm-Control", "Biofilm-Spiked"), 
    SampleWeek %in% paste("Week", 1:7), 
    Probe == "tetA" 
  ) %>% 
  mutate( 
    Concentration = as.numeric(Concentration), 
    Concentration_safe = ifelse(Concentration <= 0, 0.001, Concentration) 
  ) %>% 
  # Group by Description to calculate statistics separately for each reactor type
  group_by(Description) %>% 
  summarise( 
    Log10_Mean = mean(log10(Concentration_safe), na.rm = TRUE), 
    Log10_SD   = sd(log10(Concentration_safe), na.rm = TRUE),
    .groups    = "drop"
  )

# Loop through the rows to dynamically print the formatted sentences
for(i in 1:nrow(tetA_check_all)) {
  msg <- sprintf(
    "tetA in %s remained at %.3f ± %.3f log10 copies/mL", 
    tetA_check_all$Description[i], 
    tetA_check_all$Log10_Mean[i], 
    tetA_check_all$Log10_SD[i]
  )
  print(msg)
}

###############################################################################################
#################################################################################################
##############################################################################################
# ==============================================================================
# Verify vanA Log10 Mean and SD for all groups
# ==============================================================================
vanA_check_all <- ddpcr_all %>% 
  filter( 
    Description %in% c("Influent", "Bulkwater-Control", "Bulkwater-Spiked", "Biofilm-Control", "Biofilm-Spiked"), 
    SampleWeek %in% paste("Week", 1:7), 
    Probe == "vanA" 
  ) %>% 
  mutate( 
    Concentration = as.numeric(Concentration), 
    Concentration_safe = ifelse(Concentration <= 0, 0.001, Concentration) 
  ) %>% 
  group_by(Description) %>% 
  summarise( 
    Log10_Mean = mean(log10(Concentration_safe), na.rm = TRUE), 
    Log10_SD   = sd(log10(Concentration_safe), na.rm = TRUE),
    .groups    = "drop"
  )

# Loop through the rows to dynamically print the formatted sentences
for(i in 1:nrow(vanA_check_all)) {
  msg <- sprintf(
    "vanA in %s remained at %.3f ± %.3f log10 copies/mL", 
    vanA_check_all$Description[i], 
    vanA_check_all$Log10_Mean[i], 
    vanA_check_all$Log10_SD[i]
  )
  print(msg)
}
###############################################################################################
#################################################################################################
##############################################################################################
# ==============================================================================
# Verify NDM Log10 Mean and SD for all groups
# ==============================================================================
NDM_check_all <- ddpcr_all %>% 
  filter( 
    Description %in% c("Influent", "Bulkwater-Control", "Bulkwater-Spiked", "Biofilm-Control", "Biofilm-Spiked"), 
    SampleWeek %in% paste("Week", 1:7), 
    Probe == "NDM" # TARGETS NDM
  ) %>% 
  mutate( 
    Concentration = as.numeric(Concentration), 
    Concentration_safe = ifelse(Concentration <= 0, 0.001, Concentration) 
  ) %>% 
  group_by(Description) %>% 
  summarise( 
    Log10_Mean = mean(log10(Concentration_safe), na.rm = TRUE), 
    Log10_SD   = sd(log10(Concentration_safe), na.rm = TRUE),
    .groups    = "drop"
  )

# Loop through the rows to dynamically print the formatted sentences
for(i in 1:nrow(NDM_check_all)) {
  msg <- sprintf(
    "NDM in %s remained at %.3f ± %.3f log10 copies/mL", 
    NDM_check_all$Description[i], 
    NDM_check_all$Log10_Mean[i], 
    NDM_check_all$Log10_SD[i]
  )
  print(msg)
}
###############################################################################################
#################################################################################################
##############################################################################################
# ==============================================================================
# Verify YFP Log10 Mean and SD for all groups
# ==============================================================================
YFP_check_all <- ddpcr_all %>% 
  filter( 
    Description %in% c("Influent", "Bulkwater-Control", "Bulkwater-Spiked", "Biofilm-Control", "Biofilm-Spiked"), 
    SampleWeek %in% paste("Week", 8:10), 
    Probe == "YFP" # TARGETS YFP
  ) %>% 
  mutate( 
    Concentration = as.numeric(Concentration), 
    Concentration_safe = ifelse(Concentration <= 0, 0.001, Concentration) 
  ) %>% 
  group_by(Description) %>% 
  summarise( 
    Log10_Mean = mean(log10(Concentration_safe), na.rm = TRUE), 
    Log10_SD   = sd(log10(Concentration_safe), na.rm = TRUE),
    .groups    = "drop"
  )

# Loop through the rows to dynamically print the formatted sentences
for(i in 1:nrow(YFP_check_all)) {
  msg <- sprintf(
    "YFP in %s remained at %.3f ± %.3f log10 copies/mL", 
    YFP_check_all$Description[i], 
    YFP_check_all$Log10_Mean[i], 
    YFP_check_all$Log10_SD[i]
  )
  print(msg)
}

################################################################
################################################################
# ============================================================
# BULK WATER: PAIRED WILCOXON SIGNED-RANK TEST
# Control vs Spiked, Weeks 1–7
# NDM, YFP, vanA, tetA
# Bonferroni correction for four comparisons
# ============================================================
# ============================================================
library(dplyr)
library(tidyr)

# Step 1: Filter data and assign replicate IDs (instead of averaging)
bulk_paired <- ddpcr_all %>%
  filter(
    Description %in% c("Bulkwater-Control", "Bulkwater-Spiked"),
    SampleWeek %in% paste("Week", 1:7),
    Probe %in% c("NDM", "YFP", "vanA", "tetA")
  ) %>%
  mutate(
    Concentration = as.numeric(Concentration),
    Week = as.numeric(gsub("Week ", "", SampleWeek))
  ) %>%
  # Group by Week, Probe, and Description to number individual replicates sequentially
  group_by(Probe, Week, Description) %>%
  mutate(Replicate = row_number()) %>%
  ungroup() %>%
  # Select only the columns needed for pivoting
  select(Probe, Week, Replicate, Description, Concentration) %>%
  # Pivot wider to create side-by-side columns for your pairs
  pivot_wider(
    names_from = Description,
    values_from = Concentration
  ) %>%
  arrange(Probe, Week, Replicate)

# Check individual replicate paired data
print(bulk_paired, n = Inf)

# Step 2: Paired Wilcoxon signed-rank test for each gene (using all replicates)
bulk_wilcox_results <- bulk_paired %>%
  group_by(Probe) %>%
  group_modify(~ {
    
    # Drop rows where either the Control or Spiked replicate is missing
    paired_data <- .x %>%
      filter(
        !is.na(`Bulkwater-Control`),
        !is.na(`Bulkwater-Spiked`)
      )
    
    # Paired Wilcoxon requires at least 2 complete pairs
    if (nrow(paired_data) >= 2) {
      
      test <- wilcox.test(
        paired_data$`Bulkwater-Control`,
        paired_data$`Bulkwater-Spiked`,
        paired = TRUE,
        exact = FALSE
      )
      
      tibble(
        n_pairs = nrow(paired_data), # This will now reflect total replicate pairs, not weeks
        V = unname(test$statistic),
        p_value = test$p.value
      )
      
    } else {
      
      tibble(
        n_pairs = nrow(paired_data),
        V = NA_real_,
        p_value = NA_real_
      )
    }
  }) %>%
  ungroup()

# Step 3: Bonferroni correction for four genes
bulk_wilcox_results <- bulk_wilcox_results %>%
  mutate(
    p_adjusted = p.adjust(
      p_value,
      method = "bonferroni"
    ),
    significance = case_when(
      is.na(p_adjusted) ~ "NA",
      p_adjusted < 0.001 ~ "***",
      p_adjusted < 0.01  ~ "**",
      p_adjusted < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

# Final results
print(bulk_wilcox_results)

# ============================================================
# BIOFILM: PAIRED WILCOXON SIGNED-RANK TEST (REPLICATE-LEVEL)
# Control vs Spiked, Weeks 2–7 (pre-exchange)
# NDM, YFP, vanA, tetA
# Bonferroni correction for four comparisons
# ============================================================

library(dplyr)
library(tidyr)

# Step 1: Filter biofilm data and assign replicate IDs (instead of averaging)
biofilm_paired <- ddpcr_all %>%
  filter(
    Description %in% c("Biofilm-Control", "Biofilm-Spiked"),
    SampleWeek %in% paste("Week", 2:7),
    Probe %in% c("NDM", "YFP", "vanA", "tetA")
  ) %>%
  mutate(
    Concentration = as.numeric(Concentration),
    Week = as.numeric(gsub("Week ", "", SampleWeek))
  ) %>%
  # Group by Week, Probe, and Description to number individual replicates sequentially
  group_by(Probe, Week, Description) %>%
  mutate(Replicate = row_number()) %>%
  ungroup() %>%
  # Select only the columns needed for pivoting
  select(Probe, Week, Replicate, Description, Concentration) %>%
  # Pivot wider to create side-by-side columns for your pairs
  pivot_wider(
    names_from = Description,
    values_from = Concentration
  ) %>%
  arrange(Probe, Week, Replicate)

# Check individual replicate paired data
print(biofilm_paired, n = Inf)

# Step 2: Paired Wilcoxon signed-rank test for each gene (using all replicates)
biofilm_wilcox_results <- biofilm_paired %>%
  group_by(Probe) %>%
  group_modify(~ {
    
    # Drop rows where either the Control or Spiked replicate is missing
    paired_data <- .x %>%
      filter(
        !is.na(`Biofilm-Control`),
        !is.na(`Biofilm-Spiked`)
      )
    
    # Paired Wilcoxon requires at least 2 complete pairs
    if (nrow(paired_data) >= 2) {
      
      test = wilcox.test(
        paired_data$`Biofilm-Control`,
        paired_data$`Biofilm-Spiked`,
        paired = TRUE,
        exact = FALSE
      )
      
      tibble(
        n_pairs = nrow(paired_data), # Reflects total replicate pairs across weeks 2-7
        V = unname(test$statistic),
        p_value = test$p.value
      )
      
    } else {
      
      tibble(
        n_pairs = nrow(paired_data),
        V = NA_real_,
        p_value = NA_real_
      )
    }
  }) %>%
  ungroup()

# Step 3: Bonferroni correction for four genes
biofilm_wilcox_results <- biofilm_wilcox_results %>%
  mutate(
    p_adjusted = p.adjust(
      p_value,
      method = "bonferroni"
    ),
    significance = case_when(
      is.na(p_adjusted) ~ "NA",
      p_adjusted < 0.001 ~ "***",
      p_adjusted < 0.01  ~ "**",
      p_adjusted < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

# Final results
print(biofilm_wilcox_results)


# ============================================================
# INFLUENT: SPEARMAN TEMPORAL TREND
# Weeks 1–10
# NDM, YFP, vanA, tetA
# ============================================================

library(dplyr)

# Step 1: Average technical replicates for each gene and week
influent_weekly <- ddpcr_all %>%
  filter(
    Description == "Influent",
    SampleWeek %in% paste("Week", 1:10),
    Probe %in% c("NDM", "YFP", "vanA", "tetA")
  ) %>%
  mutate(
    Concentration = as.numeric(Concentration),
    Week = as.numeric(gsub("Week ", "", SampleWeek))
  ) %>%
  group_by(Probe, Week) %>%
  summarise(
    Mean_Concentration = ifelse(
      all(is.na(Concentration)),
      NA_real_,
      mean(Concentration, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  arrange(Probe, Week)

print(influent_weekly)


# Step 2: Spearman correlation for each gene
influent_spearman <- influent_weekly %>%
  group_by(Probe) %>%
  group_modify(~ {
    
    dat <- .x %>%
      filter(!is.na(Mean_Concentration))
    
    if(nrow(dat) >= 3) {
      
      test <- cor.test(
        dat$Week,
        dat$Mean_Concentration,
        method = "spearman",
        exact = FALSE
      )
      
      tibble(
        n = nrow(dat),
        rho = unname(test$estimate),
        p_value = test$p.value
      )
      
    } else {
      
      tibble(
        n = nrow(dat),
        rho = NA_real_,
        p_value = NA_real_
      )
    }
  }) %>%
  ungroup()


# Step 3: BH correction for testing four genes
influent_spearman <- influent_spearman %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH"),
    significance = case_when(
      p_adjusted < 0.001 ~ "***",
      p_adjusted < 0.01  ~ "**",
      p_adjusted < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

print(influent_spearman)
###################################################

# ============================================================
# INFLUENT WASTEWATER: TEMPORAL PERMANOVA
# Weeks 1–10 | NDM, YFP, vanA, tetA
# Individual ddPCR technical replicate concentrations
# Restricted time-series permutations
# ============================================================

library(dplyr)
library(vegan)
library(permute)

# 1. Prepare individual replicate concentrations
influent_replicates <- influent_wide %>%
  mutate(Week = as.numeric(as.character(Week))) %>%
  arrange(Week)

# 2. Log10-transform individual concentrations
gene_matrix_log <- influent_replicates %>%
  select(NDM, YFP, vanA, tetA) %>%
  mutate(across(everything(), ~ log10(.x + 1)))

# 3. Calculate Euclidean distances
gene_distance <- vegdist(
  gene_matrix_log,
  method = "euclidean"
)

# 4. Define restricted temporal permutations
# Each week is a block containing three technical replicates
# Weeks are shifted cyclically, preserving their sequence

time_control <- how(
  plots = Plots(
    strata = factor(influent_replicates$Week),
    type = "series",
    mirror = FALSE
  ),
  within = Within(type = "none"),
  nperm = 999
)

# 5. Run temporal PERMANOVA
set.seed(123)

permanova_temporal <- adonis2(
  gene_distance ~ Week,
  data = influent_replicates,
  permutations = time_control
)

# 6. Display results
print(permanova_temporal)

# 7. Save results
write.csv(
  as.data.frame(permanova_temporal),
  "Influent_Temporal_PERMANOVA_Individual.csv"
)



