#install package CRAN
#install.packages('pcr')
#library(pcr)
#install.packages("tidyverse")
library(tidyverse)
## locate within the pcr package
#fl <- system.file('extdata', 'ct1.csv', package = 'pcr')
#Set working directory
setwd("H:/Downloads/QPCRAnalysisFiles")
#copy the ct1.csv file in the pcr package to my computer
#destination_path <- file.path(getwd(), 'ct1_extracted.csv')
#file.copy(fl, destination_path)
#Read file and clean data
qpcr_BWdata <- read.csv("16srRNA_06_07_24_Results_BioreactorBulkwaterSamples.csv", skip = 7, stringsAsFactors = FALSE)
#check/inspect data
head(qpcr_BWdata)
# To remove rows with zero or negative Concentration values
qpcr_BWdata <- qpcr_BWdata %>%
  filter(Quantity > 0)
##Representing the data in form of the well plate##
#ggplot(qpcr_BWdata, aes(x = column, y = row, label = SampleName)) +
#  geom_tile(colour = "black") +
#  geom_text()
#ggplot(qpcr_BWdata, aes(x = column, y = row, fill = Task, label = SampleName)) +
#  geom_tile(colour = "black") +
#  geom_text(aes(size =0)) +
#  geom_text() +
#  scale_y_discrete(limits = c("H", "G", "F", "E", "D", "C", "B", "A")) +
#  scale_x_continuous(position = "bottom", guide = guide_axis(position = "top"), breaks = 1:12)
#view(qpcr_BWdata)
#calculate the log10 of quantity
qpcr_BWdata <-qpcr_BWdata %>%
  mutate(log10_quantity = log10(Quantity))
# Check for and handle any missing Ct values
qpcr_BWdata <- qpcr_BWdata %>%
  filter(!is.na(Ct))
# Fit a linear model
lm_fit <- lm(Ct ~ log10_quantity, data = qpcr_BWdata)
#standard_data <- qpcr_BWdata %>% filter(Task == "STANDARD")
#lm_fit <- lm(Ct ~ log10_quantity, data = standard_data)
# Summary of the linear model
summary(lm_fit)
#calculate efficiency
slope <- coef(lm_fit)["log10_quantity"]
efficiency <- (10^(-1/slope) - 1) * 100
cat("Efficiency: ", efficiency, "%\n")
# Calculate R²
r_squared <- summary(lm_fit)$r.squared
# Create the standard curve plot
P1 <- ggplot(qpcr_BWdata, aes(x = log10_quantity, y = Ct, color=Task,)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "qPCR Standard Curve",
    x = "Log10(quantity)",
    y = "Ct Value", 
  ) +
  annotate("text", x = max(qpcr_BWdata$log10_quantity), y = min(qpcr_BWdata$Ct),
           label = paste0("R² = ", round(r_squared, 3)),
           hjust = 2.8, vjust = -1.7) +
  annotate("text", x = max(qpcr_BWdata$log10_quantity), y = min(qpcr_BWdata$Ct) + 1,
           label = paste0("Efficiency = ", round(efficiency, 2), "%"),
           hjust = 1.5, vjust = 1) +
  theme_minimal() +
  theme(aspect.ratio = 7/7) +
  theme(axis.title.y = element_text(color = "black"), axis.title.x = element_text(color = "black"))
# Print the plot
P1
#Absolute Abundance- Prediction of log10 concentration for samples using the linear model
qpcr_BWdata <- qpcr_BWdata %>%
  mutate(predicted_log10_quantity = ifelse(Task == 'UNKNOWN', (Ct - coef(lm_fit)["(Intercept)"]) / slope, NA),
         predicted_quantity = 10^predicted_log10_quantity,
         predicted_concentration_with_units = ifelse(Task == 'UNKNOWN',
                                                     paste(round(predicted_quantity, 2), "copies/μL"),
                                                     NA))
# Plot predicted concentrations for samples
Pp1 <- ggplot(qpcr_BWdata %>% filter(Task == "UNKNOWN"), aes(x = predicted_quantity, y = Ct)) +
  geom_point(size = 3, color = "red") +
  labs(
    title = "Sample Abundances",
    x = "Predicted Concentration, copies/μL ",
    y = "Ct Value"
  ) +
  theme_minimal()
# Print the sample abundance plot
print(Pp1)
##Calculating gene concentration Values in gc/ml##
qpcr_BWdata$predicted_quantity_copies_mL <- NA
#view(qpcr_BWdata)
qpcr_BWdata$predicted_quantity_copies_mL <- ((qpcr_BWdata$predicted_quantity/10) * (10/3) * 100* (100/25))
qpcr_BWdata <-qpcr_BWdata %>%
  mutate(log10_quantity_copies_mL = log10(predicted_quantity_copies_mL))
# Verify the data preparation for boxplot creation
# Specify the order of sample names
sample_order <- c('Day6S', 'Day6C', 'Day13', 'Day20', 'Day27', 'Day34', 'Day41', 'Day48')
# Convert SampleName to a factor with specified levels
qpcr_BWdata$Period <- factor(qpcr_BWdata$Period, levels = sample_order)
print(qpcr_BWdata %>% filter(Task == "UNKNOWN"))
# Create the box plot for sample concentrations against sample names
box_plot <- ggplot(qpcr_BWdata %>% filter(Task == "UNKNOWN"), aes(x = Period, y = log10_quantity_copies_mL)) +
  geom_boxplot() +
  labs(
    title = "Concentrations of Reactor Bulkwater",
    x = "Sample Day",
    y = "log10 16SrRNA Gene (copies/mL)"
  ) +
  theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", colour = "white"),
      plot.background = element_rect(fill = "white", colour = "white"),   # White background
      panel.grid.major = element_blank(),
      axis.line = element_line(colour = "black"),  # Keep axis lines
      axis.ticks = element_line(colour = "black"),  # Keep axis ticks
      plot.title = element_text(hjust = 0.5)  # Center the title
    )

# Print the box plot
print(box_plot)
# Save the plots with specified width and height
ggsave("qPCR_standard_curve.png", plot = p, width = 10, height = 5)
ggsave("Sample_abundances.png", plot = Pp, width = 10, height = 5)
ggsave("Sample Concentrations by Sample Name.png", plot = , box_plot, width = 10, height = 5)
#################################################################
#################################################################
#For the PLWSA Biofilm Samples
#copy the ct1.csv file in the pcr package to my computer
#Read file and clean data
qpcr_pwdata <- read.csv("16SrRNA_06_20_24PLWSABiofilmSamples.csv", skip = 7, stringsAsFactors = FALSE)
#check/inspect data
head(qpcr_pwdata)
# To remove rows with zero or negative Concentration values
qpcr_pwdata <- qpcr_pwdata %>%
  filter(Quantity > 0)
#calculate the log10 of quantity
qpcr_pwdata <-qpcr_pwdata %>%
  mutate(log10_quantity = log10(Quantity))
# Check for and handle any missing Ct values
qpcr_pwdata <- qpcr_pwdata %>%
  filter(!is.na(Ct))
# Fit a linear model
lm_fit <- lm(Ct ~ log10_quantity, data = qpcr_pwdata)
#standard_data <- qpcr_pwdata %>% filter(Task == "STANDARD")
#lm_fit <- lm(Ct ~ log10_quantity, data = standard_data)
# Summary of the linear model
summary(lm_fit)
#calculate efficiency
slope <- coef(lm_fit)["log10_quantity"]
efficiency <- (10^(-1/slope) - 1) * 100
cat("Efficiency: ", efficiency, "%\n")
# Calculate R²
r_squared <- summary(lm_fit)$r.squared
# Create the standard curve plot
P2 <- ggplot(qpcr_pwdata, aes(x = log10_quantity, y = Ct, color=Task,)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "qPCR Standard Curve",
    x = "Log10(quantity)",
    y = "Ct Value", 
  ) +
  annotate("text", x = max(qpcr_pwdata$log10_quantity), y = min(qpcr_pwdata$Ct),
           label = paste0("R² = ", round(r_squared, 3)),
           hjust = 2.8, vjust = -1.7) +
  annotate("text", x = max(qpcr_pwdata$log10_quantity), y = min(qpcr_pwdata$Ct) + 1,
           label = paste0("Efficiency = ", round(efficiency, 2), "%"),
           hjust = 1.5, vjust = 1) +
  theme_minimal() +
  theme(aspect.ratio = 1) +
  theme(axis.title.y = element_text(color = "black"), axis.title.x = element_text(color = "black"))
# Print the plot
print(P2)
#Absolute Abundance- Prediction of log10 concentration for samples using the linear model
qpcr_pwdata <- qpcr_pwdata %>%
  mutate(predicted_log10_quantity = ifelse(Task == 'UNKNOWN', (Ct - coef(lm_fit)["(Intercept)"]) / slope, NA),
         predicted_quantity = 10^predicted_log10_quantity,
         predicted_concentration_with_units = ifelse(Task == 'UNKNOWN',
                                                     paste(round(predicted_quantity, 2), "copies/μL"),
                                                     NA))
# Plot predicted concentrations for samples
Pp2 <- ggplot(qpcr_pwdata %>% filter(Task == "UNKNOWN"), aes(x = predicted_quantity, y = Ct)) +
  geom_point(size = 3, color = "red") +
  labs(
    title = "Sample Abundances",
    x = "Predicted Concentration, copies/μL ",
    y = "Ct Value"
  ) +
  theme_minimal()

# Print the sample abundance plot
print(Pp2)
# Verify the data preparation for boxplot creation
# Calculate the mean predicted concentration for samples
mean_predicted_concentration <- mean(qpcr_pwdata$predicted_log10_quantity, na.rm = TRUE)
print(mean_predicted_concentration)
print(qpcr_pwdata %>% filter(Task == "UNKNOWN"))
# Create the box plot for sample concentrations against sample names
box_plot <- ggplot(qpcr_pwdata %>% filter(Task == "UNKNOWN"), aes(x = Sample.Name, y = predicted_log10_quantity)) +
  geom_boxplot(fill = "skyblue") +  # Set a single fill color
  geom_hline(yintercept = mean_predicted_concentration, color = "red", linetype = "dashed", size = 1) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.6) +  # Add jittered data points
    labs(
    title = "PLWSA Biofilm 16S Gene Concentration",
    x = "Sample",
    y = "log10 16SrRNA Gene (copies/cotton swab)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),
    plot.background = element_rect(fill = "white", colour = "white"),   # White background
    panel.grid.major = element_blank(),
    axis.line = element_line(colour = "black"),  # Keep axis lines
    axis.ticks = element_line(colour = "black"),  # Keep axis ticks
    plot.title = element_text(hjust = 0.5),  # Center the title
  )+
  scale_y_continuous(
    limits = c(5, 6),  # Start y-axis from 0, no upper limit
    breaks = seq(5, 6, by = 0.25),  # Set specific graduations on the y-axis (every 1 unit from 0 to 10)
    expand = c(0, 0)  # Remove padding around the y-axis limits
  )
print(box_plot)
# Save the plots with specified width and height
ggsave("qPCR_standard_curve.png", plot = Pp, width = 10, height = 5)
ggsave("Sample_abundances.png", plot = Pp, width = 10, height = 5)
ggsave("PLWSA Biofilm 16S Concentration.png", plot =box_plot, width = 10, height = 5)
#scale_y continuous or log
# Print the box plot
#################################################################
#################################################################
#For the Bioreactor Biofilm Samples
#copy the ct1.csv file in the pcr package to my computer
#Read file and clean data
qpcr_bbdata <- read.csv("16srRNA_06_29_24BioreactorBiofilmSamples.csv", skip = 7, stringsAsFactors = FALSE)
#check/inspect data
head(qpcr_bbdata)
# To remove rows with zero or negative Concentration values
qpcr_bbdata <- qpcr_bbdata %>% 
  filter(Quantity > 0)
##Representing the data in form of the well plate##
#ggplot(qpcr_bbdata, aes(x = column, y = row, label = SampleName)) +
#  geom_tile(colour = "black") +
#  geom_text()
#ggplot(qpcr_bbdata, aes(x = column, y = row, fill = Task, label = SampleName)) +
#  geom_tile(colour = "black") +
#  geom_text(aes(size =0)) +
#  geom_text() +
#  scale_y_discrete(limits = c("H", "G", "F", "E", "D", "C", "B", "A")) +
#  scale_x_continuous(position = "bottom", guide = guide_axis(position = "top"), breaks = 1:12)
#view(qpcr_BWdata)
#calculate the log10 of quantity
qpcr_bbdata <- qpcr_bbdata %>%
  mutate(log10_quantity = log10(Quantity))
# Check for and handle any missing Ct values
qpcr_bbdata <- qpcr_bbdata %>%
  filter(!is.na(Ct))
# Fit a linear model
lm_fit <- lm(Ct ~ log10_quantity, data = qpcr_bbdata)
#standard_data <- qpcr_BWdata %>% filter(Task == "STANDARD")
#lm_fit <- lm(Ct ~ log10_quantity, data = standard_data)
# Summary of the linear model
summary(lm_fit)
#calculate efficiency
slope <- coef(lm_fit)["log10_quantity"]
efficiency <- (10^(-1/slope) - 1) * 100
cat("Efficiency: ", efficiency, "%\n")
# Calculate R²
r_squared <- summary(lm_fit)$r.squared
# Create the standard curve plot
p3 <- ggplot(qpcr_bbdata, aes(x = log10_quantity, y = Ct, color=Task,)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "qPCR Standard Curve",
    x = "Log10(quantity)",
    y = "Ct Value", 
  ) +
  annotate("text", x = max(qpcr_bbdata$log10_quantity), y = min(qpcr_bbdata$Ct),
           label = paste0("R² = ", round(r_squared, 3)),
           hjust = 2.8, vjust = -1.7) +
  annotate("text", x = max(qpcr_bbdata$log10_quantity), y = min(qpcr_bbdata$Ct) + 1,
           label = paste0("Efficiency = ", round(efficiency, 2), "%"),
           hjust = 1.5, vjust = 1) +
  theme_minimal() +
  theme(aspect.ratio = 7/7) +
  theme(axis.title.y = element_text(color = "black"), axis.title.x = element_text(color = "black"))
# Print the plot
print (p3)
#Absolute Abundance- Prediction of log10 concentration for samples using the linear model
qpcr_bbdata <- qpcr_bbdata %>%
  mutate(predicted_log10_quantity = ifelse(Task == 'UNKNOWN', (Ct - coef(lm_fit)["(Intercept)"]) / slope, NA),
         predicted_quantity = 10^predicted_log10_quantity,
         predicted_concentration_with_units = ifelse(Task == 'UNKNOWN',
                                                     paste(round(predicted_quantity, 2), "copies/μL"),
                                                     NA))
# Plot predicted concentrations for samples
pp3 <- ggplot(qpcr_bbdata %>% filter(Task == 'UNKNOWN'), aes(x = predicted_quantity, y = Ct)) +
  geom_point(size = 3, color = "red") +
  labs(
    title = "Sample Abundances",
    x = "Predicted Concentration, copies/μL ",
    y = "Ct Value"
  ) +
  theme_minimal()

# Print the sample abundance plot
print(pp3)
##Calculating gene concentration Values in gc/ml##
qpcr_bbdata$log10_quantity_copies_mm <- NA
#view(qpcr_bbdata)
#Convert from copies/uL to copies/mL
qpcr_bbdata$predicted_quantity_copies_mm <- ((qpcr_bbdata$Quantity/10) * (10/3) * 100* (100/625))
qpcr_bbdata <-qpcr_bbdata %>%
  mutate(log10_quantity_copies_mm = log10(predicted_quantity_copies_mm))
# Verify the data preparation for boxplot creation
# Calculate the mean predicted concentration for samples
#mean_predicted_concentration <- mean(qpcr_bbdata$predicted_quantity, na.rm = TRUE)
#print(mean_predicted_concentration)
# Specify the order of sample names
sample_order <- c('DayS11', 'Day13', 'Day20', 'Day27', 'Day34', 'Day41', 'Day 48')
# Convert SampleName to a factor with specified levels
qpcr_bbdata$Period <- factor(qpcr_bbdata$Period, levels = sample_order)
print(qpcr_ARGdata %>% filter(SampleType == "Sample"))
# Create the box plot for sample concentrations against sample names
box_plot <- ggplot(qpcr_bbdata %>% filter(Task == "UNKNOWN"), aes(x = Period, y = log10_quantity_copies_mm)) +
  geom_boxplot() +
  labs(
    title = "Concentration of Bioreactor Biofilm",
    x = "Sample Day",
    y = "log10 16SrRNA Gene (copies/mm^2)"
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),
    plot.background = element_rect(fill = "white", colour = "white"),   # White background
    panel.grid.major = element_blank(),
    axis.line = element_line(colour = "black"),  # Keep axis lines
    axis.ticks = element_line(colour = "black"),  # Keep axis ticks
    plot.title = element_text(hjust = 0.5)  # Center the title
  )+
  scale_y_continuous(
    limits = c(0, NA),  # Start y-axis from 0, no upper limit
    breaks = seq(0, 10, by = 1)  # Set specific graduations on the y-axis (every 1 unit from 0 to 10)
  )
# Print the box plot
print(box_plot)
# Save the plots with specified width and height
ggsave("qPCR_standard_curve.png", plot = p, width = 10, height = 5)
ggsave("Sample_abundances.png", plot = Pp, width = 10, height = 5)
ggsave("Concentration of Bioreactor Biofilm.png", plot =box_plot, width = 10, height = 5)
##############################################################################################
#################################################################################################
############################################################################################
# Specify the file path
#file_path <- "my_data.csv"

# Save the dataframe to a CSV file
#write.csv(my_data, file = file_path, row.names = FALSE)

# This will save 'my_data' as 'my_data.csv' in your working directory.
# You can specify a different file path if needed.
file_path <- "qpcr_BWdata.csv"
write.csv(qpcr_BWdata, file = file_path, row.names = FALSE)
file_path <- "qpcr_bbdata.csv"
write.csv(qpcr_bbdata, file = file_path, row.names = FALSE)
##############################################################################################
#################################################################################################
############################################################################################
#plotting BWdata and bbdata boxplots on one plot
# Sample Data
#pcr_data <- data.frame(
  #SampleName = rep(c('Sample1', 'Sample2', 'Sample3', 'Sample4', 'Sample5'), each = 2),
  #ValueType = rep(c('Type1', 'Type2'), times = 5),
  #Value = c(5, 6, 3, 4, 6, 7, 4, 5, 7, 8)
#)
qpcr_combineddata <- read.csv("qpcr-bbdata_and_BWdataComibined.csv", stringsAsFactors = FALSE)
sample_order <- c('D6', 'D13', 'D20', 'D27', 'D34', 'D41', 'D48')
# Convert SampleName to a factor with specified levels
qpcr_combineddata$Period <- factor(qpcr_combineddata$Period, levels = sample_order)
# Create the side-by-side boxplot
box_plot_side_by_side <- ggplot(qpcr_combineddata, aes(x = Period, y = log10_quantity_copies_mL , fill = Sample))  +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.6) +  # Add jittered data points
  labs(
    title = "16S Concentration of Bulk water and Biofilm in Reactor 1 (Control)",
    x = "Sample Day",
    y = "log10 16SrRNA Gene"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),  # White background
    plot.background = element_rect(fill = "white", colour = "white"),   # White background
    panel.grid.major = element_blank(),
    axis.line = element_line(colour = "black"),  # Keep axis lines
    axis.ticks = element_line(colour = "black"),  # Keep axis ticks
    plot.title = element_text(hjust = 0.1)  # Center the title
    )+
  
  scale_y_continuous(
    limits = c(0.0, 8),  # Start y-axis from 0, no upper limit
    breaks = seq(0, 8, by = 2),  # Set specific graduations on the y-axis (every 1 unit from 0 to 10)
    expand = c(0, 0)  # Remove padding around the y-axis limits
  )

# Print the boxplot
print(box_plot_side_by_side)

# Save the plot
ggsave("box_plot_side_by_side.png", plot = box_plot_side_by_side, width = 8, height = 6)
#################################################################
#################################################################
#For the ARGs in Bioreactor Bulkwater Samples
#copy the ct1.csv file in the pcr package to my computer
#Read file and clean data
qpcr_ARGdata <- read.csv("Israel08232024_FirstBioreactorDataddPCRARGConcentration.csv", stringsAsFactors = FALSE)
#qpcr_ARGdata <- read_excel("Israel08232024_FirstBioreactorDataddPCRARGConcentration.xlsx")
#check/inspect data
head(qpcr_ARGdata)
qpcr_ARGdata$log10_quantity_copies <- NA
#view(qpcr_bbdata)
qpcr_ARGdata <-qpcr_ARGdata %>%
  mutate(log10_quantity_copies = log10(Converted.Copy.Number))
# Verify the data preparation for boxplot creation
# Specify the order of sample names
sample_order <- c('D13', 'D20', 'D27', 'D34', 'D41',  'D48')
# Convert SampleName to a factor with specified levels
qpcr_ARGdata$Period <- factor(qpcr_ARGdata$Period, levels = sample_order)
#print(qpcr_ARGdata %>% filter(SampleType == "Sample"))
# Assume the lowest detection concentration in copies/µL
lowest_detection_copies_per_uL <- min(qpcr_ARGdata$Conc.copies.µL[qpcr_ARGdata$Conc.copies.µL > 0])
print(lowest_detection_copies_per_uL)
# Calculate 1% of the lowest detection concentration
detection_limit_copies_per_uL <- 0.01 * lowest_detection_copies_per_uL
# Convert to copies/mL 
detection_limit_copies_per_mL <- detection_limit_copies_per_uL/40 * (40/5) * 10* (100/625)
# Print the detection limit for verification
print(paste("Detection limit in copies/mL:", detection_limit_copies_per_mL))
# Create the box plot for sample concentrations against sample names
box_plot <- ggplot(qpcr_ARGdata, aes(x = Period, y = log10_quantity_copies,)) +
  geom_boxplot(aes(fill = SampleType), alpha = 0.5, outlier.shape = NA) +
  #geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.6) +  # Add jittered data points
  geom_jitter(aes(shape = SampleType, color = ARG), width = 0.2, height = 0, size = 2) +  # Shape/Color for Quality
  geom_hline(yintercept = log10(detection_limit_copies_per_mL), linetype = "dashed", color = "red") +  # Add detection limit line
  labs(
    title = "ARG Concentration in Reactor samples",
    x = "Sample Day",
    y = "log10 ARG Concentration"
  ) +
  scale_fill_manual(values = c("purple", "darkgreen")) +  # Custom colors for DataType legend
  scale_color_manual(values = c("red", "blue")) +  # Custom colors for Quality legend
  scale_shape_manual(values = c(15, 17, 18)) +  # Assigning shapes to DataType (16 for circle, 17 for triangle)
  theme_minimal(base_size = 18) +
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),
    plot.background = element_rect(fill = "white", colour = "white"),   # White background
    panel.grid.major = element_blank(),
    axis.line = element_line(colour = "black"),  # Keep axis lines
    axis.ticks = element_line(colour = "black"),  # Keep axis ticks
    plot.title = element_text(hjust = 0.1),  # Center the title
    panel.spacing = unit(1, "lines"),         # Increase space between facets
    axis.text.x = element_text(angle = 55, hjust = 1),  # Rotate x-axis labels
    axis.text = element_text(size = 12)  # Adjust text size
  )+
  guides(
    fill = guide_legend(title = "Sample Type"),  # First legend for DataType
    shape = guide_legend(title = "Sample Type"),   # Second legend for Quality
    color = guide_legend(title = "ARG")    # Ensure that color matches Quality legend
  )+
  scale_y_continuous(
    limits = c(0, 4),  # Start y-axis from 0, no upper limit
    breaks = seq(0, 10, by = 1),  # Set specific graduations on the y-axis (every 1 unit from 0 to 10)
    expand = c(0, 0)  # Remove padding around the y-axis limits
  )+
  facet_wrap(SampleType~ ARG, scales = "free") #we want cols equal to work and also colour  by sample type. Cols by ARGs
# Print the box plot
print(box_plot)
# Save the plots with specified width and height
ggsave("ARG.png", plot =box_plot, width = 10, height = 5)
################################################################################
################################################################################
################################################################################
#For the ARGs in the PLWSA Influent Samples
#copy the ct1.csv file in the pcr package to my computer
#Read file and clean data
qpcr_ARGdata <- read.csv("Israel08232024_FirstBioreactorPLWSAInfDataddPCRARGConcentration.csv", stringsAsFactors = FALSE)
#qpcr_ARGdata <- read_excel("Israel08232024_FirstBioreactorDataddPCRARGConcentration.xlsx")
#check/inspect data
head(qpcr_ARGdata)
qpcr_ARGdata$log10_quantity_copies <- NA
#view(qpcr_bbdata)
qpcr_ARGdata <-qpcr_ARGdata %>%
  mutate(log10_quantity_copies = log10(Converted.Copy.Number))
# Verify the data preparation for boxplot creation
# Specify the order of sample names
sample_order <- c('D7', 'D14', 'D21', 'D28', 'D35',  'D42')
# Convert SampleName to a factor with specified levels
qpcr_ARGdata$Period <- factor(qpcr_ARGdata$Period, levels = sample_order)
#print(qpcr_ARGdata %>% filter(SampleType == "Sample"))
# Assume the lowest detection concentration in copies/µL
lowest_detection_copies_per_uL <- min(qpcr_ARGdata$Conc.copies.µL[qpcr_ARGdata$Conc.copies.µL > 0])
print(lowest_detection_copies_per_uL)
# Calculate 1% of the lowest detection concentration
detection_limit_copies_per_uL <- 0.01 * lowest_detection_copies_per_uL
# Convert to copies/mL 
detection_limit_copies_per_mL <- detection_limit_copies_per_uL/40 * (40/5) * 100* (100/25)
# Print the detection limit for verification
print(paste("Detection limit in copies/mL:", detection_limit_copies_per_mL))
# Create the box plot for sample concentrations against sample names
box_plot <- ggplot(qpcr_ARGdata, aes(x = Period, y = log10_quantity_copies,)) +
  geom_boxplot(aes(fill = SampleType), alpha = 0.5, outlier.shape = NA) +
  #geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.6) +  # Add jittered data points
  geom_jitter(aes(shape = SampleType, color = ARG), width = 0.2, height = 0, size = 2) +  # Shape/Color for Quality
  geom_hline(yintercept = log10(detection_limit_copies_per_mL), linetype = "dashed", color = "red") +  # Add detection limit line
  labs(
    title = "ARG Concentration in Influent Samples",
    x = "Sample Day",
    y = "log10 ARG Concentration"
  ) +
  scale_fill_manual(values = c("grey", "darkgreen")) +  # Custom colors for DataType legend
  scale_color_manual(values = c("red", "blue")) +  # Custom colors for Quality legend
  scale_shape_manual(values = c(15, 17, 18)) +  # Assigning shapes to DataType (16 for circle, 17 for triangle)
  theme_minimal(base_size = 18) +
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),
    plot.background = element_rect(fill = "white", colour = "white"),   # White background
    panel.grid.major = element_blank(),
    axis.line = element_line(colour = "black"),  # Keep axis lines
    axis.ticks = element_line(colour = "black"),  # Keep axis ticks
    plot.title = element_text(hjust = 0.1),  # Center the title
    panel.spacing = unit(1, "lines"),         # Increase space between facets
    axis.text.x = element_text(angle = 55, hjust = 1),  # Rotate x-axis labels
    axis.text = element_text(size = 12)  # Adjust text size
  )+
  guides(
    fill = guide_legend(title = "Sample Type"),  # First legend for DataType
    shape = guide_legend(title = "Sample Type"),   # Second legend for Quality
    color = guide_legend(title = "ARG")    # Ensure that color matches Quality legend
  )+
  scale_y_continuous(
    limits = c(0, 4),  # Start y-axis from 0, no upper limit
    breaks = seq(0, 10, by = 1),  # Set specific graduations on the y-axis (every 1 unit from 0 to 10)
    expand = c(0, 0)  # Remove padding around the y-axis limits
  )+
  facet_wrap(~ ARG, scales = "free") #we want cols equal to work and also colour  by sample type. Cols by ARGs
# Print the box plot
print(box_plot)
# Save the plots with specified width and height
ggsave("ARGPLWSAInf.png", plot =box_plot, width = 10, height = 5)
###################################################################
####################################################################
###################################################################
#Using the Excel sheet instead for Reactor Biofilm and Bulkwater ARG concentration
#qpcr_ARGdata <- read.csv("Israel08232024_FirstBioreactorDataddPCRARGConcentration.csv", stringsAsFactors = FALSE)
#install.packages("readxl")
#library(readxl)
qpcr_ARGdata <- read_excel("Israel08232024_FirstBioreactorDataddPCRARGConcentration.xlsx")
#check/inspect data
head(qpcr_ARGdata)
qpcr_ARGdata$log10_quantity_copies <- NA
# Filter for only TypeA and TypeB
qpcr_ARGdata <- qpcr_ARGdata %>%
  filter(SampleType %in% c("Biofilm", "Bulkwater"))
#view(qpcr_bbdata)
qpcr_ARGdata <-qpcr_ARGdata %>%
  mutate(log10_quantity_copies = log10(Converted_Copy_Number))
# Verify the data preparation for boxplot creation
# Specify the order of sample names
sample_order <- c('D13', 'D20', 'D27', 'D34', 'D41',  'D48')
# Convert SampleName to a factor with specified levels
qpcr_ARGdata$Period <- factor(qpcr_ARGdata$Period, levels = sample_order)
#print(qpcr_ARGdata %>% filter(SampleType == "Sample"))
# Assume the lowest detection concentration in copies/µL
lowest_detection_copies_per_uL <- min(qpcr_ARGdata$'Conc(copies/µL)'[qpcr_ARGdata$'Conc(copies/µL)' > 0])
print(lowest_detection_copies_per_uL)
# Calculate 1% of the lowest detection concentration
detection_limit_copies_per_uL <- 0.01 * lowest_detection_copies_per_uL
# Convert to copies/mL 
detection_limit_copies_per_mL <- detection_limit_copies_per_uL/40 * (40/5) * 10* (100/625)
# Print the detection limit for verification
print(paste("Detection limit in copies/mL:", detection_limit_copies_per_mL))
# Create the box plot for sample concentrations against sample names
box_plot <- ggplot(qpcr_ARGdata, aes(x = Period, y = log10_quantity_copies,)) +
  geom_boxplot(aes(fill = SampleType), alpha = 0.5, outlier.shape = NA) +
  #geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.6) +  # Add jittered data points
  geom_jitter(aes(shape = SampleType, color = ARG), width = 0.2, height = 0, size = 2) +  # Shape/Color for Quality
  geom_hline(yintercept = log10(detection_limit_copies_per_mL), linetype = "dashed", color = "red") +  # Add detection limit line
  labs(
    title = "ARG Concentration in Reactor samples",
    x = "Sample Day",
    y = "log10 ARG Concentration"
  ) +
  scale_fill_manual(values = c("purple", "darkgreen")) +  # Custom colors for DataType legend
  scale_color_manual(values = c("red", "blue")) +  # Custom colors for Quality legend
  scale_shape_manual(values = c(15, 17, 18)) +  # Assigning shapes to DataType (16 for circle, 17 for triangle)
  theme_minimal(base_size = 18) +
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),
    plot.background = element_rect(fill = "white", colour = "white"),   # White background
    panel.grid.major = element_blank(),
    axis.line = element_line(colour = "black"),  # Keep axis lines
    axis.ticks = element_line(colour = "black"),  # Keep axis ticks
    plot.title = element_text(hjust = 0.1),  # Center the title
    panel.spacing = unit(1, "lines"),         # Increase space between facets
    axis.text.x = element_text(angle = 55, hjust = 1),  # Rotate x-axis labels
    axis.text = element_text(size = 12)  # Adjust text size
  )+
  guides(
    fill = guide_legend(title = "Sample Type"),  # First legend for DataType
    shape = guide_legend(title = "Sample Type"),   # Second legend for Quality
    color = guide_legend(title = "ARG")    # Ensure that color matches Quality legend
  )+
  scale_y_continuous(
    limits = c(0, 4),  # Start y-axis from 0, no upper limit
    breaks = seq(0, 10, by = 1),  # Set specific graduations on the y-axis (every 1 unit from 0 to 10)
    expand = c(0, 0)  # Remove padding around the y-axis limits
  )+
  facet_wrap(~ ARG, scales = "free") #we want cols equal to work and also colour  by sample type. Cols by ARGs
# Print the box plot
print(box_plot)
####################################################################################
#####################################################################################
#####################################################################################
# Uing excel format instead of CSV For the ARGs in the PLWSA Influent Samples
#copy the ct1.csv file in the pcr package to my computer
#install.packages("readxl")
#library(readxl)
qpcr_ARGdata <- read_excel("Israel08232024_FirstBioreactorDataddPCRARGConcentration.xlsx")
#check/inspect data
head(qpcr_ARGdata)
qpcr_ARGdata$log10_quantity_copies <- NA
# Filter for only TypeA and TypeB
qpcr_ARGdata <- qpcr_ARGdata %>%
  filter(SampleType %in% c("WTTP Inf."))
#view(qpcr_bbdata)
head(qpcr_ARGdata)
qpcr_ARGdata$log10_quantity_copies <- NA
#view(qpcr_bbdata)
qpcr_ARGdata <-qpcr_ARGdata %>%
  mutate(log10_quantity_copies = log10(Converted_Copy_Number))
# Verify the data preparation for boxplot creation
# Specify the order of sample names
sample_order <- c('D7', 'D14', 'D21', 'D28', 'D35',  'D42')
# Convert SampleName to a factor with specified levels
qpcr_ARGdata$Period <- factor(qpcr_ARGdata$Period, levels = sample_order)
#print(qpcr_ARGdata %>% filter(SampleType == "Sample"))
# Assume the lowest detection concentration in copies/µL
lowest_detection_copies_per_uL <- min(qpcr_ARGdata$'Conc(copies/µL)'[qpcr_ARGdata$'Conc(copies/µL)'> 0])
print(lowest_detection_copies_per_uL)
# Calculate 1% of the lowest detection concentration
detection_limit_copies_per_uL <- 0.01 * lowest_detection_copies_per_uL
# Convert to copies/mL 
detection_limit_copies_per_mL <- detection_limit_copies_per_uL/40 * (40/5) * 100* (100/25)
# Print the detection limit for verification
print(paste("Detection limit in copies/mL:", detection_limit_copies_per_mL))
# Create the box plot for sample concentrations against sample names
box_plot <- ggplot(qpcr_ARGdata, aes(x = Period, y = log10_quantity_copies,)) +
  geom_boxplot(aes(fill = SampleType), alpha = 0.5, outlier.shape = NA) +
  #geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.6) +  # Add jittered data points
  geom_jitter(aes(shape = SampleType, color = ARG), width = 0.2, height = 0, size = 2) +  # Shape/Color for Quality
  geom_hline(yintercept = log10(detection_limit_copies_per_mL), linetype = "dashed", color = "red") +  # Add detection limit line
  labs(
    title = "ARG Concentration in Influent Samples",
    x = "Sample Day",
    y = "log10 ARG Concentration"
  ) +
  scale_fill_manual(values = c("grey", "darkgreen")) +  # Custom colors for DataType legend
  scale_color_manual(values = c("red", "blue")) +  # Custom colors for Quality legend
  scale_shape_manual(values = c(15, 17, 18)) +  # Assigning shapes to DataType (16 for circle, 17 for triangle)
  theme_minimal(base_size = 18) +
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),
    plot.background = element_rect(fill = "white", colour = "white"),   # White background
    panel.grid.major = element_blank(),
    axis.line = element_line(colour = "black"),  # Keep axis lines
    axis.ticks = element_line(colour = "black"),  # Keep axis ticks
    plot.title = element_text(hjust = 0.1),  # Center the title
    panel.spacing = unit(1, "lines"),         # Increase space between facets
    axis.text.x = element_text(angle = 55, hjust = 1),  # Rotate x-axis labels
    axis.text = element_text(size = 12)  # Adjust text size
  )+
  guides(
    fill = guide_legend(title = "Sample Type"),  # First legend for DataType
    shape = guide_legend(title = "Sample Type"),   # Second legend for Quality
    color = guide_legend(title = "ARG")    # Ensure that color matches Quality legend
  )+
  scale_y_continuous(
    limits = c(0, 4),  # Start y-axis from 0, no upper limit
    breaks = seq(0, 10, by = 1),  # Set specific graduations on the y-axis (every 1 unit from 0 to 10)
    expand = c(0, 0)  # Remove padding around the y-axis limits
  )+
  facet_wrap(~ ARG, scales = "free") #we want cols equal to work and also colour  by sample type. Cols by ARGs
# Print the box plot
print(box_plot)
# Save the plots with specified width and height
ggsave("ARGPLWSAInf.png", plot =box_plot, width = 10, height = 5)
###################################################################
####################################################################
#####################################################
##############
# Load the ggplot2 library
library(ggplot2)
install.packages("readxl")
library(readxl)
library(RColorBrewer)
# Install tidytext package if not already installed
if (!requireNamespace("tidytext", quietly = TRUE)) install.packages("tidytext")
library(tidytext)
bar_chart <- read.csv("AntibioticResistantGeneReview2.csv", skip = 1)
#check/inspect data
# Create a stacked bar chart
# Reorder Taxa levels within each Sample by RelativeAbundance
bar_chart <- bar_chart %>%
  mutate(Genus = factor(Genus, levels = names(sort(tapply(RelativeAbundance, Genus, mean), decreasing = TRUE))))
# Create custom color palette
#custom_colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(unique(bar_chart$Genus)))
bar_chart <- ggplot(bar_chart, aes(x = Sample, y = RelativeAbundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  #scale_fill_brewer(palette = "Set2") +  # Use a color palette for unique colors
  scale_x_discrete(labels = function(x) str_wrap(x, width = 12)) +  # Wrap x-axis text
  labs(
    title = NULL, # Remove plot title
    x = "Sample",
    y = "Relative Abundance (%)",
    fill = NULL # Remove legend title
  ) +
  theme_minimal() +
  theme(
    # Remove grid lines
    panel.grid = element_blank(),
    # Keep axis lines bold
    axis.line = element_line(color = "black", size = 1.2),
    # Add axis ticks and graduations
    axis.ticks = element_line(color = "black", size = 1),
    axis.ticks.length = unit(0.25, "cm"),  # Length of tick marks
    plot.title = element_text(hjust = 0.5, size = 16),  # Center and size title
    axis.title = element_text(size = 14),              # Axis title size
    axis.text = element_text(size = 11),               # Axis text size
    legend.title = element_text(size = 14),            # Legend title size
    legend.text = element_text(size = 12)              # Legend text size
  )

# Print the bar chart
print(bar_chart)

# Save the plot
ggsave("microbial_relative_abundance.png", plot = bar_chart, width = 10, height = 6)
###################################################################
####################################################################
#####################################################
##############
#Making a Relative Abundance figure with the abundance values known
# Load the ggplot2 library
library(ggplot2)
install.packages("readxl")
library(readxl)
library(RColorBrewer)
# Install tidytext package if not already installed
if (!requireNamespace("tidytext", quietly = TRUE)) install.packages("tidytext")
library(tidytext)
bar_chart <- read.csv("AntibioticResistantGeneReview2.csv", skip = 1)
#check/inspect data
# Create a stacked bar chart
# Reorder Taxa levels within each Sample by RelativeAbundance
bar_chart <- bar_chart %>%
  mutate(Genus = factor(Genus, levels = names(sort(tapply(RelativeAbundance, Phyla, mean), decreasing = TRUE))))
# Create custom color palette
#custom_colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(unique(bar_chart$Genus)))
bar_chart <- ggplot(bar_chart, aes(x = Sample, y = RelativeAbundance, fill = Phyla)) +
  geom_bar(stat = "identity", position = "stack") +
  #scale_fill_brewer(palette = "Set2") +  # Use a color palette for unique colors
  scale_x_discrete(labels = function(x) str_wrap(x, width = 12)) +  # Wrap x-axis text
  labs(
    title = NULL, # Remove plot title
    x = "Sample",
    y = "Relative Abundance (%)",
    fill = NULL # Remove legend title
  ) +
  theme_minimal() +
  theme(
    # Remove grid lines
    panel.grid = element_blank(),
    # Keep axis lines bold
    axis.line = element_line(color = "black", size = 1.2),
    # Add axis ticks and graduations
    axis.ticks = element_line(color = "black", size = 1),
    axis.ticks.length = unit(0.25, "cm"),  # Length of tick marks
    plot.title = element_text(hjust = 0.5, size = 16),  # Center and size title
    axis.title = element_text(size = 14),              # Axis title size
    axis.text = element_text(size = 11),               # Axis text size
    legend.title = element_text(size = 14),            # Legend title size
    legend.text = element_text(size = 12)              # Legend text size
  )

# Print the bar chart
print(bar_chart)

# Save the plot
ggsave("microbial_relative_abundance.png", plot = bar_chart, width = 10, height = 6)
###################################################################
####################################################################
#####################################################
##############
#ReviewPaper with Dr. Keenum
library(ggplot2)
install.packages("readxl")
library(readxl)
install.packages("tidyverse")  # For data manipulation and visualization
library(tidyverse)
library(RColorBrewer)
# Install tidytext package if not already installed
if (!requireNamespace("tidytext", quietly = TRUE)) install.packages("tidytext")
library(tidytext)
# Step 1: Load and prepare the data
data <- read_excel("MetaData for DrKeenumLabAntibiotic ResistantReviewPaper.xlsx")

# Ensure all columns except the first one (assumed to be Taxon or identifier) are numeric
numeric_data <- data %>%
  mutate(across(-1, ~as.numeric(.)))  # Ensure numeric format for all columns except the first one

# Calculate the total abundance for each sample (column sums excluding the Genus column)
total_abundance <- colSums(numeric_data[,-1], na.rm = TRUE)  # Ensure all columns are summed, skipping Genus
# Normalize data to relative abundance (keeping the first column as is)
relative_abundance <- numeric_data %>%
  mutate(across(-1, ~.x / total_abundance[match(cur_column(), names(total_abundance))]))  # Proper normalization
#relative_abundance <- numeric_data %>%
  #mutate(across(-1, ~ .x / total_abundance[colnames(.)])) # Normalize by sample totals

# Calculate total abundance across all samples for each taxon (row sums)
numeric_data$total_abundance <- rowSums(numeric_data[,-1], na.rm = TRUE)

# Identify the top 15 taxa by total abundance
top_15_taxa <- numeric_data %>%
  arrange(desc(total_abundance)) %>%
  slice(1:15) %>%
  pull(Genus)  # Replace 'Genus' with the column name for taxa or identifier

# Add a column categorizing taxa as "Top 15" or "Others"
relative_abundance <- relative_abundance %>%
  mutate(Genus = ifelse(Genus %in% top_15_taxa, Genus, "Others")) %>%
  group_by(Genus) %>%
  summarise(across(-1, sum, na.rm = TRUE)) %>%  # Summarize others
  ungroup()

# Ensure 'Others' comes last in the factor levels
relative_abundance$Genus <- factor(relative_abundance$Genus, levels = c(top_15_taxa, "Others"))

# Convert from wide to long format for ggplot
plot_data <- relative_abundance %>%
  pivot_longer(-Genus, names_to = "Sample", values_to = "Relative_Abundance")

# Set the Sample order explicitly based on the Excel sheet
sample_order <- names(numeric_data)[-1]  # Get sample order directly from the Excel sheet (excluding Genus)
plot_data$Sample <- factor(plot_data$Sample, levels = sample_order, ordered = TRUE)

# Generate 16 unique colors
set3_colors <- brewer.pal(12, "Set3")  # Generate 12 colors from Set3
extra_colors <- colorRampPalette(c("#FF5733", "#33FF57", "#3357FF", "#FFFF33"))(4)  # Add 4 more colors
color_palette <- c(set3_colors, extra_colors)

# Plot the bar chart
bar_plot <- ggplot(plot_data, aes(x = Sample, y = Relative_Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +  # Black bar boundaries
  scale_y_continuous(labels = scales::number_format(scale = 1)) +  # Show values between 0.0 and 1
  scale_fill_manual(values = color_palette) +  # Apply the unique colors
  theme_minimal() +
  labs(title = "Relative Abundance of Top 15 Genera",
       x = "Samples",
       y = "Relative Abundance",
       fill = "Genus") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  # Bold x-axis labels
    axis.text.y = element_text(face = "bold"),  # Bold y-axis labels
    axis.ticks = element_line(color = "black", size = 0.8),  # Bold ticklines
    axis.line = element_line(color = "black", size = 0.8),  # Bold axis lines
    panel.grid = element_blank(),  # Remove gridlines
    panel.border = element_blank()  # Remove panel border
  )
print(bar_plot)
# Save the plot
ggsave(
  filename = "Relative_Abundance_Bar_Plot.png",  # File name
  plot = bar_plot,                               # ggplot object
  width = 8,                                     # Width in inches
  height = 5,                                    # Height in inches
  dpi = 300                                      # Resolution in dots per inch
)
###################################################################
####################################################################
#####################################################
##############
library(ggplot2)
install.packages("readxl")
library(readxl)
install.packages("tidyverse")  # For data manipulation and visualization
library(tidyverse)
library(RColorBrewer)
# Install tidytext package if not already installed
if (!requireNamespace("tidytext", quietly = TRUE)) install.packages("tidytext")
library(tidytext)
# Step 1: Load and prepare the data
#data <- read_excel("MetaData.xlsx")
data <- read.csv("MetaData.csv")

# Ensure all columns except the first one (assumed to be Taxon or identifier) are numeric
numeric_data <- data %>%
  mutate(across(-1, ~as.numeric(.)))  # Ensure numeric format for all columns except the first one

# Calculate the total abundance for each sample (column sums excluding the Genus column)
total_abundance <- colSums(numeric_data[,-1], na.rm = TRUE)  # Ensure all columns are summed, skipping Genus
# Normalize data to relative abundance (keeping the first column as is)
relative_abundance <- numeric_data %>%
  mutate(across(-1, ~.x / total_abundance[match(cur_column(), names(total_abundance))]))  # Proper normalization
#relative_abundance <- numeric_data %>%
#mutate(across(-1, ~ .x / total_abundance[colnames(.)])) # Normalize by sample totals

# Calculate total abundance across all samples for each taxon (row sums)
numeric_data$total_abundance <- rowSums(numeric_data[,-1], na.rm = TRUE)

# Identify the top 15 taxa by total abundance
top_15_taxa <- numeric_data %>%
  arrange(desc(total_abundance)) %>%
  slice(1:15) %>%
  pull(Genus)  # Replace 'Genus' with the column name for taxa or identifier

# Add a column categorizing taxa as "Top 15" or "Others"
relative_abundance <- relative_abundance %>%
  mutate(Genus = ifelse(Genus %in% top_15_taxa, Genus, "Others")) %>%
  group_by(Genus) %>%
  summarise(across(-1, sum, na.rm = TRUE)) %>%  # Summarize others
  ungroup()

# Ensure 'Others' comes last in the factor levels
relative_abundance$Genus <- factor(relative_abundance$Genus, levels = c(top_15_taxa, "Others"))

# Convert from wide to long format for ggplot
plot_data <- relative_abundance %>%
  pivot_longer(-Genus, names_to = "Sample", values_to = "Relative_Abundance")

# Set the Sample order explicitly based on the Excel sheet
sample_order <- names(numeric_data)[-1]  # Get sample order directly from the Excel sheet (excluding Genus)
plot_data$Sample <- factor(plot_data$Sample, levels = sample_order, ordered = TRUE)

# Generate 16 unique colors
set3_colors <- brewer.pal(12, "Set3")  # Generate 12 colors from Set3
extra_colors <- colorRampPalette(c("#FF5733", "#33FF57", "#3357FF", "#FFFF33"))(4)  # Add 4 more colors
color_palette <- c(set3_colors, extra_colors)

# Plot the bar chart
bar_plot <- ggplot(plot_data, aes(x = Sample, y = Relative_Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +  # Black bar boundaries
  scale_y_continuous(labels = scales::number_format(scale = 1)) +  # Show values between 0.0 and 1
  scale_fill_manual(values = color_palette) +  # Apply the unique colors
  theme_minimal() +
  labs(title = "Relative Abundance of Top 15 Genera",
       x = "Samples",
       y = "Relative Abundance",
       fill = "Genus") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  # Bold x-axis labels
    axis.text.y = element_text(face = "bold"),  # Bold y-axis labels
    axis.ticks = element_line(color = "black", size = 0.8),  # Bold ticklines
    axis.line = element_line(color = "black", size = 0.8),  # Bold axis lines
    panel.grid = element_blank(),  # Remove gridlines
    panel.border = element_blank()  # Remove panel border
  )
print(bar_plot)
# Save the plot
ggsave(
  filename = "Relative_Abundance_Bar_Plot.png",  # File name
  plot = bar_plot,                               # ggplot object
  width = 8,                                     # Width in inches
  height = 5,                                    # Height in inches
  dpi = 300                                      # Resolution in dots per inch
)
###################################################################
####################################################################
#####################################################
##############
# Load necessary libraries
# Reorder 'Genus' so "Others" comes first, and the rest are ordered by abundance (smallest to largest)
library(ggplot2)
install.packages("readxl")
library(readxl)
install.packages("tidyverse")  # For data manipulation and visualization
library(tidyverse)
library(RColorBrewer)

# Install tidytext package if not already installed
if (!requireNamespace("tidytext", quietly = TRUE)) install.packages("tidytext")
library(tidytext)

# Step 1: Load and prepare the data
data <- read_excel("MetaData for DrKeenumLabAntibiotic ResistantReviewPaper.xlsx")

# Ensure all columns except the first one (assumed to be Taxon or identifier) are numeric
numeric_data <- data %>%
  mutate(across(-1, ~as.numeric(.)))  # Ensure numeric format for all columns except the first one

# Calculate the total abundance for each sample (column sums excluding the Genus column)
total_abundance <- colSums(numeric_data[,-1], na.rm = TRUE)  # Ensure all columns are summed, skipping Genus

# Normalize data to relative abundance (keeping the first column as is)
relative_abundance <- numeric_data %>%
  mutate(across(-1, ~.x / total_abundance[match(cur_column(), names(total_abundance))]))  # Proper normalization

# Calculate total abundance across all samples for each taxon (row sums)
numeric_data$total_abundance <- rowSums(numeric_data[,-1], na.rm = TRUE)

# Identify the top 15 taxa by total abundance
top_15_taxa <- numeric_data %>%
  arrange(desc(total_abundance)) %>%
  slice(1:15) %>%
  pull(Genus)  # Replace 'Genus' with the column name for taxa or identifier

# Add a column categorizing taxa as "Top 15" or "Others"
relative_abundance <- relative_abundance %>%
  mutate(Genus = ifelse(Genus %in% top_15_taxa, Genus, "Others")) %>%
  group_by(Genus) %>%
  summarise(across(-1, sum, na.rm = TRUE)) %>%  # Summarize others
  ungroup()

# Reorder 'Genus' so "Others" comes first, and the rest are ordered by abundance (smallest to largest)
ordered_genera <- relative_abundance %>%
  filter(Genus != "Others") %>%  # Exclude "Others" for sorting
  mutate(Total = rowSums(across(-Genus), na.rm = TRUE)) %>%  # Calculate total abundance per genus
  arrange(Total) %>%  # Sort by total abundance in ascending order
  pull(Genus)  # Extract the sorted genus names

relative_abundance$Genus <- factor(
  relative_abundance$Genus,
  levels = c("Others", ordered_genera)  # Place "Others" first, followed by sorted genera
)

# Convert from wide to long format for ggplot
plot_data <- relative_abundance %>%
  pivot_longer(-Genus, names_to = "Sample", values_to = "Relative_Abundance")

# Set the Sample order explicitly based on the Excel sheet
sample_order <- names(numeric_data)[-1]  # Get sample order directly from the Excel sheet (excluding Genus)
plot_data$Sample <- factor(plot_data$Sample, levels = sample_order, ordered = TRUE)

# Generate 16 unique colors
set3_colors <- brewer.pal(12, "Set3")  # Generate 12 colors from Set3
extra_colors <- colorRampPalette(c("#FF5733", "#33FF57", "#3357FF", "#FFFF33"))(4)  # Add 4 more colors
color_palette <- c(set3_colors, extra_colors)

# Create a named color palette with "Others" explicitly assigned to light grey
genus_colors <- c(
  setNames(color_palette, ordered_genera[ordered_genera != "Others"]),  # Keep the existing colors for genera
  "Others" = "grey100"  # Explicitly assign light grey to "Others"
)

# Plot the bar chart
bar_plot <- ggplot(plot_data, aes(x = Sample, y = Relative_Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +  # Black bar boundaries
  scale_y_continuous(labels = scales::number_format(scale = 1)) +  # Show values between 0.0 and 1
  scale_fill_manual(values = genus_colors) +  # Apply the custom colors
  theme_minimal() +
  labs(title = "Relative Abundance of Top 15 Genera",
       x = "Samples",
       y = "Relative Abundance",
       fill = "Genus") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  # Bold x-axis labels
    axis.text.y = element_text(face = "bold"),  # Bold y-axis labels
    axis.ticks = element_line(color = "black", size = 0.8),  # Bold ticklines
    axis.line = element_line(color = "black", size = 0.8),  # Bold axis lines
    panel.grid = element_blank(),  # Remove gridlines
    panel.border = element_blank()  # Remove panel border
  )
# Print the plot
print(bar_plot)

# Save the plot
ggsave(
  filename = "Relative_Abundance_Bar_Plot.png",  # File name
  plot = bar_plot,                               # ggplot object
  width = 8,                                     # Width in inches
  height = 5,                                    # Height in inches
  dpi = 300                                      # Resolution in dots per inch
)
###################################################################
####################################################################
#####################################################
##############
# Load necessary libraries
# Reorder 'Genus' so "Others" comes last, and the rest are ordered by abundance (largest to smallest)
library(ggplot2)
install.packages("readxl")
library(readxl)
install.packages("tidyverse")  # For data manipulation and visualization
library(tidyverse)
library(RColorBrewer)

# Install tidytext package if not already installed
if (!requireNamespace("tidytext", quietly = TRUE)) install.packages("tidytext")
library(tidytext)

# Step 1: Load and prepare the data
data <- read_excel("Supplementary File 3.xlsx")

# Ensure all columns except the first one (assumed to be Taxon or identifier) are numeric
numeric_data <- data %>%
  mutate(across(-1, ~as.numeric(.)))  # Ensure numeric format for all columns except the first one

# Calculate the total abundance for each sample (column sums excluding the Genus column)
total_abundance <- colSums(numeric_data[,-1], na.rm = TRUE)  # Ensure all columns are summed, skipping Genus

# Normalize data to relative abundance (keeping the first column as is)
relative_abundance <- numeric_data %>%
  mutate(across(-1, ~.x / total_abundance[match(cur_column(), names(total_abundance))]))  # Proper normalization

# Calculate total abundance across all samples for each taxon (row sums)
numeric_data$total_abundance <- rowSums(numeric_data[,-1], na.rm = TRUE)

# Identify the top 15 taxa by total abundance
top_15_taxa <- numeric_data %>%
  arrange(desc(total_abundance)) %>%
  slice(1:15) %>%
  pull(Genus)  # Replace 'Genus' with the column name for taxa or identifier

# Add a column categorizing taxa as "Top 15" or "Others"
relative_abundance <- relative_abundance %>%
  mutate(Genus = ifelse(Genus %in% top_15_taxa, Genus, "Others")) %>%
  group_by(Genus) %>%
  summarise(across(-1, sum, na.rm = TRUE)) %>%  # Summarize others
  ungroup()

# Reorder 'Genus' so "Others" comes last, and the rest are ordered by abundance (largest to smallest)
ordered_genera <- relative_abundance %>%
  filter(Genus != "Others") %>%  # Exclude "Others" for sorting
  mutate(Total = rowSums(across(-Genus), na.rm = TRUE)) %>%  # Calculate total abundance per genus
  arrange(desc(Total)) %>%  # Sort by total abundance in descending order
  pull(Genus)  # Extract the sorted genus names

relative_abundance$Genus <- factor(
  relative_abundance$Genus,
  levels = c(ordered_genera, "Others")  # Place sorted genera first, followed by "Others" 
)

# Convert from wide to long format for ggplot
plot_data <- relative_abundance %>%
  pivot_longer(-Genus, names_to = "Sample", values_to = "Relative_Abundance")

# Set the Sample order explicitly based on the Excel sheet
sample_order <- names(numeric_data)[-1]  # Get sample order directly from the Excel sheet (excluding Genus)
plot_data$Sample <- factor(plot_data$Sample, levels = sample_order, ordered = TRUE)

# Generate 16 unique colors
set3_colors <- brewer.pal(12, "Set3")  # Generate 12 colors from Set3
extra_colors <- colorRampPalette(c("#FF5733", "#33FF57", "#3357FF", "#FFFF33"))(4)  # Add 4 more colors
color_palette <- c(set3_colors, extra_colors)

# Create a named color palette with "Others" explicitly assigned to light grey
genus_colors <- c(
  setNames(color_palette, ordered_genera[ordered_genera != "Others"]),  # Keep the existing colors for genera
  "Others" = "grey100"  # Explicitly assign light grey to "Others"
)

# Plot the bar chart
bar_plot <- ggplot(plot_data, aes(x = Sample, y = Relative_Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +  # Black bar boundaries
  scale_y_continuous(labels = scales::number_format(scale = 1)) +  # Show values between 0.0 and 1
  scale_fill_manual(values = genus_colors) +  # Apply the custom colors
  theme_minimal() +
  labs(title = "Relative Abundance of Top 15 Genera",
       x = "Samples",
       y = "Relative Abundance",
       fill = "Genus") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  # Bold x-axis labels
    axis.text.y = element_text(face = "bold"),  # Bold y-axis labels
    axis.ticks = element_line(color = "black", linewidth = 0.8),  # Bold ticklines
    axis.line = element_line(color = "black", linewidth = 0.8),  # Bold axis lines
    panel.grid = element_blank(),  # Remove gridlines
    panel.border = element_blank()  # Remove panel border
  )
# Print the plot
print(bar_plot)

# Save the plot
ggsave(
  filename = "Relative_Abundance_Bar_Plot.png",  # File name
  plot = bar_plot,                               # ggplot object
  width = 8,                                     # Width in inches
  height = 5,                                    # Height in inches
  dpi = 300                                      # Resolution in dots per inch
)
##############################################################################################
#################################################################################################
############################################################################################
#Plotting the CONTROL-BULKWATER AND SPIKED-BULKWATER FOR SECOND BIROEACTOR PROJECT
#USING SCATTERRED PLOT. ACTUAL PLOT USED FOR PRESENTATION
#plotting "Bulk Wwater-C" vs "Bulk Wwater-S"
#Load required libraries
#library(ggplot2)
#library(dplyr)
#library(readxl)
# ==============================
qpcr_combineddata <- read_excel(
  "qPCR-BWdata-Biofilm-Data-SecondReactor.xlsx",
  sheet = "New_for_bar_Arrow"
)
# Define the sample order (Days)
# ==============================
sample_order <- c(
  "D6", "D13", "D20", "D27", "D34",
  "D41", "D48", "D55", "D62", "D69"
)
# Define Sample_Names to include
# ==============================
sample_names_of_interest <- c(
  "Bulk Wwater-C",
  "Bulk Wwater-S"
)
# Filter the dataset
# ==============================
filtered_data <- qpcr_combineddata %>%
  filter(
    Sample_Name %in% sample_names_of_interest,
    Period %in% sample_order
  )

# ==============================
# Convert Period (Day) to factor
# ==============================
filtered_data$Period <- factor(
  filtered_data$Period,
  levels = sample_order
)

# ==============================
# ADD: Day → Week mapping
# ==============================
day_to_week <- c(
  "D6"  = "1",
  "D13" = "2",
  "D20" = "3",
  "D27" = "4",
  "D34" = "5",
  "D41" = "6",
  "D48" = "7",
  "D55" = "8",
  "D62" = "9",
  "D69" = "10"
)

filtered_data <- filtered_data %>%
  mutate(Sample_Week = day_to_week[Period])

filtered_data$Sample_Week <- factor(
  filtered_data$Sample_Week,
  levels = paste(1:10)
)

# ==============================
# Create the scatter plot (by WEEK)
# ==============================
scatter_plot <- ggplot(
  filtered_data,
  aes(
    x = Sample_Week,
    y = log10_quantity_copies_mL,
    color = Sample_Name,
    shape = Sample_Name,
    group = Sample_Name
  )
) +
  geom_point(size = 6) +
  geom_line(size = 1) +
  labs(
    x = "Sample Week",
    y = expression(log[10] * " 16S rRNA Gene (Copies/mL)")
  ) +
  theme_minimal(base_size = 22) +
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),
    plot.background = element_rect(fill = "white", colour = "white"),
    panel.grid.major = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 32),
    axis.title.y = element_text(size = 32),
    axis.text.x = element_text(size = 32),
    axis.text.y = element_text(size = 32),
    legend.title = element_text(size = 32),
    legend.text = element_text(size = 32)
  ) +
  scale_color_manual(
    values = c(
      "Bulk Wwater-C" = "blue",
      "Bulk Wwater-S" = "brown"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Bulk Wwater-C" = 17,
      "Bulk Wwater-S" = 17
    )
  ) +
  scale_y_continuous(
    limits = c(0, 9),
    breaks = seq(0, 9, by = 2),
    expand = c(0, 0)
  )

# ==============================
# Display the plot
# ==============================
print(scatter_plot)
# Save the plot
#ggsave("scatter_plot.png", plot = scatter_plot, width = 8, height = 6)
#For poster presentation
ggsave("scatter_plot.png", plot = scatter_plot, 
       width = 16, height = 10, dpi = 300)
ggsave("scatter_plot.pdf", plot = scatter_plot, 
       width = 16, height = 10, dpi = 300)
######################################################################################################
######################################################################################################
# TotalBacteria in PLWSA Influent
# Load required libraries
library(readxl)
library(dplyr)
library(ggplot2)

# 1. Load the Excel data
qpcr_combineddata <- read_excel("qPCR-BWdata-Biofilm-Data-SecondReactor.xlsx")

# 2. Define the NEW sample order and names of interest
sample_order <- c('D0', 'D7', 'D14', 'D21', 'D28', 'D35', 'D42', 'D49', 'D56', 'D63')
sample_names_of_interest <- c("PLWSA Influent")

# 3. Filter and process the dataset
filtered_data <- qpcr_combineddata %>%
  filter(
    Sample_Name %in% sample_names_of_interest, 
    Period %in% sample_order
  ) %>%
  mutate(
    Period = factor(Period, levels = sample_order),
    Sample_Week = case_match(
      Period,
      "D0"  ~ "1", "D7"  ~ "2", "D14" ~ "3", "D21" ~ "4", "D28" ~ "5",
      "D35" ~ "6", "D42" ~ "7", "D49" ~ "8", "D56" ~ "9", "D63" ~ "10"
    ),
    Sample_Week = factor(Sample_Week, levels = as.character(1:10))
  )

# 4. Create the plot with Summary Statistics (Replaces geom_point/geom_line)
scatter_plot <- ggplot(filtered_data, aes(x = Sample_Week, y = log10_quantity_copies_mL, 
                                          color = Sample_Name, shape = Sample_Name, group = Sample_Name)) +
  # Adds error bars (Standard Deviation)
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), 
               geom = "errorbar", width = 0.2, linewidth = 1) +
  # Adds a line connecting the MEAN of each week
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  # Adds a large point at the MEAN of each week
  stat_summary(fun = mean, geom = "point", size = 8) + 
  labs(
    x = "Sample Week",
    y = "log10 16S rRNA Gene (Copies/mL)"
  ) +
  theme_minimal(base_size = 22) +
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),
    plot.background = element_rect(fill = "white", colour = "white"),
    panel.grid.major = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 32, margin = margin(t = 20)),
    axis.title.y = element_text(size = 32, margin = margin(r = 20)),
    axis.text.x = element_text(size = 32),
    axis.text.y = element_text(size = 32),
    legend.title = element_text(size = 32),
    legend.text = element_text(size = 32),
    legend.position = "right"
  ) +
  scale_color_manual(values = c("PLWSA Influent" = "blue")) +
  scale_shape_manual(values = c("PLWSA Influent" = 17)) +
  scale_y_continuous(
    limits = c(0.0, 9),  
    breaks = seq(0, 9, by = 3),
    expand = c(0, 0)
  )

# 5. Print and Save
print(scatter_plot)
ggsave("scatter_PLWSA_plot_clean.png", plot = scatter_plot, 
       width = 16, height = 10, dpi = 300)
######################################################################################################
#USING SCATTERRED PLOT
#plotting "Biofilm-C", "Biofilm-S"
#Load required libraries
#library(ggplot2)
#library(dplyr)
#library(readxl)

# Load the Excel data
qpcr_combineddata <- read_excel("qPCR-BWdata-Biofilm-Data-SecondReactor.xlsx", sheet = "New_for_bar_Arrow")

# Define the sample order
sample_order <- c('D6', 'D13', 'D20','D27', 'D34', 'D41', 'D48','D55','D62','D69')

# Define the list of Sample_Name values to include
sample_names_of_interest <- c("Biofilm-C", "Biofilm-S")

# Filter the dataset for the selected Sample_Names and Periods
filtered_data <- qpcr_combineddata %>%
  filter(Sample_Name %in% sample_names_of_interest, Period %in% sample_order)

# Convert Period to a factor with specified levels
filtered_data$Period <- factor(filtered_data$Period, levels = sample_order)

# ==============================
# ADD: Day → Week mapping
# ==============================
day_to_week <- c(
  "D6"  = "1",
  "D13" = "2",
  "D20" = "3",
  "D27" = "4",
  "D34" = "5",
  "D41" = "6",
  "D48" = "7",
  "D55" = "8",
  "D62" = "9",
  "D69" = "10"
)

filtered_data <- filtered_data %>%
  mutate(Sample_Week = day_to_week[Period])

filtered_data$Sample_Week <- factor(
  filtered_data$Sample_Week,
  levels = paste(1:10)
)

# Create the scatter plot with different shapes and colors
scatter_plot <- ggplot(filtered_data, aes(x = Sample_Week, y = log10_quantity_copies_mL, 
                                          color = Sample_Name, shape = Sample_Name, group = Sample_Name)) +
  geom_point(size = 6) +  # Adjust point size for better visibility
  geom_line(size = 1) +  # Add connecting line
  labs(
    #title = "16S rRNA Concentration in Reactor Biofilms",
    x = "Sample Week",
    y = expression("log"[10]~"16S rRNA Gene (Copies/cm"^2*")")
  ) +
  theme_minimal(base_size = 34) +
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),  # White background
    plot.background = element_rect(fill = "white", colour = "white"),   # White background
    panel.grid.major = element_blank(),
    axis.line = element_line(colour = "black"),  # Keep axis lines
    axis.ticks = element_line(colour = "black"),  # Keep axis ticks
    plot.title = element_text(hjust = 0.5),  # Center the title
    axis.title.x = element_text(size = 34),  # X-axis label size
    axis.title.y = element_text(size = 34),  # Y-axis label size
    axis.text.x = element_text(size = 34),  # X-axis tick labels
    axis.text.y = element_text(size = 34),  # Y-axis tick labels
    legend.title = element_text(size = 34),  # Legend title
    legend.text = element_text(size = 34)   # Legend text siz
  ) +
  scale_color_manual(values = c("Biofilm-C" = "blue", 
                                "Biofilm-S" = "brown")) +  # Set custom colors
  scale_shape_manual(values = c("Biofilm-C" = 16,  # Different shapes
                                "Biofilm-S" = 16)) +
  scale_y_continuous(
    limits = c(0.0, 9),  
    breaks = seq(0, 9, by = 2),  # Set specific graduations on the y-axis
    expand = c(0, 0)
  )

# Print the scatter plot
print(scatter_plot)
# Save the plot
#ggsave("scatter_plot.png", plot = scatter_Biofilm_plot, width = 8, height = 6)
#For poster presentation
ggsave("scatter_Biofilm_plot.png", plot = scatter_plot, 
       width = 16, height = 10, dpi = 300)
ggsave("scatter_Biofilm_plot.pdf", plot = scatter_plot, 
       width = 16, height = 10, dpi = 300)
######################################################################################################
#USING SCATTERRED PLOT
#scatter_Biofilm_Bulkwater_Spikedplot
#Load required libraries
#library(ggplot2)
#library(dplyr)
#library(readxl)

# Load the Excel data
# 1. Read Data
qpcr_combineddata <- read_excel("qPCR-BWdata-Biofilm-Data-SecondReactor.xlsx", sheet = "New_2_Bar_Arrow")

# Define the sample order
sample_order <- c('D6', 'D13', 'D20','D27', 'D34', 'D41', 'D48','D55','D62','D69')

# Define the list of Sample_Name values to include
sample_names_of_interest <- c("Bulk Wwater-S  (copies/mL)", "Biofilm-S (copies/cm^2)")

# Filter the dataset for the selected Sample_Names and Periods
filtered_data <- qpcr_combineddata %>%
  filter(Sample_Name %in% sample_names_of_interest,
         Period %in% sample_order,
         !is.na(log10_quantity_copies_mL))  # Remove rows with NA values

# Convert Period to a factor with specified levels
filtered_data$Period <- factor(filtered_data$Period, levels = sample_order)

# Set Sample_Name as a factor in the correct plotting order (to fix color mismatch)
filtered_data$Sample_Name <- factor(
  filtered_data$Sample_Name,
  levels = c("Bulk Wwater-S  (copies/mL)", "Biofilm-S (copies/cm^2)")
)

# Create the scatter plot with corrected legend order and colors
scatter_plot <- ggplot(filtered_data, aes(x = Period, y = log10_quantity_copies_mL, 
                                          color = Sample_Name, shape = Sample_Name, group = Sample_Name)) +
  geom_point(size = 6) +  # Adjust point size for better visibility
  geom_line(size = 1) +   # Add connecting line
  labs(
    x = "Sample Day",
    y = expression(log[10]*" 16S rRNA Gene")
  ) +
  theme_minimal(base_size = 32) +
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),
    plot.background = element_rect(fill = "white", colour = "white"),
    panel.grid.major = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 32),
    axis.title.y = element_text(size = 32),
    axis.text.x = element_text(size = 32),
    axis.text.y = element_text(size = 32),
    legend.title = element_text(size = 32),
    legend.text = element_text(size = 32)
  ) +
  scale_color_manual(
    values = c("Bulk Wwater-S  (copies/mL)" = "black", 
               "Biofilm-S (copies/cm^2)" = "brown"),
    labels = c("Bulk Wastewater-S (copies/mL)", 
               expression("Biofilm-S (copies/cm"^2*")"))
  ) +
  scale_shape_manual(
    values = c("Bulk Wwater-S  (copies/mL)" = 17, 
               "Biofilm-S (copies/cm^2)" = 16),
    labels = c("Bulk Wastewater-S (copies/mL)", 
               expression("Biofilm-S (copies/cm"^2*")"))
  ) +
  scale_y_continuous(
    limits = c(0.0, 9),  
    breaks = seq(0, 9, by = 3),
    expand = c(0, 0)
  )

# Print the plot
print(scatter_plot)
# Save the plot
#ggsave("scatter_plot.png", plot = scatter_Biofilm_plot, width = 8, height = 6)
#For poster presentation
ggsave("scatter_Biofilm_Bulkwater_Spikedplot.png", plot = scatter_plot, 
       width = 19, height = 10, dpi = 300)
######################################################################################################
#USING SCATTERRED PLOT
#scatter_Biofilm_Bulkwater_Control
#Load required libraries
#library(ggplot2)
#library(dplyr)
#library(readxl)
#setwd("H:/Downloads/QPCRAnalysisFiles/QPCR_SecondReactorProjectSamples")
# Load the Excel data
qpcr_combineddata <- read_excel("qPCR-BWdata-Biofilm-Data-SecondReactor.xlsx", sheet = "New_2_Bar_Arrow")

# Define the sample order
sample_order <- c('D6', 'D13', 'D20','D27', 'D34', 'D41', 'D48','D55','D62','D69')

# Define the list of Sample_Name values to include
sample_names_of_interest <- c("Bulk Wwater-C (copies/mL)", "Biofilm-C (copies/cm^2)")

# Filter the dataset for the selected Sample_Names and Periods
filtered_data <- qpcr_combineddata %>%
  filter(Sample_Name %in% sample_names_of_interest,
         Period %in% sample_order,
         !is.na(log10_quantity_copies_mL))  # Remove rows with NA values

# Convert Period to a factor with specified levels
filtered_data$Period <- factor(filtered_data$Period, levels = sample_order)

# Set Sample_Name as a factor in the correct plotting order (to fix color mismatch)
filtered_data$Sample_Name <- factor(
  filtered_data$Sample_Name,
  levels = c("Bulk Wwater-C (copies/mL)", "Biofilm-C (copies/cm^2)")
)

# Create the scatter plot with corrected legend order and colors
scatter_plot <- ggplot(filtered_data, aes(x = Period, y = log10_quantity_copies_mL, 
                                          color = Sample_Name, shape = Sample_Name, group = Sample_Name)) +
  geom_point(size = 6) +  # Adjust point size for better visibility
  geom_line(size = 1) +   # Add connecting line
  labs(
    x = "Sample Day",
    y = expression(log[10]*" 16S rRNA Gene")
  ) +
  theme_minimal(base_size = 32) +
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),
    plot.background = element_rect(fill = "white", colour = "white"),
    panel.grid.major = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 32),
    axis.title.y = element_text(size = 32),
    axis.text.x = element_text(size = 32),
    axis.text.y = element_text(size = 32),
    legend.title = element_text(size = 32),
    legend.text = element_text(size = 32)
  ) +
  scale_color_manual(
    values = c("Bulk Wwater-C (copies/mL)" = "black", 
               "Biofilm-C (copies/cm^2)" = "blue"),
    labels = c("Bulk Wastewater-C (copies/mL)", 
               expression("Biofilm-C (copies/cm"^2*")"))
  ) +
  scale_shape_manual(
    values = c("Bulk Wwater-C (copies/mL)" = 17, 
               "Biofilm-C (copies/cm^2)" = 16),
    labels = c("Bulk Wastewater-C (copies/mL)", 
               expression("Biofilm-C (copies/cm"^2*")"))
  ) +
  scale_y_continuous(
    limits = c(0.0, 9),  
    breaks = seq(0, 9, by = 3),
    expand = c(0, 0)
  )

# Print the scatter plot
print(scatter_plot)
# Save the plot
#ggsave("scatter_plot.png", plot = scatter_Biofilm_plot, width = 8, height = 6)
#For poster presentation
ggsave("scatter_Biofilm_Bulkwater_Controlplot.png", plot = scatter_plot, 
       width = 19, height = 10, dpi = 300)

###############################################################################
###############################################################################
# Paired Wilcoxon signed-rank test
# Biofilm Control vs Spiked — Pre-exchange (Weeks 2–7)
###############################################################################

library(readxl)
library(dplyr)
library(tidyr)

# ============================================================
# Step 1: Read the correct Excel sheet
# ============================================================

qpcr_combineddata <- read_excel(
  "qPCR-BWdata-Biofilm-Data-SecondReactor.xlsx",
  sheet = "For Biofilm Concentration"
)


# ============================================================
# Step 2: Prepare Biofilm Control vs Spiked data
# Pre-exchange = D13–D48 (Weeks 2–7)
# ============================================================

biofilm_wilcox_data <- qpcr_combineddata %>%
  
  filter(
    Sample_Name %in% c(
      "Biofilm-C (copies/cm^2)",
      "Biofilm-S (copies/cm^2)"
    ),
    Period %in% c(
      "D13", "D20", "D27",
      "D34", "D41", "D48"
    )
  ) %>%
  
  # Average technical replicates for each sampling point
  group_by(Period, Sample_Name) %>%
  summarise(
    Mean_16S = mean(
      log10_quantity_copies_mL,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  
  # Put Control and Spiked into separate columns
  pivot_wider(
    names_from = Sample_Name,
    values_from = Mean_16S
  )

# View data
print(biofilm_wilcox_data)


# ============================================================
# Step 3: Keep only complete Control-Spiked pairs
# ============================================================

biofilm_paired <- biofilm_wilcox_data %>%
  filter(
    !is.na(`Biofilm-C (copies/cm^2)`),
    !is.na(`Biofilm-S (copies/cm^2)`)
  )

# View paired observations
print(biofilm_paired)

# Number of paired sampling points
nrow(biofilm_paired)


# ============================================================
# Step 4: Paired Wilcoxon signed-rank test
# ============================================================

biofilm_wilcox <- wilcox.test(
  biofilm_paired$`Biofilm-C (copies/cm^2)`,
  biofilm_paired$`Biofilm-S (copies/cm^2)`,
  paired = TRUE,
  exact = FALSE
)

print(biofilm_wilcox)


# ============================================================
# Step 5: Calculate mean ± standard deviation
# ============================================================

# Control biofilm
mean_control_biofilm <- mean(
  biofilm_paired$`Biofilm-C (copies/cm^2)`,
  na.rm = TRUE
)

sd_control_biofilm <- sd(
  biofilm_paired$`Biofilm-C (copies/cm^2)`,
  na.rm = TRUE
)

# Spiked biofilm
mean_spiked_biofilm <- mean(
  biofilm_paired$`Biofilm-S (copies/cm^2)`,
  na.rm = TRUE
)

sd_spiked_biofilm <- sd(
  biofilm_paired$`Biofilm-S (copies/cm^2)`,
  na.rm = TRUE
)

# ============================================================
# Display results as Mean ± SD
# ============================================================

cat(sprintf(
  "Control biofilm = %.1f ± %.3f log10 copies/cm²\n",
  mean_control_biofilm,
  sd_control_biofilm
))

cat(sprintf(
  "Spiked biofilm = %.1f ± %.3f log10 copies/cm²\n",
  mean_spiked_biofilm,
  sd_spiked_biofilm
))
###############################################################################
##############################################################################
# ============================================================
# Paired Wilcoxon signed-rank test
# Bulk water Control vs Spiked
# Pre-exchange: Weeks 1–7
# Sheet: "For bulkwater Concentration"
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)

# Read correct worksheet
qpcr_combineddata <- read_excel(
  "qPCR-BWdata-Biofilm-Data-SecondReactor.xlsx",
  sheet = "For bulkwater Concentration"
)

# ============================================================
# Step 1: Average technical replicates for each sampling point
# ============================================================

bulk_wilcox_data <- qpcr_combineddata %>%
  filter(
    Sample_Name %in% c("Bulk water-C", "Bulk water-S"),
    Period %in% c(
      "D6", "D13", "D20", "D27",
      "D34", "D41", "D48"
    )
  ) %>%
  group_by(Period, Sample_Name) %>%
  summarise(
    Mean_16S = mean(log10_quantity_copies_mL, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Sample_Name,
    values_from = Mean_16S
  )

print(bulk_wilcox_data)


# ============================================================
# Step 2: Keep complete Control-Spiked pairs
# ============================================================

bulk_paired <- bulk_wilcox_data %>%
  filter(
    !is.na(`Bulk water-C`),
    !is.na(`Bulk water-S`)
  )

print(bulk_paired)

# Number of paired sampling points
cat("Number of pairs =", nrow(bulk_paired), "\n")


# ============================================================
# Step 3: Paired Wilcoxon signed-rank test
# ============================================================

bulk_wilcox <- wilcox.test(
  bulk_paired$`Bulk water-C`,
  bulk_paired$`Bulk water-S`,
  paired = TRUE,
  exact = FALSE
)

print(bulk_wilcox)


# ============================================================
# Step 4: Calculate mean ± standard deviation
# ============================================================

# Control bulk water
mean_control_bulk <- mean(
  bulk_paired$`Bulk water-C`,
  na.rm = TRUE
)

sd_control_bulk <- sd(
  bulk_paired$`Bulk water-C`,
  na.rm = TRUE
)

# Spiked bulk water
mean_spiked_bulk <- mean(
  bulk_paired$`Bulk water-S`,
  na.rm = TRUE
)

sd_spiked_bulk <- sd(
  bulk_paired$`Bulk water-S`,
  na.rm = TRUE
)

# ============================================================
# Display results as Mean ± SD
# ============================================================

cat(sprintf(
  "Control bulk-water = %.1f ± %.3f log10 copies/mL\n",
  mean_control_bulk,
  sd_control_bulk
))

cat(sprintf(
  "Spiked bulk-water = %.1f ± %.3f log10 copies/mL\n",
  mean_spiked_bulk,
  sd_spiked_bulk
))

################
###############
# Corrected - Bonferroni adjustment showing original vs corrected p-values

# 1. Collect the original p-values
raw_p_values <- c(
  Bulk_Water = bulk_wilcox$p.value,
  Biofilm    = biofilm_wilcox$p.value
)

# 2. Apply the Bonferroni correction
corrected_p_values <- p.adjust(raw_p_values, method = "bonferroni")

# 3. Combine into a clean table to show both
p_value_summary <- data.frame(
  Original_p  = raw_p_values,
  Adjusted_p  = corrected_p_values
)

# Print the results
print(p_value_summary)

# ============================================================
# Mean and SD for PLWSA INFLUENT (Weeks 1–10)
# ============================================================

library(readxl)
library(dplyr)

# 1. Read data and calculate average per sampling day
influent_data <- read_excel("qPCR-BWdata-Biofilm-Data-SecondReactor.xlsx", sheet = "For bulkwater Concentration") %>%
  filter(
    Sample_Name == "PLWSA Influent",
    Period %in% c("D0", "D7", "D14", "D21", "D28", "D35", "D42", "D49", "D56", "D63")
  ) %>%
  group_by(Period) %>%
  summarise(Mean_16S = mean(log10_quantity_copies_mL, na.rm = TRUE), .groups = "drop")

# 2. Calculate overall Mean and SD across the 10 weeks
overall_mean <- mean(influent_data$Mean_16S, na.rm = TRUE)
overall_sd   <- sd(influent_data$Mean_16S, na.rm = TRUE)

# 3. Print Summary Output
cat(sprintf("PLWSA Influent = %.1f ± %.3f log10 copies/mL\n", overall_mean, overall_sd))


