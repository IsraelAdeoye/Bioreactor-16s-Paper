# Load necessary library
library(ggplot2)

# Sample Data
pcr_data <- data.frame(
  SampleName = rep(c('Sample1', 'Sample2', 'Sample3', 'Sample4', 'Sample5'), each = 2),
  ValueType = rep(c('Type1', 'Type2'), times = 5),
  Value = c(5, 6, 3, 4, 6, 7, 4, 5, 7, 8)
)
# Create the side-by-side boxplot
box_plot_side_by_side <- ggplot(pcr_data, aes(x = SampleName, y = Value, fill = ValueType)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  labs(
    title = "Boxplots of Different Value Types by Sample Name",
    x = "Sample Name",
    y = "Value"
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),  # White background
    plot.background = element_rect(fill = "white", colour = "white"),   # White background
    axis.line = element_line(colour = "black"),  # Keep axis lines
    axis.ticks = element_line(colour = "black"),  # Keep axis ticks
    plot.title = element_text(hjust = 0.5)  # Center the title
  )

# Print the boxplot
print(box_plot_side_by_side)

# Save the plot
ggsave("box_plot_side_by_side.png", plot = box_plot_side_by_side, width = 8, height = 6)