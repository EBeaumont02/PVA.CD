############################################
## Common dolphin population projection model
############################################
## Bycatch Senario Graphs

# Load necessary library
library(ggplot2)


# CONSTANT HARVEST GRAPHS
##Population Size 

# Read the data
constant_bycatch_pop_size0 <- read.table("constant_bycatch_pop_size0.csv", header = TRUE, sep = ",")
constant_bycatch_pop_size0.25 <- read.table("constant_bycatch_pop_size0.25.csv", header = TRUE, sep = ",")
constant_bycatch_pop_size0.5 <- read.table("constant_bycatch_pop_size0.5.csv", header = TRUE, sep = ",")
constant_bycatch_pop_size0.75 <- read.table("constant_bycatch_pop_size0.75.csv", header = TRUE, sep = ",")
constant_bycatch_pop_size1 <- read.table("constant_bycatch_pop_size1.csv", header = TRUE, sep = ",")

# Combine all data into one data frame with an additional column to identify the bycatch level
combined_data_bycatch_cps <- rbind(
  data.frame(yplot.n = constant_bycatch_pop_size0$yplot.n, n.mn = constant_bycatch_pop_size0$n.mn, Bycatch = "0"),
  data.frame(yplot.n = constant_bycatch_pop_size0.25$yplot.n, n.mn = constant_bycatch_pop_size0.25$n.mn, Bycatch = "0.25"),
  data.frame(yplot.n = constant_bycatch_pop_size0.5$yplot.n, n.mn = constant_bycatch_pop_size0.5$n.mn, Bycatch = "0.5"),
  data.frame(yplot.n = constant_bycatch_pop_size0.75$yplot.n, n.mn = constant_bycatch_pop_size0.75$n.mn, Bycatch = "0.75"),
  data.frame(yplot.n = constant_bycatch_pop_size1$yplot.n, n.mn = constant_bycatch_pop_size1$n.mn, Bycatch = "1")
)


# Define your custom color palette
# Define custom color palette and labels for bycatch levels
custom_colors <- c("0" = "#585591", "0.25" = "#4689AA", "0.5" = "#5FA160", "0.75" = "#DAAE3A", "1" = "#C82B1F")
custom_labels <- c("0" = "0", "0.25" = "25", "0.5" = "50", "0.75" = "75", "1" = "100")

ggplot(combined_data_bycatch_cps, aes(x = yplot.n, y = n.mn, color = Bycatch)) +
  geom_line(size = 1.1) +
  labs(x = "year",
       y = "population size (N)",
       color = "Bycatch increase (%)") +  # This sets the legend title
  scale_color_manual(values = custom_colors, labels = custom_labels) +  # Apply custom colors and labels
  scale_x_continuous(limits = c(2037, 2100), breaks = seq(2040, 2100, by = 10),expand = c(0, 0)) +  # Remove gap between data and x-axis
  scale_y_continuous(limits = c(0, 1000), breaks = seq(0, 900, by = 100),expand = c(0, 0)) +  # Remove gap between data and y-axis
  theme_minimal(base_size = 12) +
  theme(
    axis.line = element_line(color = "black"),  # Add axis lines
    axis.ticks = element_line(color = "black"),  # Ensure tick marks are visible
    axis.ticks.length = unit(0.25, "cm"),  # Adjust tick mark length
    axis.text.x = element_text(size = 12, color = "black"),  # Increase x-axis tick mark label size
    axis.text.y = element_text(size = 12, color = "black"),  # Increase y-axis tick mark label size
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    plot.title = element_blank()  # Remove title
  )

# Save the plot with increased dimensions
ggsave("bycatch_rate_lines_effect_on_population_size(constant bycatch).jpeg", width = 10, height = 6)


##Growth Rate

# Read the data
constant_bycatch_growth_rate0 <- read.table("constant_bycatch_growth_rate0.csv", header = TRUE, sep = ",")
constant_bycatch_growth_rate0.25 <- read.table("constant_bycatch_growth_rate0.25.csv", header = TRUE, sep = ",")
constant_bycatch_growth_rate0.5 <- read.table("constant_bycatch_growth_rate0.5.csv", header = TRUE, sep = ",")
constant_bycatch_growth_rate0.75 <- read.table("constant_bycatch_growth_rate0.75.csv", header = TRUE, sep = ",")
constant_bycatch_growth_rate1 <- read.table("constant_bycatch_growth_rate1.csv", header = TRUE, sep = ",")

# Replace NA and NaN values with 0
constant_bycatch_growth_rate0[is.na(constant_bycatch_growth_rate0)] <- 0
constant_bycatch_growth_rate0.25[is.na(constant_bycatch_growth_rate0.25)] <- 0
constant_bycatch_growth_rate0.5[is.na(constant_bycatch_growth_rate0.5)] <- 0
constant_bycatch_growth_rate0.75[is.na(constant_bycatch_growth_rate0.75)] <- 0
constant_bycatch_growth_rate1[is.na(constant_bycatch_growth_rate1)] <- 0

# Combine all data into one data frame with an additional column to identify the harvest level
combined_data_bycatch_cgr <- rbind(
  data.frame(yplot.r = constant_bycatch_growth_rate0$yplot.r, r.mn = constant_bycatch_growth_rate0$r.mn, Bycatch = "0"),
  data.frame(yplot.r = constant_bycatch_growth_rate0.25$yplot.r, r.mn = constant_bycatch_growth_rate0.25$r.mn, Bycatch = "0.25"),
  data.frame(yplot.r = constant_bycatch_growth_rate0.5$yplot.r, r.mn = constant_bycatch_growth_rate0.5$r.mn, Bycatch = "0.5"),
  data.frame(yplot.r = constant_bycatch_growth_rate0.75$yplot.r, r.mn = constant_bycatch_growth_rate0.75$r.mn, Bycatch = "0.75"),
  data.frame(yplot.r = constant_bycatch_growth_rate1$yplot.r, r.mn = constant_bycatch_growth_rate1$r.mn, Bycatch = "1")
)

# Load the dplyr package for the pipe operator and filter function
library(dplyr)

# Define custom color palette and labels for bycatch levels
custom_colors <- c("0" = "#585591", "0.25" = "#4689AA", "0.5" = "#5FA160", "0.75" = "#DAAE3A", "1" = "#C82B1F")
custom_labels <- c("0" = "0", "0.25" = "25", "0.5" = "50", "0.75" = "75", "1" = "100")

# Filter out rows where growth rate (r.mn) is zero
filtered_data <- combined_data_bycatch_cgr %>%
  filter(r.mn != 0)

# Plot the filtered data
ggplot(filtered_data, aes(x = yplot.r, y = r.mn, color = Bycatch)) +
  geom_line(size = 1.1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.7) +  # Add dashed red line at y = 0
  labs(x = "year",
       y = "growth rate (r)",
       color = "Bycatch increase (%)") +  # This sets the legend title
  scale_color_manual(values = custom_colors, labels = custom_labels) +  # Apply custom colors and labels
  scale_x_continuous(limits = c(2037, 2100), breaks = seq(2040, 2100, by = 10),expand = c(0, 0)) +  # Remove gap between data and x-axis
  scale_y_continuous(limits = c(-1.3, 0.1), breaks = seq(-1.4, 0, by = 0.1), expand = c(0, 0.05)) +  # Add buffer and set y-axis increments
  theme_minimal(base_size = 12) +
  theme(
    axis.line = element_line(color = "black"),  # Add axis lines
    axis.ticks = element_line(color = "black"),  # Ensure tick marks are visible
    axis.ticks.length = unit(0.25, "cm"),  # Adjust tick mark length
    axis.text.x = element_text(size = 12, color = "black"),  # Increase x-axis tick mark label size
    axis.text.y = element_text(size = 12, color = "black"),  # Increase y-axis tick mark label size
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    plot.title = element_blank()  # Remove title
  )

# Save the plot with increased dimensions
ggsave("bycatch_rate_lines_effect_on_population_size(constant_bycatch_filtered).jpeg", width = 10, height = 6)

# PROPORTIONAL HARVEST GRAPHS
##Population Size 

# Read the data
prop_bycatch_pop_size0 <- read.table("prop_bycatch_pop_size0.csv", header = TRUE, sep = ",")
prop_bycatch_pop_size0.25 <- read.table("prop_bycatch_pop_size0.25.csv", header = TRUE, sep = ",")
prop_bycatch_pop_size0.5 <- read.table("prop_bycatch_pop_size0.5.csv", header = TRUE, sep = ",")
prop_bycatch_pop_size0.75 <- read.table("prop_bycatch_pop_size0.75.csv", header = TRUE, sep = ",")
prop_bycatch_pop_size1 <- read.table("prop_bycatch_pop_size1.csv", header = TRUE, sep = ",")

# Combine all data into one data frame with an additional column to identify the bycatch level
combined_data_bycatch_pps <- rbind(
  data.frame(yplot.n = prop_bycatch_pop_size0$yplot.n, n.mn = prop_bycatch_pop_size0$n.mn, Bycatch = "0"),
  data.frame(yplot.n = prop_bycatch_pop_size0.25$yplot.n, n.mn = prop_bycatch_pop_size0.25$n.mn, Bycatch = "0.25"),
  data.frame(yplot.n = prop_bycatch_pop_size0.5$yplot.n, n.mn = prop_bycatch_pop_size0.5$n.mn, Bycatch = "0.5"),
  data.frame(yplot.n = prop_bycatch_pop_size0.75$yplot.n, n.mn = prop_bycatch_pop_size0.75$n.mn, Bycatch = "0.75"),
  data.frame(yplot.n = prop_bycatch_pop_size1$yplot.n, n.mn = prop_bycatch_pop_size1$n.mn, Bycatch = "1")
)

# Define your custom color palette
# Define custom color palette and labels for bycatch levels
custom_colors <- c("0" = "#585591", "0.25" = "#4689AA", "0.5" = "#5FA160", "0.75" = "#DAAE3A", "1" = "#C82B1F")
custom_labels <- c("0" = "0", "0.25" = "25", "0.5" = "50", "0.75" = "75", "1" = "100")

ggplot(combined_data_bycatch_pps, aes(x = yplot.n, y = n.mn, color = Bycatch)) +
  geom_line(size = 1.1) +
  labs(x = "year",
       y = "population size (N)",
       color = "bycatch increase (%)") +  # This sets the legend title
  scale_color_manual(values = custom_colors, labels = custom_labels) +  # Apply custom colors and labels
  scale_x_continuous(limits = c(2037, 2100), breaks = seq(2040, 2100, by = 10),expand = c(0, 0)) +  # Remove gap between data and x-axis
  scale_y_continuous(limits = c(0, 1000), breaks = seq(0, 900, by = 100),expand = c(0, 0)) +  # Remove gap between data and y-axis
  theme_minimal(base_size = 12) +
  theme(
    axis.line = element_line(color = "black"),  # Add axis lines
    axis.ticks = element_line(color = "black"),  # Ensure tick marks are visible
    axis.ticks.length = unit(0.25, "cm"),  # Adjust tick mark length
    axis.text.x = element_text(size = 12, color = "black"),  # Increase x-axis tick mark label size
    axis.text.y = element_text(size = 12, color = "black"),  # Increase y-axis tick mark label size
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    plot.title = element_blank()  # Remove title
  )

# Save the plot with increased dimensions
ggsave("bycatch_rate_lines_effect_on_population_size(proportional bycatch).jpeg", width = 10, height = 6)



##Growth Rate

# Read the data
prop_bycatch_growth_rate0 <- read.table("prop_bycatch_growth_rate0.csv", header = TRUE, sep = ",")
prop_bycatch_growth_rate0.25 <- read.table("prop_bycatch_growth_rate0.25.csv", header = TRUE, sep = ",")
prop_bycatch_growth_rate0.5 <- read.table("prop_bycatch_growth_rate0.5.csv", header = TRUE, sep = ",")
prop_bycatch_growth_rate0.75 <- read.table("prop_bycatch_growth_rate0.75.csv", header = TRUE, sep = ",")
prop_bycatch_growth_rate1 <- read.table("prop_bycatch_growth_rate1.csv", header = TRUE, sep = ",")

# Combine all data into one data frame with an additional column to identify the bycatch level
combined_data_bycatch_pgr <- rbind(
  data.frame(yplot.r = prop_bycatch_growth_rate0$yplot.r, r.mn = prop_bycatch_growth_rate0$r.mn, Bycatch = "0"),
  data.frame(yplot.r = prop_bycatch_growth_rate0.25$yplot.r, r.mn = prop_bycatch_growth_rate0.25$r.mn, Bycatch = "0.25"),
  data.frame(yplot.r = prop_bycatch_growth_rate0.5$yplot.r, r.mn = prop_bycatch_growth_rate0.5$r.mn, Bycatch = "0.5"),
  data.frame(yplot.r = prop_bycatch_growth_rate0.75$yplot.r, r.mn = prop_bycatch_growth_rate0.75$r.mn, Bycatch = "0.75"),
  data.frame(yplot.r = prop_bycatch_growth_rate1$yplot.r, r.mn = prop_bycatch_growth_rate1$r.mn, Bycatch = "1")
)


# Define your custom color palette
# Define custom color palette and labels for bycatch levels
custom_colors <- c("0" = "#585591", "0.25" = "#4689AA", "0.5" = "#5FA160", "0.75" = "#DAAE3A", "1" = "#C82B1F")
custom_labels <- c("0" = "0", "0.25" = "25", "0.5" = "50", "0.75" = "75", "1" = "100")

ggplot(combined_data_bycatch_pgr, aes(x = yplot.r, y = r.mn, color = Bycatch)) +
  geom_line(size = 1.1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.7) +  # Add dashed red line at y = 0
  labs(x = "year",
       y = "growth rate (r)",
       color = "Bycatch increase (%)") +  # This sets the legend title
  scale_color_manual(values = custom_colors, labels = custom_labels) +  # Apply custom colors and labels
  scale_x_continuous(limits = c(2037, 2100), breaks = seq(2040, 2100, by = 10),expand = c(0, 0)) +  # Remove gap between data and x-axis
  scale_y_continuous(limits = c(-0.11, 0), breaks = seq(-0.11, 0, by = 0.01), expand = c(0, 0.01)) +  # Add buffer and set y-axis increments
  theme_minimal(base_size = 12) +
  theme(
    axis.line = element_line(color = "black"),  # Add axis lines
    axis.ticks = element_line(color = "black"),  # Ensure tick marks are visible
    axis.ticks.length = unit(0.25, "cm"),  # Adjust tick mark length
    axis.text.x = element_text(size = 12, color = "black"),  # Increase x-axis tick mark label size
    axis.text.y = element_text(size = 12, color = "black"),  # Increase y-axis tick mark label size
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    plot.title = element_blank()  # Remove title
  )

# Save the plot with increased dimensions
ggsave("bycatch_rate_lines_effect_on_growth_rate(proportional bycatch).jpeg", width = 10, height = 6)

# Load the patchwork library
library(patchwork)

# Save the combined plot
ggsave("aligned_bycatch_plots.jpeg", plot = final_combined_plot, width = 12, height = 10)




# CONSTANT HARVEST GRAPHS
## Population Size (Constant Bycatch) - Graph A
constant_pop_size_plot <- ggplot(combined_data_bycatch_cps, aes(x = yplot.n, y = n.mn, color = Bycatch)) +
  geom_line(size = 1.1) +
  labs(x = NULL,  # Remove x-axis label
       y = "population size (N)",  # Keep y-axis label for this plot
       color = "Bycatch increase (%)") +  # This sets the legend title
  scale_color_manual(values = custom_colors, labels = custom_labels) +
  scale_x_continuous(limits = c(2037, 2100), breaks = seq(2040, 2100, by = 10), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1000), breaks = seq(0, 900, by = 100), expand = c(0, 0)) +
  theme_minimal(base_size = 12) +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),  # Add y-axis ticks for A
    axis.text.y = element_text(size = 12, color = "black"),  # Add y-axis labels for A
    axis.title.y = element_text(size = 14),  # Increase y-axis label size
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    axis.text.x = element_blank(),  # Remove x-axis labels
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_blank()
  )

## Growth Rate (Constant Bycatch) - Graph C
constant_growth_rate_plot <- ggplot(filtered_data, aes(x = yplot.r, y = r.mn, color = Bycatch)) +
  geom_line(size = 1.1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.7) +  # Add dashed red line at y = 0
  labs(x = NULL, y = "growth rate (r)", color = "Bycatch increase (%)") +
  scale_color_manual(values = custom_colors, labels = custom_labels) +
  scale_x_continuous(limits = c(2037, 2100), breaks = seq(2040, 2100, by = 10), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1.3, 0), breaks = seq(-1.3, 0, by = 0.1), expand = c(0.05, 0.05)) +
  theme_minimal(base_size = 12) +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    axis.ticks.x = element_line(color = "black"),  # Keep x-axis ticks
    axis.ticks.y = element_line(color = "black"),  # Keep y-axis ticks
    axis.text.x = element_text(size = 12, color = "black"),  # Keep x-axis labels
    axis.text.y = element_text(size = 12, color = "black"),  # Keep y-axis labels
    axis.title.x = element_text(size = 14),  # Increase x-axis label size
    axis.title.y = element_text(size = 14),  # Increase y-axis label size
    plot.margin = unit(c(0, 0.5, 0, 0), "cm"),  # Adjust margin to match horizontally
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_blank()
  )

# PROPORTIONAL HARVEST GRAPHS
## Population Size (Proportional Bycatch) - Graph B
proportional_pop_size_plot <- ggplot(combined_data_bycatch_pps, aes(x = yplot.n, y = n.mn, color = Bycatch)) +
  geom_line(size = 1.1) +
  labs(x = NULL,  # Remove x-axis label
       y = NULL,  # Remove y-axis label for this plot
       color = "Bycatch increase (%)") +  # This sets the legend title
  scale_color_manual(values = custom_colors, labels = custom_labels) +
  scale_x_continuous(limits = c(2037, 2100), breaks = seq(2040, 2100, by = 10), expand = c(0, 0)) +  # Ensure alignment with Graph D
  scale_y_continuous(limits = c(0, 1000), breaks = seq(0, 900, by = 100), expand = c(0, 0)) +  # No expansion
  theme_minimal(base_size = 12) +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    axis.ticks.y = element_blank(),  # Remove y-axis ticks
    axis.text.x = element_blank(),  # Remove x-axis labels
    axis.text.y = element_blank(),  # Remove y-axis labels
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0.75), "cm")  # Reduced margin to adjust alignment
  )

## Growth Rate (Proportional Bycatch) - Graph D
proportional_growth_rate_plot <- ggplot(combined_data_bycatch_pgr, aes(x = yplot.r, y = r.mn_scaled, color = Bycatch)) +
  geom_line(size = 1.1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.7) +
  labs(x = NULL,  # Keep x-axis label
       y = NULL,  # Remove y-axis label
       color = "Bycatch increase (%)") +
  scale_color_manual(values = custom_colors, labels = custom_labels) +
  scale_x_continuous(limits = c(2037, 2100), breaks = seq(2040, 2100, by = 10), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1.3, 0), breaks = seq(-1.3, 0, by = 0.1), expand = c(0.05, 0.05)) +  # No buffer for alignment with C
  theme_minimal(base_size = 12) +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks.y = element_blank(),  # Remove y-axis ticks for D
    axis.text.y = element_blank(),  # Remove y-axis labels for D
    axis.ticks.x = element_line(color = "black"),  # Keep x-axis ticks
    axis.text.x = element_text(size = 12, color = "black"),  # Keep x-axis labels
    axis.title.x = element_text(size = 14),
    plot.margin = unit(c(0, 0, 0, 0.4), "cm"),  # Adjust left margin for alignment
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_blank()
  )

# Combine the plots into one figure using patchwork
final_combined_plot <- (constant_pop_size_plot + proportional_pop_size_plot) / 
  (constant_growth_rate_plot + proportional_growth_rate_plot) +
  plot_layout(guides = "collect", widths = unit(c(1, 1.05), "null"))  # Slightly reduce width of D

# Display the final combined plot
print(final_combined_plot)
