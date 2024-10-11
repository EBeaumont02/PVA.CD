#### GRAPHS

source ("CD_PVA_MATRIX.R")

# proportion mature:
library(ggplot2)

# Create a data frame for the observed and predicted data
observed_data <- data.frame(Age = f.age.maturity, ProportionMature = f.prop.mature)
predicted_data <- data.frame(Age = age.vec.cont, ProportionMature = pred.p.matt)

library(ggplot2)

# Create a data frame for the observed and predicted data
observed_data <- data.frame(Age = f.age.maturity, ProportionMature = f.prop.mature)
predicted_data <- data.frame(Age = age.vec.cont, ProportionMature = pred.p.matt)

# Create the ggplot
ggplot() + 
  # Plot observed data points with larger size
  geom_point(data = observed_data, aes(x = Age, y = ProportionMature), color = "black", shape = 19, size = 3) +
  
  # Plot predicted line
  geom_line(data = predicted_data, aes(x = Age, y = ProportionMature), linetype = "dashed", size = 1.2, color = "black") +
  
  # Customize axis labels
  labs(x = "age (yrs)", y = "proportion mature") +
  
  # Set axis limits to ensure both axes start at 0, and adjust x and y limits to prevent data from being cut off
  scale_x_continuous(expand = c(0, 0), limits = c(0, max(f.age.maturity) + 1), breaks = seq(0, max(f.age.maturity) + 1, by = 2)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.05), breaks = seq(0, 1, by = 0.1)) +  # Increase upper y-limit to 1.05 to avoid cutting data off
  
  # Remove gridlines and set a minimal theme
  theme_minimal() +
  theme(
    panel.grid = element_blank(),                # Remove gridlines
    axis.line = element_line(color = "black"),   # Black axis lines
    axis.ticks = element_line(color = "black"),  # Black tick marks
    axis.text = element_text(size = 12, color = "black"),         # Adjust tick label size and ensure black color
    axis.title = element_text(size = 14, color = "black"),        # Adjust axis title size and ensure black color
    plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm")  # Add some margin around the plot
  )

# litter size:

library(ggplot2)

# Create a data frame for the observed and predicted data
observed_data <- size.litter  # The observed data (age and litter_size)
predicted_data <- data.frame(Age = litt.age.vec.cont, LitterSize = litt.pred.p.matt)

# Determine reasonable y-axis limits (based on observed and predicted data)
y_max <- max(c(observed_data$litter_size, predicted_data$LitterSize)) * 1.05  # Add a 5% buffer to avoid cutting off

# Create the ggplot
ggplot() + 
  # Plot observed data points with larger size
  geom_point(data = observed_data, aes(x = age, y = litter_size), color = "black", shape = 19, size = 3) +
  
  # Plot predicted line with specified preferences (black dashed line)
  geom_line(data = predicted_data, aes(x = Age, y = LitterSize), linetype = "dashed", size = 1.2, color = "black") +
  
  # Customize axis labels
  labs(x = "age (yrs)", y = "litter size") +
  
  # Set axis limits to ensure both axes start at 0, and add space to prevent cutting off data
  scale_x_continuous(expand = c(0, 0), limits = c(0, max(observed_data$age) + 1), breaks = seq(0, max(observed_data$age) + 1, by = 2)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, y_max), breaks = seq(0, y_max, by = 0.1)) +
  
  # Remove gridlines and set a minimal theme
  theme_minimal() +
  theme(
    panel.grid = element_blank(),                # Remove gridlines
    axis.line = element_line(color = "black"),   # Black axis lines
    axis.ticks = element_line(color = "black"),  # Black tick marks
    axis.text = element_text(size = 12, color = "black"),         # Adjust tick label size and ensure black color
    axis.title = element_text(size = 14, color = "black"),        # Adjust axis title size and ensure black color
    plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm")  # Add some margin around the plot
  )

# fertility:
# Create a data frame for the fertility data
data_fert <- data.frame(Age = seq_along(f.fert.vec), Fertility = f.fert.vec)

# Add Age = 0 and Fertility = 0
data_fert <- rbind(data.frame(Age = 0, Fertility = 0), data_fert)

# Create the ggplot
ggplot() + 
  # Plot observed fertility points
  #geom_point(data = data_fert, aes(x = Age, y = Fertility), color = "black", shape = 16, size = 3) +
  
  # Plot the fertility line, removing the last point and making it dashed
  geom_line(data = data_fert, aes(x = Age, y = Fertility), size = 1.4, color = "#70AD47") +
  
  # Customize axis labels
  labs(x = "age (yrs)", y = "fertility") +
  
  # Set axis limits to ensure both axes start at 0, add space at the top to prevent line cutoff
  scale_x_continuous(expand = c(0.03, 0.03), limits = c(0, length(f.fert.vec) + 1), breaks = seq(0, length(f.fert.vec) + 1, by = 2)) +
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, max(f.fert.vec) * 1.05), breaks = seq(0, max(f.fert.vec), by = 0.05)) +
  
  # Minimal theme with customizations to ensure axes and ticks are visible
  theme_minimal() +
  theme(
    panel.grid = element_blank(),                # Remove gridlines
    axis.line.x = element_line(color = "black", linewidth = 0.8), # Black x-axis line
    axis.line.y = element_line(color = "black", linewidth = 0.8), # Black y-axis line
    axis.ticks.x = element_line(color = "black", linewidth = 0.8), # Ensure ticks on x-axis
    axis.ticks.y = element_line(color = "black", linewidth = 0.8), # Ensure ticks on y-axis
    axis.text.x = element_text(size = 14, color = "black"), # X-axis text size and color
    axis.text.y = element_text(size = 14, color = "black"), # Y-axis text size and color
    axis.title.x = element_text(size = 16, color = "black"), # X-axis title size and color
    axis.title.y = element_text(size = 16, color = "black"), # Y-axis title size and color
    plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm")  # Add some margin around the plot
  )


# survival:
library(ggplot2)

# Create a data frame for the age and survival probability
data <- data.frame(Age = age.vec, Survival = surv.vec)


ggplot() + 
  # Plot the last observed data point
  geom_point(data = data[nrow(data), ], aes(x = Age, y = Survival), color = "#FF8E00", shape = 19, size = 3) +
  
  # Plot the line, removing the last point and making it dashed
  geom_line(data = data[-nrow(data), ], aes(x = Age, y = Survival), size = 1.5, color = "#FF8E00") +
  
  # Customize axis labels
  labs(x = "age (yrs)", y = "survival probability") +
  
  # Set axis limits to ensure both axes start at 0, and ensure no data is cut off on the right side
  scale_x_continuous(expand = c(0.01, 0.01), limits = c(0, max(age.vec) + 1), breaks = seq(0, max(age.vec) + 1, by = 2)) +
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  
  # Remove gridlines and set a minimal theme
  theme_minimal() +
  theme(
    panel.grid = element_blank(),                # Remove gridlines
    axis.line = element_line(color = "black", linewidth = 0.8),   # Black axis lines
    axis.ticks = element_line(color = "black", linewidth = 0.8),  # Black tick marks
    axis.text = element_text(size = 14, color = "black"),         # Adjust tick label size and ensure black color
    axis.title = element_text(size = 16, color = "black"),        # Adjust axis title size and ensure black color
    plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm")  # Add some margin around the plot
  )

# population age distribution:

# Create a data frame from the vectors
data.age.dist <- data.frame(Age = age.vec, Population = init.vec.total.pop)



# Create the ggplot
ggplot(data.age.dist, aes(x = Age, y = Population)) + 
  geom_line(color = "black") +                 # Line plot with black color
  labs(x = "age (yrs)", y = "number of individuals (N)") +  # Axis labels
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, max(age.vec), by = 2)) +  # More x-axis tick marks, start at 0
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, max(init.vec.total.pop), by = 10)) +  # More y-axis tick marks, start at 0
  theme_minimal() +                            # Clean minimal theme
  theme(
    panel.grid = element_blank(),              # Remove grid lines
    axis.line = element_line(color = "black"), # Add axis lines in black
    axis.ticks = element_line(color = "black"),# Add ticks in black
    axis.text = element_text(size = 12, color = "black"),  # Increase tick mark label font size and color
    axis.title = element_text(size = 14, color = "black")  # Increase axis title font size and color
  )


# harvesting proabbility: 
library(ggplot2)

# Create data frames for both the population and probability of becoming bycatch
data.age.dist <- data.frame(Age = age.vec, Population = init.vec.total.pop)
data.prob <- data.frame(Age = age.vec, Probability = prob.vec)

# Create the ggplot for population with a secondary axis for bycatch probability
ggplot(data.age.dist, aes(x = Age)) + 
  # Line plot for population with thicker line
  geom_line(aes(y = Population), color = "black", size = 1.5) +
  
  # Line plot for bycatch probability with thicker, dashed line
  geom_line(data = data.prob, aes(y = Probability * max(init.vec.total.pop) / max(prob.vec)), color = "#1E90FF", linetype = "dashed", size = 1.5) +
  
  # Primary axis labels for Population
  labs(x = "age (yrs)", y = "number of individuals (N)") +
  
  # Add a primary y-axis with a tick at 100, and a secondary y-axis for bycatch probability
  scale_y_continuous(
    name = "number of individuals (N)",
    breaks = c(seq(0, max(init.vec.total.pop), by = 10), 100),  # Ensure tick at 100
    sec.axis = sec_axis(~ . * max(prob.vec) / max(init.vec.total.pop), name = "probability of becoming bycatch")
  ) +
  
  # Set up x-axis and ticks
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, max(age.vec), by = 2)) +
  
  # Remove gridlines and use a clean minimal theme
  theme_minimal() +
  theme(
    panel.grid = element_blank(),                # Remove gridlines
    axis.line.x = element_line(color = "black"), # Black x-axis line
    axis.line.y = element_line(color = "black"), # Black y-axis lines
    axis.ticks = element_line(color = "black"),  # Black tick marks
    axis.text = element_text(size = 12, color = "black"), # Axis text size and color
    axis.title = element_text(size = 14, color = "black"), # Axis title size and color
    axis.text.y.right = element_text(color = "#1E90FF"),       # Right y-axis labels in red
    axis.title.y.right = element_text(color = "#1E90FF")       # Right y-axis title in red
  )




#Constant

library(ggplot2)

constant_bycatch_r_data <- read.table("constant_bycatch_r_data.csv",header=T,sep=",") 

# remove 0 to na
constant_bycatch_r_data$r.mn.harv[constant_bycatch_r_data$r.mn.harv == 0] <- NA
constant_bycatch_r_data$r.mn.lo[constant_bycatch_r_data$r.mn.lo == 0] <- NA
constant_bycatch_r_data$r.mn.up[constant_bycatch_r_data$r.mn.up == 0] <- NA

constant_bycatch_r_plot <- ggplot(constant_bycatch_r_data, aes(x = harv.p.vec)) +
  geom_ribbon(aes(ymin = r.mn.lo, ymax = r.mn.up), fill = "#ADD8E6", alpha = 0.8) +  # Light blue shaded area for confidence interval
  geom_line(aes(y = r.mn.harv), color = "black", size = 0.8, na.rm = TRUE) +  # Main line in black, on top of the shaded area
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.7) +  # Red dashed line at y = 0
  labs(x = "bycatch increase", y = "growth rate (r)") +  
  theme_classic() +  
  scale_y_continuous(limits = c(min(constant_bycatch_r_data$r.mn.harv, constant_bycatch_r_data$r.mn.lo, constant_bycatch_r_data$r.mn.up, na.rm = TRUE), 0),
                     breaks = pretty(c(min(constant_bycatch_r_data$r.mn.harv, constant_bycatch_r_data$r.mn.lo, constant_bycatch_r_data$r.mn.up, na.rm = TRUE), 0), n = 12),
                     expand = expansion(mult = c(0.05, 0.05))) +  
  scale_x_continuous(limits = c(0, 0.12),  # Change x-axis limits to 0.12
                     expand = expansion(mult = c(0, 0.05))) +  
  theme(axis.text = element_text(color = "black", size = 12),  
        axis.title = element_text(color = "black", size = 14),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

# Display the plot
print(constant_bycatch_r_plot)


constant_bycatch_n_data <- read.table("constant_bycatch_n_data.csv",header=T,sep=",") 

constant_bycatch_n_plot <- ggplot(constant_bycatch_n_data, aes(x = harv.p.vec)) +
  geom_ribbon(aes(ymin = n.mn.lo, ymax = n.mn.up), fill = "#ADD8E6", alpha = 0.8) +  # Light blue shaded area for confidence interval
  geom_line(aes(y = n.mn.harv), color = "black", size = 0.8) +  # Main line in black, on top of the shaded area
  labs(x = "bycatch increase", y = "population size (N)") +  
  theme_classic() +  
  scale_y_continuous(limits = c(0, max(constant_bycatch_n_data$n.mn.harv, constant_bycatch_n_data$n.mn.lo, constant_bycatch_n_data$n.mn.up, na.rm = TRUE)),
                     breaks = pretty(c(0, max(constant_bycatch_n_data$n.mn.harv, constant_bycatch_n_data$n.mn.lo, constant_bycatch_n_data$n.mn.up, na.rm = TRUE)), n = 12),
                     expand = expansion(mult = c(0, 0.1))) +  
  scale_x_continuous(limits = c(0, 0.12),  # Match x-axis limits to previous graph
                     expand = expansion(mult = c(0, 0.05))) +  
  theme(axis.text = element_text(color = "black", size = 12),  
        axis.title = element_text(color = "black", size = 14),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

# Display the plot
print(constant_bycatch_n_plot)


#Proportional

proportional_bycatch_r_data <- read.table("proportional_bycatch_r_data.csv",header=T,sep=",") 


# remove 0 to na
proportional_bycatch_r_data$r.mn.harv[proportional_bycatch_r_data$r.mn.harv == 0] <- NA
proportional_bycatch_r_data$r.mn.lo[proportional_bycatch_r_data$r.mn.lo == 0] <- NA
proportional_bycatch_r_data$r.mn.up[proportional_bycatch_r_data$r.mn.up == 0] <- NA

proportional_bycatch_r_plot <- ggplot(proportional_bycatch_r_data, aes(x = harv.p.vec)) +
  geom_ribbon(aes(ymin = r.mn.lo, ymax = r.mn.up), fill = "#ADD8E6", alpha = 10) +  # Light blue shaded area for confidence interval
  geom_line(aes(y = r.mn.harv), color = "black", size = 0.6, na.rm = TRUE) +  # Main line in black, on top of the shaded area
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.7) +  # Red dashed line at y = 0
  labs(x = "bycatch increase", y = "growth rate (r)") +  
  theme_classic() +  
  scale_y_continuous(limits = c(min(proportional_bycatch_r_data$r.mn.harv, proportional_bycatch_r_data$r.mn.lo, proportional_bycatch_r_data$r.mn.up, na.rm = TRUE), 0),
                     breaks = pretty(c(min(proportional_bycatch_r_data$r.mn.harv, proportional_bycatch_r_data$r.mn.lo, proportional_bycatch_r_data$r.mn.up, na.rm = TRUE), 0), n = 12),
                     expand = expansion(mult = c(0.05, 0.05))) +  
  scale_x_continuous(limits = c(0, 1),  # Set x-axis limits to 1
                     expand = expansion(mult = c(0, 0.05))) +  
  theme(axis.text = element_text(color = "black", size = 13),  
        axis.title = element_text(color = "black", size = 12),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

# Display the plot
print(proportional_bycatch_r_plot)


proportional_bycatch_n_data <- read.table("proportional_bycatch_n_data.csv",header=T,sep=",") 

proportional_bycatch_n_plot <- ggplot(proportional_bycatch_n_data, aes(x = harv.p.vec)) +
  geom_ribbon(aes(ymin = n.mn.lo, ymax = n.mn.up), fill = "#ADD8E6", alpha = 0.8) +  # Light blue shaded area for confidence interval
  geom_line(aes(y = n.mn.harv), color = "black", size = 0.8) +  # Main line in black, on top of the shaded area
  labs(x = "bycatch increase", y = "population size (N)") +  
  theme_classic() +  
  scale_y_continuous(limits = c(0, max(proportional_bycatch_n_data$n.mn.harv, proportional_bycatch_n_data$n.mn.lo, proportional_bycatch_n_data$n.mn.up, na.rm = TRUE)),
                     breaks = pretty(c(0, max(proportional_bycatch_n_data$n.mn.harv, proportional_bycatch_n_data$n.mn.lo, proportional_bycatch_n_data$n.mn.up, na.rm = TRUE)), n = 12),
                     expand = expansion(mult = c(0, 0.1))) +  
  scale_x_continuous(limits = c(0, 1),  # Match x-axis limits to previous graph
                     expand = expansion(mult = c(0, 0.05))) +  
  theme(axis.text = element_text(color = "black", size = 11),  
        axis.title = element_text(color = "black", size = 12),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

# Display the plot
print(proportional_bycatch_n_plot)


library(ggplot2)
library(patchwork)

# Constant Bycatch Growth Rate Plot (Graph C)
constant_bycatch_r_plot <- ggplot(constant_bycatch_r_data, aes(x = harv.p.vec)) +
  geom_ribbon(aes(ymin = r.mn.lo, ymax = r.mn.up), fill = "#ADD8E6", alpha = 0.8) +  # Light blue shaded area
  geom_line(aes(y = r.mn.harv), color = "black", size = 0.8, na.rm = TRUE) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.7) +
  labs(y = "growth rate (r)") +  # Removed x-axis label
  theme_classic() +  
  scale_y_continuous(limits = c(min(constant_bycatch_r_data$r.mn.harv, constant_bycatch_r_data$r.mn.lo, constant_bycatch_r_data$r.mn.up, na.rm = TRUE), 0),
                     breaks = pretty(c(min(constant_bycatch_r_data$r.mn.harv, constant_bycatch_r_data$r.mn.lo, constant_bycatch_r_data$r.mn.up, na.rm = TRUE), 0), n = 12),
                     expand = expansion(mult = c(0.05, 0.05))) +
  scale_x_continuous(limits = c(0, 0.12), breaks = seq(0, 0.12, by = 0.02),  # Set breaks and limits to match graph A
                     expand = expansion(mult = c(0, 0.05))) +  
  theme(axis.text = element_text(color = "black", size = 12),  
        axis.title.x = element_blank(),  # Removed x-axis label
        axis.title = element_text(color = "black", size = 14))

# Constant Bycatch Population Size Plot (Graph A)
constant_bycatch_n_plot <- ggplot(constant_bycatch_n_data, aes(x = harv.p.vec)) +
  geom_ribbon(aes(ymin = n.mn.lo, ymax = n.mn.up), fill = "#ADD8E6", alpha = 0.8) +  # Light blue shaded area
  geom_line(aes(y = n.mn.harv), color = "black", size = 0.8) +
  labs(x = NULL, y = "population size (N)") +  
  theme_classic() +  
  scale_y_continuous(limits = c(0, 900),
                     breaks = seq(0, 900, by = 100),  
                     expand = expansion(mult = c(0, 0.1))) +
  scale_x_continuous(limits = c(0, 0.12), breaks = seq(0, 0.12, by = 0.02),  # Ensure alignment with Graph C
                     expand = expansion(mult = c(0, 0.05))) +  
  theme(axis.text.x = element_blank(),  # Remove x-axis ticks and labels for Graph A
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),  # Removed x-axis label
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 14))

# Proportional Bycatch Growth Rate Plot (Graph D)
proportional_bycatch_r_plot <- ggplot(proportional_bycatch_r_data, aes(x = harv.p.vec)) +
  geom_ribbon(aes(ymin = r.mn.lo, ymax = r.mn.up), fill = "#ADD8E6", alpha = 1) +  # Light blue shaded area
  geom_line(aes(y = r.mn.harv), color = "black", size = 0.8, na.rm = TRUE) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.7) +
  labs(y = NULL) +  # Removed x-axis label
  theme_classic() +  
  scale_y_continuous(limits = c(min(proportional_bycatch_r_data$r.mn.harv, proportional_bycatch_r_data$r.mn.lo, proportional_bycatch_r_data$r.mn.up, na.rm = TRUE), 0),
                     breaks = pretty(c(min(proportional_bycatch_r_data$r.mn.harv, proportional_bycatch_r_data$r.mn.lo, proportional_bycatch_r_data$r.mn.up, na.rm = TRUE), 0), n = 12),
                     expand = expansion(mult = c(0.05, 0.05))) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +  
  theme(axis.text = element_text(color = "black", size = 12),  
        axis.title.x = element_blank(),  # Removed x-axis label
        axis.title = element_text(color = "black", size = 14),
        axis.text.y = element_blank(),  
        axis.ticks.y = element_blank())

# Proportional Bycatch Population Size Plot (Graph B)
proportional_bycatch_n_plot <- ggplot(proportional_bycatch_n_data, aes(x = harv.p.vec)) +
  geom_ribbon(aes(ymin = n.mn.lo, ymax = n.mn.up), fill = "#ADD8E6", alpha = 0.8) +  # Light blue shaded area
  geom_line(aes(y = n.mn.harv), color = "black", size = 0.8) +
  labs(x = NULL, y = NULL) +  
  theme_classic() +  
  scale_y_continuous(limits = c(0, max(proportional_bycatch_n_data$n.mn.harv, proportional_bycatch_n_data$n.mn.lo, proportional_bycatch_n_data$n.mn.up, na.rm = TRUE)),
                     breaks = seq(0, max(proportional_bycatch_n_data$n.mn.harv, na.rm = TRUE), by = 100),  
                     expand = expansion(mult = c(0, 0.1))) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +  
  theme(axis.text = element_blank(),  
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),  
        axis.title.x = element_blank(),  # Removed x-axis label
        axis.title = element_text(color = "black", size = 14))

# Combine the plots using patchwork with a single x-axis label in the middle
final_plot <- (constant_bycatch_n_plot + proportional_bycatch_n_plot) / 
  (constant_bycatch_r_plot + proportional_bycatch_r_plot) +
  plot_layout(guides = "collect") +  
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")) + 
  labs(x = "bycatch increase")  # Single x-axis label in the middle

# Display the combined plot
print(final_plot)

