############################################
## Common dolphin population projection model
############################################
## Harvest Rate Graphs

# Load necessary library
library(ggplot2)


# CONSTANT HARVEST GRAPHS
##Population Size 

# Read the data
population_size_harvest0 <- read.table("population_size_harvest0.csv", header = TRUE, sep = ",")
population_size_harvest0.05 <- read.table("population_size_harvest0.05.csv", header = TRUE, sep = ",")
population_size_harvest0.25 <- read.table("population_size_harvest0.25.csv", header = TRUE, sep = ",")
population_size_harvest0.5 <- read.table("population_size_harvest0.5.csv", header = TRUE, sep = ",")
population_size_harvest0.75 <- read.table("population_size_harvest0.75.csv", header = TRUE, sep = ",")
population_size_harvest0.95 <- read.table("population_size_harvest0.95.csv", header = TRUE, sep = ",")
population_size_harvest1 <- read.table("population_size_harvest1.csv", header = TRUE, sep = ",")

# Combine all data into one data frame with an additional column to identify the harvest level
combined_data <- rbind(
  data.frame(yplot.n = population_size_harvest0$yplot.n, n.mn = population_size_harvest0$n.mn, Harvest = "0"),
  data.frame(yplot.n = population_size_harvest0.05$yplot.n, n.mn = population_size_harvest0.05$n.mn, Harvest = "0.05"),
  data.frame(yplot.n = population_size_harvest0.25$yplot.n, n.mn = population_size_harvest0.25$n.mn, Harvest = "0.25"),
  data.frame(yplot.n = population_size_harvest0.5$yplot.n, n.mn = population_size_harvest0.5$n.mn, Harvest = "0.5"),
  data.frame(yplot.n = population_size_harvest0.75$yplot.n, n.mn = population_size_harvest0.75$n.mn, Harvest = "0.75"),
  data.frame(yplot.n = population_size_harvest0.95$yplot.n, n.mn = population_size_harvest0.95$n.mn, Harvest = "0.95"),
  data.frame(yplot.n = population_size_harvest1$yplot.n, n.mn = population_size_harvest1$n.mn, Harvest = "1")
)

# Plot the data
ggplot(combined_data, aes(x = yplot.n, y = n.mn, color = Harvest)) +
  geom_line(size = 0.9) +
  labs(x = "Year",
       y = "Mean Population Size",
       color = "Bycatch Level") +  # This sets the legend title
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),   # Remove grid
    axis.line = element_line(color = "black"), # Add axis lines
    plot.title = element_blank()    # Remove title
  )

# Save the plot with increased dimensions
ggsave("harvest_rate_lines_effect_on_population_size(constant harvest).jpeg", width = 10, height = 6)


##Growth Rate

# Read the data
population_growth_rate_harvest0 <- read.table("population_growth_rate_harvest0.csv", header = TRUE, sep = ",")
population_growth_rate_harvest0.05 <- read.table("population_growth_rate_harvest0.05.csv", header = TRUE, sep = ",")
population_growth_rate_harvest0.25 <- read.table("population_growth_rate_harvest0.25.csv", header = TRUE, sep = ",")
population_growth_rate_harvest0.5 <- read.table("population_growth_rate_harvest0.5.csv", header = TRUE, sep = ",")
population_growth_rate_harvest0.75 <- read.table("population_growth_rate_harvest0.75.csv", header = TRUE, sep = ",")
population_growth_rate_harvest0.95 <- read.table("population_growth_rate_harvest0.95.csv", header = TRUE, sep = ",")
population_growth_rate_harvest1 <- read.table("population_growth_rate_harvest1.csv", header = TRUE, sep = ",")

# Replace NA and NaN values with 0
population_growth_rate_harvest0[is.na(population_growth_rate_harvest0)] <- 0
population_growth_rate_harvest0.05[is.na(population_growth_rate_harvest0.05)] <- 0
population_growth_rate_harvest0.25[is.na(population_growth_rate_harvest0.25)] <- 0
population_growth_rate_harvest0.5[is.na(population_growth_rate_harvest0.5)] <- 0
population_growth_rate_harvest0.75[is.na(population_growth_rate_harvest0.75)] <- 0
population_growth_rate_harvest0.95[is.na(population_growth_rate_harvest0.95)] <- 0
population_growth_rate_harvest1[is.na(population_growth_rate_harvest1)] <- 0

# Combine all data into one data frame with an additional column to identify the harvest level
combined_growth_rate_data <- rbind(
  data.frame(yplot.r = population_growth_rate_harvest0$yplot.r, r.mn = population_growth_rate_harvest0$r.mn, Harvest = "0"),
  data.frame(yplot.r = population_growth_rate_harvest0.05$yplot.r, r.mn = population_growth_rate_harvest0.05$r.mn, Harvest = "0.05"),
  data.frame(yplot.r = population_growth_rate_harvest0.25$yplot.r, r.mn = population_growth_rate_harvest0.25$r.mn, Harvest = "0.25"),
  data.frame(yplot.r = population_growth_rate_harvest0.5$yplot.r, r.mn = population_growth_rate_harvest0.5$r.mn, Harvest = "0.5"),
  data.frame(yplot.r = population_growth_rate_harvest0.75$yplot.r, r.mn = population_growth_rate_harvest0.75$r.mn, Harvest = "0.75"),
  data.frame(yplot.r = population_growth_rate_harvest0.95$yplot.r, r.mn = population_growth_rate_harvest0.95$r.mn, Harvest = "0.95"),
  data.frame(yplot.r = population_growth_rate_harvest1$yplot.r, r.mn = population_growth_rate_harvest1$r.mn, Harvest = "1")
)

# Plot the data with expanded y-axis limits
ggplot(combined_growth_rate_data, aes(x = yplot.r, y = r.mn, color = Harvest)) +
  geom_line(size = 0.9) +
  labs(x = "Year",
       y = "Mean Population Growth Rate",
       color = "Bycatch Level") +  # This sets the legend title
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),   # Remove grid
    axis.line = element_line(color = "black"), # Add axis lines
    plot.title = element_blank()    # Remove title
  ) +
  expand_limits(y = c(min(combined_growth_rate_data$r.mn) - 0.1, 
                      max(combined_growth_rate_data$r.mn) + 0.1))  # Expand y-axis limits

# Save the plot with increased dimensions
ggsave("harvest_rate_lines_effect_on_growth_rate(constant harvest).jpeg", width = 10, height = 6)




# PROPORTIONAL HARVEST GRAPHS
##Population Size 

# Read the data
population_size_harvest0_p <- read.table("population_size_harvest0_p.csv", header = TRUE, sep = ",")
population_size_harvest0.05_p <- read.table("population_size_harvest0.05_p.csv", header = TRUE, sep = ",")
population_size_harvest0.25_p <- read.table("population_size_harvest0.25_p.csv", header = TRUE, sep = ",")
population_size_harvest0.5_p <- read.table("population_size_harvest0.5_p.csv", header = TRUE, sep = ",")
population_size_harvest0.75_p <- read.table("population_size_harvest0.75_p.csv", header = TRUE, sep = ",")
population_size_harvest0.95_p <- read.table("population_size_harvest0.95_p.csv", header = TRUE, sep = ",")
population_size_harvest1_p <- read.table("population_size_harvest1_p.csv", header = TRUE, sep = ",")

# Combine all data into one data frame with an additional column to identify the harvest level
combined_data_p <- rbind(
  data.frame(yplot.n = population_size_harvest0_p$yplot.n, n.mn = population_size_harvest0_p$n.mn, Harvest = "0"),
  data.frame(yplot.n = population_size_harvest0.05_p$yplot.n, n.mn = population_size_harvest0.05_p$n.mn, Harvest = "0.05"),
  data.frame(yplot.n = population_size_harvest0.25_p$yplot.n, n.mn = population_size_harvest0.25_p$n.mn, Harvest = "0.25"),
  data.frame(yplot.n = population_size_harvest0.5_p$yplot.n, n.mn = population_size_harvest0.5_p$n.mn, Harvest = "0.5"),
  data.frame(yplot.n = population_size_harvest0.75_p$yplot.n, n.mn = population_size_harvest0.75_p$n.mn, Harvest = "0.75"),
  data.frame(yplot.n = population_size_harvest0.95_p$yplot.n, n.mn = population_size_harvest0.95_p$n.mn, Harvest = "0.95"),
  data.frame(yplot.n = population_size_harvest1_p$yplot.n, n.mn = population_size_harvest1_p$n.mn, Harvest = "1")
)

# Plot the data
ggplot(combined_data_p, aes(x = yplot.n, y = n.mn, color = Harvest)) +
  geom_line(size = 0.9) +
  labs(x = "Year",
       y = "Mean Population Size",
       color = "Bycatch Level") +  # This sets the legend title
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),   # Remove grid
    axis.line = element_line(color = "black"), # Add axis lines
    plot.title = element_blank()    # Remove title
  )

# Save the plot with increased dimensions
ggsave("harvest_rate_lines_effect_on_population_size(proportional harvest).jpeg", width = 10, height = 6)


##Growth Rate

# Read the data
population_growth_rate_harvest0_p <- read.table("population_growth_rate_harvest0_p.csv", header = TRUE, sep = ",")
population_growth_rate_harvest0.05_p <- read.table("population_growth_rate_harvest0.05_p.csv", header = TRUE, sep = ",")
population_growth_rate_harvest0.25_p <- read.table("population_growth_rate_harvest0.25_p.csv", header = TRUE, sep = ",")
population_growth_rate_harvest0.5_p <- read.table("population_growth_rate_harvest0.5_p.csv", header = TRUE, sep = ",")
population_growth_rate_harvest0.75_p <- read.table("population_growth_rate_harvest0.75_p.csv", header = TRUE, sep = ",")
population_growth_rate_harvest0.95_p <- read.table("population_growth_rate_harvest0.95_p.csv", header = TRUE, sep = ",")
population_growth_rate_harvest1_p <- read.table("population_growth_rate_harvest1_p.csv", header = TRUE, sep = ",")

# Combine all data into one data frame with an additional column to identify the harvest level
combined_growth_rate_data_p <- rbind(
  data.frame(yplot.r = population_growth_rate_harvest0_p$yplot.r, r.mn = population_growth_rate_harvest0_p$r.mn, Harvest = "0"),
  data.frame(yplot.r = population_growth_rate_harvest0.05_p$yplot.r, r.mn = population_growth_rate_harvest0.05_p$r.mn, Harvest = "0.05"),
  data.frame(yplot.r = population_growth_rate_harvest0.25_p$yplot.r, r.mn = population_growth_rate_harvest0.25_p$r.mn, Harvest = "0.25"),
  data.frame(yplot.r = population_growth_rate_harvest0.5_p$yplot.r, r.mn = population_growth_rate_harvest0.5_p$r.mn, Harvest = "0.5"),
  data.frame(yplot.r = population_growth_rate_harvest0.75_p$yplot.r, r.mn = population_growth_rate_harvest0.75_p$r.mn, Harvest = "0.75"),
  data.frame(yplot.r = population_growth_rate_harvest0.95_p$yplot.r, r.mn = population_growth_rate_harvest0.95_p$r.mn, Harvest = "0.95"),
  data.frame(yplot.r = population_growth_rate_harvest1_p$yplot.r, r.mn = population_growth_rate_harvest1_p$r.mn, Harvest = "1")
)

# Plot the data with expanded y-axis limits
ggplot(combined_growth_rate_data_p, aes(x = yplot.r, y = r.mn, color = Harvest)) +
  geom_line(size = 0.9) +
  labs(x = "Year",
       y = "Mean Population Growth Rate",
       color = "Bycatch Level") +  # This sets the legend title
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),   # Remove grid
    axis.line = element_line(color = "black"), # Add axis lines
    plot.title = element_blank()    # Remove title
  ) +
  expand_limits(y = c(min(combined_growth_rate_data_p$r.mn) - 0.1, 
                      max(combined_growth_rate_data_p$r.mn) + 0.1))  # Expand y-axis limits

# Save the plot with increased dimensions
ggsave("harvest_rate_lines_effect_on_growth_rate(proportional harvest).jpeg", width = 10, height = 6)

##saving to csv
yplot.r <- yrs[(gen.l + 1):(t)]
r.mn_harvest1_p <- data.frame(yplot.r, r.mn, r.lo, r.up)
write.csv(r.mn_harvest1_p, file = 'population_growth_rate_harvest1_p.csv', row.names = F )



# add on to growth rate code:


##saving to csv
yplot.n <- yrs[(gen.l + 1):(t + 1)]
n.mn_harvest1_p <- data.frame(yplot.n, n.mn, n.lo, n.up)
write.csv(n.mn_harvest1_p, file = 'population_size_harvest1_p.csv', row.names = F )


