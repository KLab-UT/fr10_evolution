# Loading in packages 
install.packages("ggplot2")
install.packages("reshape2")
install.packages("pheatmap")
library(ggplot2)
library(reshape2)
library(pheatmap)

# loading in the csv and viewing the csv
heatmap<- read.csv("rmsd_values.csv")
View(heatmap) 


data_matrix <- as.matrix(heatmap)  

# Reshape data for ggplot2
melted_data <- melt(heatmap)

# Create heat map using ggplot2
ggplot(data = heatmap, aes(x = P1, y = P2, fill = rmsd)) +
  geom_tile() +
  scale_fill_gradient(low = "Blue", high = "red") +
  theme_minimal() +
  labs(x = "Rows", y = "Columns", fill = "Value")


# Load the data
heatmap <- read.csv("rmsd_values.csv")





