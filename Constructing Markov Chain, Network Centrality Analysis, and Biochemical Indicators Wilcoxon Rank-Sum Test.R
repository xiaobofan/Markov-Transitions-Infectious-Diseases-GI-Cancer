install.packages("tidyr")
install.packages("tidyverse")
install.packages("markovchain")
library(tidyr)
library(tidyverse)
library(markovchain)

### Construct the Markov chain ###
# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Set input and output file paths
# args[1] - Input file path
# args[2] - Output file path
input_file <- ifelse(length(args) >= 1, args[1], "test data mapping.txt")
output_file <- ifelse(length(args) >= 2, args[2], "Colorectal Cancer transitionMatrix.csv")

# Set the working directory.
setwd(dirname(input_file))

# Load the data.
SC_data <- read.csv(file = input_file, header = FALSE, fill = FALSE)

# Construct the Markov chain
data_colorectal <- SC_data %>%
  pull(V1) %>%
  strsplit(" ") %>%
  unlist()
fit_markov_colorectal <- markovchainFit(data_colorectal)

# Extract the transition matrix and save it to a CSV file.
data <- fit_markov_colorectal$estimate@transitionMatrix
write.csv(data, file = output_file)

### Network graph analysis ###
install.packages("igraph")
library(igraph)
# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Set input and output file paths
input_file <- ifelse(length(args) >= 1, args[1], "0.05.csv") #The 0.05.csv file contains specific transition data selected from the transition matrix obtained in the previous step, including only transitions with probabilities greater than 0.05 that are related to cancer.
output_file <- ifelse(length(args) >= 2, args[2], "Colorectal_cancer_network_centrality_results_without_Prion.csv")

# Set the working directory dynamically based on input file directory
setwd(dirname(input_file))

# Read in the input file
GD_file <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)
GD <- graph_from_data_frame(GD_file, directed = TRUE)

# Add edge weight as width
E(GD)$width <- E(GD)$weight

# Calculating node centrality, excluding self-loops, as they can affect the overall network centrality
de <- degree(GD, normalize = TRUE) # Vertex degree 
de_out <- degree(GD, mode = "out", normalize = TRUE) # Node out-degree
de_in <- degree(GD, mode = "in", normalized = TRUE) # Node in-degree
be <- betweenness(GD, normalized = TRUE) # Node Betweenness
cl <- closeness(GD, normalized = TRUE) # Node Closeness
cl_out <- closeness(GD, mode = "out", normalized = TRUE) # Node out-closeness
cl_in <- closeness(GD, mode = "in", normalized = TRUE) # Node in-closeness
ed <- eigen_centrality(GD) # Node Eigen-centrality
ed_sort <- sort(ed$vector, decreasing = TRUE) # Sort data in descending order

# Create a data frame with the centrality measures
GD_network <- data.frame(degree = de, 
                         degree_out = de_out, 
                         degree_in = de_in, 
                         betweenness = be, 
                         closeness = cl, 
                         closeness_out = cl_out, 
                         closeness_in = cl_in, 
                         eigen_centrality = ed$vector)

# Write the results to a CSV file
write.csv(GD_network, file = output_file)


### Biochemical indicators Wilcoxon rank-sum test analysis ###
install.packages("effectsize")
library(effectsize)

x <- c(5.1, 7.3, 6.8, 8.0, 6.2) #For testing purposes only. For data requests, please contact the respective database administrator
y <- c(4.9, 5.5, 6.1, 5.7, 5.3) #For testing purposes only. For data requests, please contact the respective database administrator

result <- wilcox.test(x, y, exact = FALSE)
r_value <- rank_biserial(result)
print(r_value)

