

# Load necessary libraries
library(ggplot2)
library(pheatmap)
library(viridis)
# Function to calculate Jaccard similarity between two sets
jaccard_similarity <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union_size <- length(union(set1, set2))
  jaccard_index <- intersection / union_size
  return(jaccard_index)
}
# Function to create a correlation matrix from gene lists
create_correlation_matrix <- function(group1_genes, group2_genes) {
  num_genes_group1 <- length(group1_genes)
  num_genes_group2 <- length(group2_genes)
  # Initialize matrix to store correlation coefficients
  correlation_matrix <- matrix(0, nrow = num_genes_group1, ncol = num_genes_group2)
  # Compute Jaccard similarity for each pair of gene sets
  for (i in 1:num_genes_group1) {
    for (j in 1:num_genes_group2) {
      correlation_matrix[i, j] <- jaccard_similarity(group1_genes[[i]], group2_genes[[j]])
    }
  }
  rownames(correlation_matrix) <- names(group1_genes)
  colnames(correlation_matrix) <- names(group2_genes)
  return(correlation_matrix)
}
# Example gene lists for two groups with multiple subgroups

# Example gene lists for two groups
group1_genes <- read.csv('Tirosh_MPs.csv')
#group1_genes <- read.csv('Neural_sig.csv')
#group1_genes <- read.csv('Neural_sig_Bernstein.csv')
#group1_genes <- read.csv('Neural_sig_Bernstein_overlap.csv')
#group1_genes <- read.csv('~/Documents/Izar_Group/melanoma_bm_cell/Bernstein_signatures.csv')
#group1_genes <- read.csv('~/Documents/Izar_Group/melanoma_bm_cell/Bernstein_signatures.csv')
row.names(group1_genes) <- group1_genes[,1]
group1_genes<-group1_genes[,-1]
rownames(group1_genes)


#group2_genes <- read.csv('~/Documents/Izar_Group/NSCLC/Revision/kinomo/Metaprograms.csv')
#group2_genes <- read.csv('~/Documents/Izar_Group/melanoma_bm_cell/Melanoma_Cell_MPs.csv')
group2_genes <- read.csv('~/Documents/Izar_Group/melanoma_ribas/table_metaprograms_best_top100_complete_7MPs_stouffer.csv')
row.names(group2_genes) <- group2_genes[,1]
group2_genes<-group2_genes[,-1]
rownames(group2_genes)

# group1_genes <- list(
#   subgroup1 = c("gene1", "gene2", "gene3"),
#   subgroup2 = c("gene4", "gene5"),
#   subgroup3 = c("gene6", "gene7", "gene8")
# )
# group2_genes <- list(
#   subgroup1 = c("gene2", "gene3", "gene9"),
#   subgroup2 = c("gene4", "gene10", "gene11", "gene12")
# )
# Calculate correlation matrix
correlation_matrix <- create_correlation_matrix(group1_genes, group2_genes)
#write.csv(correlation_matrix,file='~/Documents/Izar_Group/melanoma_bm_cell/Melanoma_MPs_vs_Bernstein_overlap_correlation_matrix.csv')
#write.csv(correlation_matrix,file='~/Documents/Izar_Group/NSCLC/NSCLC_MPs_vs_Neural_sig_correlation_matrix.csv')
#write.csv(correlation_matrix,file='~/Documents/Izar_Group/NSCLC/NSCLC_MPs_vs_Bernstein_overlap_correlation_matrix.csv')
# Plot correlation matrix
#correlation_matrix <- read.csv('~/Documents/Izar_Group/NSCLC/NSCLC_MPs_vs_Neural_sig_correlation_matrix.csv')
# correlation_matrix <- read.csv('~/Documents/Izar_Group/melanoma_bm_cell/Melanoma_MPs_vs_Bernstein_overlap_correlation_matrix.csv')
# row.names(correlation_matrix) <- correlation_matrix[,1]
# correlation_matrix<-correlation_matrix[,-1]
# #rownames(group1_genes)
library(readr)
#data <- read_csv("Melanoma_MPs_vs_Bernstein_overlap_correlation_matrix.csv")

# Function to convert Jaccard similarity to correlation
jaccard_to_correlation <- function(jaccard_similarity) {
  # Assuming Jaccard similarity is between 0 and 1
  correlation <- (2 * jaccard_similarity) - 1
  return(correlation)
}

# data$correlation_values <- jaccard_to_correlation(data$)
# data$correlation_values
# Example usage:
# A vector of Jaccard similarities
jaccard_similarities <- jaccard_to_correlation(correlation_matrix)

# Convert to correlation values
correlation_values <- jaccard_to_correlation(jaccard_similarities)

# Print the correlation values
print(correlation_values)

# correlation_matrix <- read.csv('~/Documents/Izar_Group/melanoma_bm_cell/Melanoma_MPs_vs_Bernstein_overlap_correlation_matrix_correlation.csv')
# row.names(correlation_matrix) <- correlation_matrix[,1]
# correlation_matrix<-correlation_matrix[,-1]

plot_correlation_matrix <- function(correlation_matrix) {
  # Plot correlation matrix using pheatmap
  pheatmap(as.matrix(correlation_matrix),
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           color = colorRampPalette(c("white", "gold"))(10),
           #color = scale_color_viridis(discrete = TRUE),
           #color = scale_color_viridis(), 
          # border_color = NA,
         #  main = "Correlation: NSCLC MPs\n Neural-like sig Bernstein")
           main = "Similarity (Jaccard): Melanoma Ribas MPs\n Gavish et al 2023 MPs")
           # main = "Correlation: NSCLC \n Gavish et al 2023 MPs")
}
# Plot the correlation matrix
#pdf('plot_correlation_matrix_kinomo_MPs_Tirosh_MPs.pdf',height=10,width=5)
pdf('plot_correlation_matrix_ribas_MPs_Tirosh_MPs.pdf',height=10,width=5)
#pdf('plot_correlation_matrix_kinomo_MPs_Tirosh_MPs_1.pdf',height=10,wi2dth=5)
#pdf('plot_correlation_matrix_kinomo_MPs_Tirosh_MPs_2.pdf',height=10,width=5)
#pdf('~/Documents/Izar_Group/NSCLC/NSCLC_MPs_vs_Neural_sig_Bernstein.pdf',height=3,width=7)
#pdf('~/Documents/Izar_Group/NSCLC/NSCLC_MPs_vs_Bernstein_overlap.pdf',height=3,width=7)
#pdf('~/Documents/Izar_Group/melanoma_bm_cell/Melanoma_MPs_vs_Bernstein_overlap.pdf',height=3,width=7)
plot_correlation_matrix(correlation_matrix)
dev.off()



#### pvalues

# Load necessary libraries
library(ggplot2)
library(pheatmap)
library(viridis)

# Function to calculate Jaccard similarity between two sets
jaccard_similarity <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union_size <- length(union(set1, set2))
  jaccard_index <- intersection / union_size
  return(jaccard_index)
}

# Function to create a correlation matrix from gene lists and compute p-values
create_correlation_matrix_with_pvalues <- function(group1_genes, group2_genes, num_permutations = 1000) {
  num_genes_group1 <- length(group1_genes)
  num_genes_group2 <- length(group2_genes)
  
  # Initialize matrix to store correlation coefficients
  correlation_matrix <- matrix(0, nrow = num_genes_group1, ncol = num_genes_group2)
  pvalue_matrix <- matrix(0, nrow = num_genes_group1, ncol = num_genes_group2)
  
  # Compute Jaccard similarity and p-values for each pair of gene sets
  for (i in 1:num_genes_group1) {
    for (j in 1:num_genes_group2) {
      # Original Jaccard similarity
      original_similarity <- jaccard_similarity(group1_genes[[i]], group2_genes[[j]])
      correlation_matrix[i, j] <- original_similarity
      
      # Permutation test
      permuted_similarities <- numeric(num_permutations)
      combined_genes <- c(group1_genes[[i]], group2_genes[[j]])
      for (perm in 1:num_permutations) {
        # Permute the combined genes and split into two sets
        permuted_genes <- sample(combined_genes)
        permuted_set1 <- permuted_genes[1:length(group1_genes[[i]])]
        permuted_set2 <- permuted_genes[(length(group1_genes[[i]]) + 1):length(permuted_genes)]
        permuted_similarities[perm] <- jaccard_similarity(permuted_set1, permuted_set2)
      }
      
      # Calculate p-value: proportion of permuted similarities >= original similarity
      pvalue_matrix[i, j] <- mean(permuted_similarities >= original_similarity)
    }
  }
  
  rownames(correlation_matrix) <- names(group1_genes)
  colnames(correlation_matrix) <- names(group2_genes)
  rownames(pvalue_matrix) <- names(group1_genes)
  colnames(pvalue_matrix) <- names(group2_genes)
  
  list(correlation_matrix = correlation_matrix, pvalue_matrix = pvalue_matrix)
}

# Example of running the function
# Assume group1_genes and group2_genes are your lists of gene sets for two groups
result <- create_correlation_matrix_with_pvalues(group1_genes, group2_genes, num_permutations = 1000)

# To access the matrices

correlation_matrix <- result$correlation_matrix
pvalue_matrix <- result$pvalue_matrix

### For hard-coded correlation matrix as csv file

# Load necessary libraries
library(ggplot2)
library(pheatmap)
library(viridis)

# Function to calculate Jaccard similarity between two sets
jaccard_similarity <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union_size <- length(union(set1, set2))
  jaccard_index <- intersection / union_size
  return(jaccard_index)
}

# Function to calculate p-values for provided correlation matrix
calculate_pvalues_from_matrix <- function(correlation_matrix, group1_genes, group2_genes, num_permutations = 1000) {
  num_genes_group1 <- length(group1_genes)
  num_genes_group2 <- length(group2_genes)
  
  # Initialize matrix to store p-values
  pvalue_matrix <- matrix(0, nrow = num_genes_group1, ncol = num_genes_group2)
  
  # Perform permutation test to compute p-values
  for (i in 1:num_genes_group1) {
    for (j in 1:num_genes_group2) {
      original_similarity <- correlation_matrix[i, j]
      
      # Permutation test
      permuted_similarities <- numeric(num_permutations)
      combined_genes <- c(group1_genes[[i]], group2_genes[[j]])
      for (perm in 1:num_permutations) {
        # Permute the combined genes and split into two sets
        permuted_genes <- sample(combined_genes)
        permuted_set1 <- permuted_genes[1:length(group1_genes[[i]])]
        permuted_set2 <- permuted_genes[(length(group1_genes[[i]]) + 1):length(permuted_genes)]
        permuted_similarities[perm] <- jaccard_similarity(permuted_set1, permuted_set2)
      }
      
      # Calculate p-value: proportion of permuted similarities >= original similarity
      pvalue_matrix[i, j] <- mean(permuted_similarities >= original_similarity)
    }
  }
  
  rownames(pvalue_matrix) <- rownames(correlation_matrix)
  colnames(pvalue_matrix) <- colnames(correlation_matrix)
  
  return(pvalue_matrix)
}

# Example usage:
# Load the correlation matrix from CSV
correlation_matrix <- as.matrix(read.csv("~/Documents/Izar_Group/melanoma_bm_cell/Melanoma_MPs_vs_Bernstein_overlap_correlation_matrix_correlation.csv", row.names = 1))

# Assume group1_genes and group2_genes are your lists of gene sets
pvalue_matrix <- calculate_pvalues_from_matrix(correlation_matrix, group1_genes, group2_genes, num_permutations = 1000)

# Save the p-value matrix as a CSV file
# write.csv(pvalue_matrix, "pvalue_matrix.csv")

