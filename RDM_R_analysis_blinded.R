# Correlation Matrix Comparison Analysis
# Last updated: 8/14/2025

#libraries & helpful functions
read.csv3 <- function(x){read.csv(x, header = F)}
library(proxy)
library(tibble)
library(magrittr)
library(resemble)


#Pull in data into a set of lists
asd_list <- list.files(path = "...\\individual_RDMs\\ASD", pattern = "\\.csv$", full.names = TRUE)
asd_original <- lapply(asd_list, read.csv3)
names(asd_original) <- tools::file_path_sans_ext(basename(asd_list))

nt_list <- list.files( path = "...\\individual_RDMs\\NT", pattern = "\\.csv$", full.names = TRUE)
nt_original <- lapply(nt_list, read.csv3)
names(nt_original) <- tools::file_path_sans_ext(basename(nt_list))

#note about order: house1 house2 house3 person1 person2 person3 invperson1 invperson2 invperson3
#note about order: person 1 = white, person 2 = black, person 3 = asian

asd_original <- lapply(asd_original, as.matrix)
nt_original <- lapply(nt_original, as.matrix)

#Create dissimilarity matrices based on original correlation matrices, where dissimilarity = 1 - cor
asd <- lapply(asd_original, function(mat) 1 - mat)
nt <- lapply(nt_original, function(mat) 1 - mat)
asd <- lapply(asd, as.matrix)
nt <- lapply(nt, as.matrix)


########Frobenius (Euclidian) distance approach
# Frobenius norm of matrix difference
frobenius_distance <- function(mat1, mat2) {
        norm(mat1 - mat2, type = "F") #equivalent to sqrt(sum((mat1 - mat2)^2))
}

# Compute mean matrix for a list of matrices
mean_matrix <- function(matrix_list) {
        Reduce("+", matrix_list) / length(matrix_list)
}

# Elementwise mean matrix (returns array)
stack_matrices <- function(matrix_list) {
        abind::abind(matrix_list, along = 3)
}
elementwise_between_perm_test <- function(asd, nt, n_perm = 10000, alpha = 0.05, correction = "fdr") {
        stopifnot(length(asd) > 0, length(nt) > 0)
        
        rdm_dim <- dim(asd[[1]])
        stopifnot(!is.null(rdm_dim), length(rdm_dim) == 2)
        
        all_data <- c(asd, nt)
        n_asd <- length(asd)
        n_total <- length(all_data)
        
        mean_asd <- Reduce("+", asd) / n_asd
        mean_nt  <- Reduce("+", nt)  / length(nt)
        observed_diff <- mean_asd - mean_nt
        
        # Store p-values in a flat vector
        p_vec <- numeric(length = prod(rdm_dim))
        
        for (k in seq_along(p_vec)) {
                i <- ((k - 1) %% rdm_dim[1]) + 1
                j <- ((k - 1) %/% rdm_dim[1]) + 1
                
                asd_vals <- sapply(asd, function(r) r[i, j])
                nt_vals  <- sapply(nt,  function(r) r[i, j])
                obs_diff <- mean(asd_vals) - mean(nt_vals)
                
                all_vals <- c(asd_vals, nt_vals)
                
                perm_diffs <- replicate(n_perm, {
                        idx_asd <- sample(n_total, n_asd)
                        idx_nt  <- setdiff(1:n_total, idx_asd)
                        mean(all_vals[idx_asd]) - mean(all_vals[idx_nt])
                })
                
                p_vec[k] <- mean(abs(perm_diffs) >= abs(obs_diff))
        }
        
        # Reshape into matrix
        p_matrix <- matrix(p_vec, nrow = rdm_dim[1], ncol = rdm_dim[2])
        
        # Multiple comparisons corrections
        p_fdr_matrix  <- matrix(p.adjust(p_vec, method = "fdr"), nrow = rdm_dim[1])
        p_bonf_matrix <- matrix(p.adjust(p_vec, method = "bonferroni"), nrow = rdm_dim[1])
        
        # Determine which to use for significance thresholding
        if (correction == "fdr") {
                significant_mask <- p_fdr_matrix < alpha
                p_used <- p_fdr_matrix
        } else if (correction == "bonferroni") {
                significant_mask <- p_bonf_matrix < alpha
                p_used <- p_bonf_matrix
        } else {
                significant_mask <- p_matrix < alpha
                p_used <- p_matrix
        }
        
        thresholded_diff <- observed_diff
        thresholded_diff[!significant_mask] <- 0
        
        return(list(
                observed_diff = observed_diff,
                p_matrix_raw = p_matrix,
                p_matrix_fdr = p_fdr_matrix,
                p_matrix_bonf = p_bonf_matrix,
                significant_mask = significant_mask,
                thresholded_diff = thresholded_diff,
                correction_used = correction,
                alpha = alpha,
                p_matrix_used = p_used
        ))
}


# Inputs: asd and nt (each with k correlation matrices)

# Compute group mean matrices
mean_asd <- mean_matrix(asd) %>% as.matrix()
mean_nt  <- mean_matrix(nt) %>% as.matrix()
diff_mat <- mean_asd - mean_nt %>% as.matrix()
obs_stat <- frobenius_distance(mean_asd, mean_nt)

# Set seed
set.seed(43058)
# Permutation test
n_perm <- 10000
perm_stats2 <- numeric(n_perm)
all_data <- c(asd, nt)

for (i in 1:n_perm) {
        perm_idx <- sample(1:(length(asd)+length(nt)), length(asd)) #randomly select asd length IDs
        perm_asd <- all_data[perm_idx]
        perm_nt <- all_data[-perm_idx]
        
        mean_perm_asd <- mean_matrix(perm_asd)
        mean_perm_nt  <- mean_matrix(perm_nt)
        
        perm_stats2[i] <- frobenius_distance(mean_perm_asd, mean_perm_nt)
}

# Compute empirical p-value
p_value <- mean(perm_stats2 >= obs_stat)
p_value

obs_stat
mean(perm_stats2) #average difference of permutation test
obs_stat - mean(perm_stats2)
sd(perm_stats2)

# Output results
cat("Observed Frobenius distance (Euclidian norm of a matrix) between group means:", round(obs_stat, 3), "\n")
cat("Permutation-based p-value:", round(p_value, 10), "\n")

# Plot permutation null distribution
par(mar=c(5,6,4,1.3), mgp = c(3.5, 1, 0))
hist(perm_stats2, breaks = 50, main = expression(paste("Distribution of Permutation Frobenius Distances, ", Delta[italic(F)[italic(j)]])),
     cex.main = 1.2,
     cex.lab = 1.2,
     cex.axis = 1.2,
     xaxs = "i", yaxs = "i", 
     xlab = expression(Delta[italic(F)[italic(j)]]), 
     col = "lightgrey",
     xlim = c(.2,.7),
     ylim = c(0, 900))
abline(v = obs_stat, col = "black", lwd = 2, lty = 2)
text(x = .6, y = 600, labels = bquote(Delta[italic(F)[0]] == .(round(obs_stat, 3))), col = "black", cex = 1.2)
text(x = .6, y = 530, labels = bquote(italic(p)[italic(F)[0]] == .(format(round(p_value, 3), nsmall = 3))), col = "black", cex = 1.2)


library(corrplot)
par(mfrow=c(1,3))
corrplot(mean_asd, 
         method = "color", # Other options: "square", "ellipse", "number", "shade", "color", "pie"
         type = "lower",    # Other options: "lower", "full"
         tl.col = "black",  # Color of text labels
         tl.srt = 0,       #  Rotation of text labels
         title = "Average ASD RDM", mar = c(0, 2, 0, 0), line = -5, 
         tl.pos = "l",
         cl.pos = "n",
         addCoef.col = "grey", # Add correlation coefficients as text
         diag = T       # Remove diagonal elements
)

corrplot(mean_nt, 
         method = "color", # Other options: "square", "ellipse", "number", "shade", "color", "pie"
         type = "lower",    # Other options: "lower", "full"
         tl.col = "black",  # Color of text labels
         tl.srt = 0,       # Rotation of text labels
         title = "Average NT RDM", mar = c(0, 2, 0, 0), line = -5,
         tl.pos = "n",
         cl.pos = "n",
         addCoef.col = "grey", # Add correlation coefficients as text
         diag = T       # Remove diagonal elements
)

corrplot(diff_mat, 
         method = "color", # Other options: "square", "ellipse", "number", "shade", "color", "pie"
         type = "lower",    # Other options: "lower", "full"
         tl.col = "black",  # Color of text labels
         tl.srt = 0,       # Rotation of text labels
         title = "ASD - NT RDM", mar = c(0, 2, 0, 0), line = -5,
         tl.pos = "n",
         cl.pos = "n",
         addCoef.col = "grey30", # Add correlation coefficients as text
         diag = T       # Remove diagonal elements
)

par(mfrow=c(1,1)) 


#Result: Based on the permutation test that adjusted for squared differences over all cells to test whether the groups differ in overall RDM geometry.
#The observed Frobenius distance between group mean matrices is 0.643 and a permutation test showed that this effect was much larger than the randomly computed 
#difference (.380) based on 10,000 permutations, holding group sizes to their original sample size values (p < .001). 
#Thus, while the groups overall look different, there is no overall difference effect shown.
#Issue with this result: it does not norm the differences nor adjust appropriately for Euclidian distance. 

########## Look at results for the condition_1 (first 3 elements = houses)
# Compute group mean matrices
asd_1 <- lapply(asd, function(mat) {mat[1:3, 1:3]})
nt_1 <- lapply(nt, function(mat) {mat[1:3, 1:3]})
mean_asd_1 <- mean_matrix(asd_1) %>% as.matrix()
mean_nt_1  <- mean_matrix(nt_1) %>% as.matrix()
diff_mat_1 <- mean_asd_1 - mean_nt_1 %>% as.matrix()
obs_stat_1 <- frobenius_distance(mean_asd_1, mean_nt_1)

# Set seed for reproducibility
set.seed(43058)
# Permutation test
n_perm <- 10000
perm_stats2_1 <- numeric(n_perm)
all_data_1 <- c(asd_1, nt_1)

for (i in 1:n_perm) {
        perm_idx <- sample(1:(length(asd_1)+length(nt_1)), length(asd_1)) #randomly select asd length IDs
        perm_asd <- all_data_1[perm_idx]
        perm_nt <- all_data_1[-perm_idx]
        
        mean_perm_asd <- mean_matrix(perm_asd)
        mean_perm_nt  <- mean_matrix(perm_nt)
        
        perm_stats2_1[i] <- frobenius_distance(mean_perm_asd, mean_perm_nt)
}

# Compute empirical p-value
p_value_1 <- mean(perm_stats2_1 >= obs_stat_1)
p_value_1

obs_stat_1
mean(perm_stats2_1) #average difference of permutation test
obs_stat_1 - mean(perm_stats2_1)
sd(perm_stats2_1)

# Output results
cat("Observed Frobenius distance (Euclidian norm of a matrix) between house group means:", round(obs_stat_1, 3), "\n")
cat("Permutation-based p-value:", round(p_value_1, 10), "\n")

########## Look at results for the condition_2 (middle 3 elements = faces, upright)
# Compute group mean matrices
asd_2 <- lapply(asd, function(mat) {mat[4:6, 4:6]})
nt_2 <- lapply(nt, function(mat) {mat[4:6, 4:6]})
mean_asd_2 <- mean_matrix(asd_2) %>% as.matrix()
mean_nt_2  <- mean_matrix(nt_2) %>% as.matrix()
diff_mat_2 <- mean_asd_2 - mean_nt_2 %>% as.matrix()
obs_stat_2 <- frobenius_distance(mean_asd_2, mean_nt_2)

# Set seed for reproducibility
set.seed(43058)
# Permutation test
n_perm <- 10000
perm_stats2_2 <- numeric(n_perm)
all_data_2 <- c(asd_2, nt_2)

for (i in 1:n_perm) {
        perm_idx <- sample(1:(length(asd_2)+length(nt_2)), length(asd_2)) #randomly select asd length IDs
        perm_asd <- all_data_2[perm_idx]
        perm_nt <- all_data_2[-perm_idx]
        
        mean_perm_asd <- mean_matrix(perm_asd)
        mean_perm_nt  <- mean_matrix(perm_nt)
        
        perm_stats2_2[i] <- frobenius_distance(mean_perm_asd, mean_perm_nt)
}

# Compute empirical p-value
p_value_2 <- mean(perm_stats2_2 >= obs_stat_2)
pvalue2 <- (sum(ifelse(perm_stats2_2 >= obs_stat_2, 1, 0))+1) / (length(perm_stats2_2) + 1)
p_value_2
pvalue2

obs_stat_2
mean(perm_stats2_2) #average difference of permutation test
obs_stat_2 - mean(perm_stats2_2)
sd(perm_stats2_2)

# Output results
cat("Observed Frobenius distance (Euclidian norm of a matrix) between faces group means:", round(obs_stat_2, 3), "\n")
cat("Permutation-based p-value:", round(p_value_2, 10), "\n")



########## Look at results for the condition_3 (last 3 elements: inverse faces)
# Compute group mean matrices
asd_3 <- lapply(asd, function(mat) {mat[7:9, 7:9]})
nt_3 <- lapply(nt, function(mat) {mat[7:9, 7:9]})
mean_asd_3 <- mean_matrix(asd_3) %>% as.matrix()
mean_nt_3  <- mean_matrix(nt_3) %>% as.matrix()
diff_mat_3 <- mean_asd_3 - mean_nt_3 %>% as.matrix()
obs_stat_3 <- frobenius_distance(mean_asd_3, mean_nt_3)

# Set seed for reproducibility
set.seed(43058)
# Permutation test
n_perm <- 10000
perm_stats2_3 <- numeric(n_perm)
all_data_3 <- c(asd_3, nt_3)

for (i in 1:n_perm) {
        perm_idx <- sample(1:(length(asd_3)+length(nt_3)), length(asd_3)) #randomly select asd length IDs
        perm_asd <- all_data_3[perm_idx]
        perm_nt <- all_data_3[-perm_idx]
        
        mean_perm_asd <- mean_matrix(perm_asd)
        mean_perm_nt  <- mean_matrix(perm_nt)
        
        perm_stats2_3[i] <- frobenius_distance(mean_perm_asd, mean_perm_nt)
}

# Compute empirical p-value
p_value_3 <- mean(perm_stats2_3 >= obs_stat_3)
p_value_3

obs_stat_3
mean(perm_stats2_3) #average difference of permutation test
obs_stat_3 - mean(perm_stats2_3)
sd(perm_stats2_3)

# Output results
cat("Observed Frobenius distance (Euclidian norm of a matrix) between inverse faces group means:", round(obs_stat_3, 3), "\n")
cat("Permutation-based p-value:", round(p_value_3, 10), "\n")

cat("Observed Frobenius Distance Between ASD and NT for Houses:", round(obs_stat_1, 3), "\n",
    "Permutation p = ", round(p_value_1, 3), "\n",
    "Observed Frobenius Distance Between ASD and NT for Faces:", round(obs_stat_2, 3), "\n",
    "Permutation p = ", round(p_value_2, 3), "\n",
    "Observed Frobenius Distance Between ASD and NT Inverse Faces:", round(obs_stat_3, 3), "\n",
    "Permutation p = ", round(p_value_3, 3), "\n")

delta1 = round(obs_stat_1, 3)
delta2 = round(obs_stat_2, 3)
delta3 = round(obs_stat_3, 3)
pf1 = round(p_value_1, 3)
pf2 = round(p_value_2, 3)
pf3 = round(p_value_3, 3)

par(mfrow=c(1,3), mar=c(5,6,4,1.3), mgp = c(4, 1, 0))
# Plot permutation null distribution
hist(perm_stats2_1, breaks = 50, main = "House Conditions", 
     cex.main = 1.5,
     cex.lab = 1.5,
     cex.axis = 1.5,
     xaxs = "i", yaxs = "i", 
     xlab = "", col = "grey90",
     xlim = c(0,.35),
     ylim = c(0, 550))
abline(v = obs_stat_1, col = "black", lwd = 2, lty = 2)
text(x = .246, y = 400, labels = bquote(Delta[italic(F)[0]] == .(delta1)), col = "black", cex = 1.5)
text(x = .246, y = 350, labels = bquote(italic(p)[italic(F)[0]] == .(pf1)), col = "black", cex = 1.5)

# Plot permutation null distribution
hist(perm_stats2_2, breaks = 50, main = paste("Upright Face Conditions"), 
     cex.main = 1.5,
     cex.lab = 1.5,
     cex.axis = 1.5,
     xaxs = "i", yaxs = "i", 
     xlab = expression(paste("Null Distribution of Frobenius Distances, ", Delta[italic(F)[italic(j)]])), col = "grey90",
     xlim = c(0,.35),
     ylim = c(0, 550),
     ylab = "")
abline(v = obs_stat_2, col = "black", lwd = 2, lty = 2)
text(x = .252, y = 400, labels = bquote(Delta[italic(F)[0]] == .(delta2)), col = "black", cex = 1.5)
text(x = .252, y = 350, labels = bquote(italic(p)[italic(F)[0]] == .(pf2)), col = "black", cex = 1.5)

# Plot permutation null distribution
hist(perm_stats2_3, breaks = 50, main = "Inverted Face Conditions", 
     cex.main = 1.5,
     cex.lab = 1.5,
     cex.axis = 1.5,
     xaxs = "i", yaxs = "i", 
     xlab = "", col = "grey90",
     xlim = c(0,.35),
     ylim = c(0, 550),
     ylab = "")
abline(v = obs_stat_3, col = "black", lwd = 2, lty = 2)
text(x = .248, y = 400, labels = bquote(Delta[italic(F)[0]] == .(delta3)), col = "black", cex = 1.5)
text(x = .248, y = 350, labels = bquote(italic(p)[italic(F)[0]] == .(pf3)), col = "black", cex = 1.5)
par(mfrow=c(1,1), mgp = c(3, 1, 0))