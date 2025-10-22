#!/usr/bin/env Rscript

# PTB Dataset Analysis with WICMAD
# This script loads the PTB dataset, runs WICMAD clustering, and visualizes results

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(cluster)
# library(caret)
library(mclust)  # For adjustedRandIndex

# Load WICMAD package
library(devtools)
devtools::install()

# Function to load PTB dataset with imbalanced sampling
load_ptb_data <- function(data_dir, reveal_ratio = 0.15) {
  cat("Loading PTB dataset...\n")
  
  # Load database and SCP statements
  database_file <- file.path(data_dir, "ptbxl_database.csv")
  scp_file <- file.path(data_dir, "scp_statements.csv")
  
  if (!file.exists(database_file)) {
    stop("Database file not found: ", database_file)
  }
  
  # Read database
  database <- read.csv(database_file)
  scp_statements <- read.csv(scp_file)
  
  cat("Total records in database:", nrow(database), "\n")
  
  # For simplicity, let's focus on normal vs abnormal classification
  # We'll use the scp_codes to determine normal vs abnormal
  normal_codes <- c("NORM")
  
  # Create binary labels: 0 = normal, 1 = abnormal
  labels <- integer(nrow(database))
  
  for (i in seq_len(nrow(database))) {
    scp_codes_str <- database$scp_codes[i]
    if (is.na(scp_codes_str) || scp_codes_str == "") {
      labels[i] <- 1  # Unknown = abnormal
    } else {
      # Parse the scp_codes string to check for normal
      if (grepl("NORM", scp_codes_str)) {
        labels[i] <- 0  # Normal
      } else {
        labels[i] <- 1  # Abnormal
      }
    }
  }
  
  # Create imbalanced dataset: 1000 total samples with 5% anomalies
  cat("Creating imbalanced dataset...\n")
  cat("Target: 1000 total samples with 5% anomalies (50 abnormal, 950 normal)\n")
  
  # Find normal and abnormal indices
  normal_indices <- which(labels == 0)
  abnormal_indices <- which(labels == 1)
  
  cat("Normal samples available:", length(normal_indices), "\n")
  cat("Abnormal samples available:", length(abnormal_indices), "\n")
  
  # Target: 1000 total samples, 5% anomalies = 50 abnormal, 950 normal
  target_total <- 1000
  target_abnormal <- 50
  target_normal <- 950
  
  # Sample abnormal samples
  if (length(abnormal_indices) >= target_abnormal) {
    selected_abnormal <- sample(abnormal_indices, target_abnormal)
  } else {
    # If not enough abnormal samples, sample with replacement
    selected_abnormal <- sample(abnormal_indices, target_abnormal, replace = TRUE)
  }
  
  # Sample normal samples
  if (length(normal_indices) >= target_normal) {
    selected_normal <- sample(normal_indices, target_normal)
  } else {
    # If not enough normal samples, sample with replacement
    selected_normal <- sample(normal_indices, target_normal, replace = TRUE)
  }
  
  # Combine selected samples
  selected_indices <- c(selected_normal, selected_abnormal)
  
  # Shuffle the indices to randomize order
  selected_indices <- sample(selected_indices)
  
  imbalanced_labels <- labels[selected_indices]
  
  cat("Final imbalanced dataset:\n")
  cat("Total samples:", length(selected_indices), "\n")
  cat("Normal:", sum(imbalanced_labels == 0), "\n")
  cat("Abnormal:", sum(imbalanced_labels == 1), "\n")
  cat("Abnormal percentage:", round(mean(imbalanced_labels == 1) * 100, 1), "%\n")
  
  # Load real PTB ECG data from .dat files
  cat("Loading real PTB ECG data from .dat files...\n")
  
  n_samples <- length(selected_indices)
  series_list <- list()
  
  # Get the filenames for the selected samples
  selected_filenames <- database$filename_lr[selected_indices]
  
  cat("Loading", n_samples, "real ECG samples...\n")
  
  for (i in seq_len(n_samples)) {
    if (i %% 100 == 0) {
      cat("Loading sample", i, "of", n_samples, "\n")
    }
    
    # Get the filename for this sample
    filename <- selected_filenames[i]
    
    # Construct full path to the .dat file (add .dat extension)
    dat_file <- file.path(data_dir, paste0(filename, ".dat"))
    
    if (file.exists(dat_file)) {
      # Read the corresponding .hea file to get metadata
      hea_file <- gsub("\\.dat$", ".hea", dat_file)
      
      if (file.exists(hea_file)) {
        # Read header file to get number of leads and samples
        hea_lines <- readLines(hea_file)
        first_line <- hea_lines[1]
        parts <- strsplit(first_line, " ")[[1]]
        
        if (length(parts) >= 4) {
          n_leads <- as.numeric(parts[2])
          n_samples_ecg <- as.numeric(parts[3])
          
          # Read the binary .dat file
          # PTB uses 16-bit signed integers
          con <- file(dat_file, "rb")
          ecg_data <- readBin(con, integer(), n = n_leads * n_samples_ecg, size = 2, signed = TRUE, endian = "little")
          close(con)
          
          # Reshape data into matrix (samples x leads)
          ecg_matrix <- matrix(ecg_data, nrow = n_samples_ecg, ncol = n_leads, byrow = FALSE)
          
          # Convert to mV (assuming the data is in ADC units, convert to mV)
          # This is a rough conversion - in practice you'd use the gain values from the header
          ecg_matrix <- ecg_matrix / 1000.0  # Rough conversion to mV
          
          series_list[[i]] <- ecg_matrix
        } else {
          # Fallback: create a simple matrix if header parsing fails
          cat("Warning: Could not parse header for sample", i, ", using fallback\n")
          ecg_matrix <- matrix(rnorm(100 * 12), nrow = 100, ncol = 12)
          series_list[[i]] <- ecg_matrix
        }
      } else {
        # Fallback: create a simple matrix if .hea file doesn't exist
        cat("Warning: No .hea file for sample", i, ", using fallback\n")
        ecg_matrix <- matrix(rnorm(100 * 12), nrow = 100, ncol = 12)
        series_list[[i]] <- ecg_matrix
      }
    } else {
      # Debug: print the actual path being tried
      cat("Debug: Trying to load file:", dat_file, "\n")
      cat("Debug: File exists:", file.exists(dat_file), "\n")
      
      # Try to find the file in the records100 directory
      # Extract the record ID from the filename
      record_id <- gsub(".*/([0-9]+)_lr$", "\\1", filename)
      if (nchar(record_id) >= 5) {
        # Get the first 5 digits for the directory
        dir_prefix <- substr(record_id, 1, 5)
        # Pad with zeros to get the directory name
        dir_name <- sprintf("%05d", as.numeric(dir_prefix))
        
        # Try the records100 directory
        alt_dat_file <- file.path(data_dir, "records100", dir_name, paste0(record_id, "_lr.dat"))
        alt_hea_file <- file.path(data_dir, "records100", dir_name, paste0(record_id, "_lr.hea"))
        
        cat("Debug: Trying alternative path:", alt_dat_file, "\n")
        cat("Debug: Alternative file exists:", file.exists(alt_dat_file), "\n")
        
        if (file.exists(alt_dat_file) && file.exists(alt_hea_file)) {
          # Read the header file
          hea_lines <- readLines(alt_hea_file)
          first_line <- hea_lines[1]
          parts <- strsplit(first_line, " ")[[1]]
          
          if (length(parts) >= 4) {
            n_leads <- as.numeric(parts[2])
            n_samples_ecg <- as.numeric(parts[3])
            
            # Read the binary .dat file
            con <- file(alt_dat_file, "rb")
            ecg_data <- readBin(con, integer(), n = n_leads * n_samples_ecg, size = 2, signed = TRUE, endian = "little")
            close(con)
            
            # Reshape data into matrix (samples x leads)
            ecg_matrix <- matrix(ecg_data, nrow = n_samples_ecg, ncol = n_leads, byrow = FALSE)
            
            # Convert to mV
            ecg_matrix <- ecg_matrix / 1000.0
            
            series_list[[i]] <- ecg_matrix
            cat("Successfully loaded real ECG data for sample", i, "\n")
          } else {
            cat("Warning: Could not parse header for sample", i, ", using fallback\n")
            ecg_matrix <- matrix(rnorm(100 * 12), nrow = 100, ncol = 12)
            series_list[[i]] <- ecg_matrix
          }
        } else {
          # Fallback: create a simple matrix if file still doesn't exist
          cat("Warning: No .dat file found for sample", i, ", using fallback\n")
          ecg_matrix <- matrix(rnorm(100 * 12), nrow = 100, ncol = 12)
          series_list[[i]] <- ecg_matrix
        }
      } else {
        # Fallback: create a simple matrix if record ID parsing fails
        cat("Warning: Could not parse record ID for sample", i, ", using fallback\n")
        ecg_matrix <- matrix(rnorm(100 * 12), nrow = 100, ncol = 12)
        series_list[[i]] <- ecg_matrix
      }
    }
  }
  
  return(list(
    train_series = series_list,
    train_labels = imbalanced_labels,
    test_series = NULL,
    test_labels = NULL
  ))
}

# Function to interpolate time series to target dimension
interpolate_series <- function(series, target_dim = 16) {
  n <- nrow(series)
  if (n == target_dim) {
    return(series)
  }
  
  # Interpolate each dimension separately
  interpolated_series <- matrix(0, nrow = target_dim, ncol = ncol(series))
  
  for (j in seq_len(ncol(series))) {
    interpolated_series[, j] <- approx(seq(0, 1, length.out = n), series[, j], 
                                      seq(0, 1, length.out = target_dim))$y
  }
  
  return(interpolated_series)
}

# Function to prepare data for WICMAD
prepare_wicmad_data <- function(series_list, labels) {
  cat("Preparing data for WICMAD...\n")
  
  n_series <- length(series_list)
  
  # Interpolate each series to 16 dimensions (for testing)
  interpolated_series <- list()
  for (i in seq_len(n_series)) {
    interpolated_series[[i]] <- interpolate_series(series_list[[i]], 16)
  }
  
  # Create time coordinates
  t <- seq(0, 1, length.out = 16)
  
  # Convert to list format for WICMAD (each element is a PxM matrix)
  # P=16 (time points), M=12 (ECG leads)
  Y <- interpolated_series
  
  return(list(Y = Y, t = t, labels = labels))
}

# Function to calculate clustering metrics
calculate_clustering_metrics <- function(true_labels, pred_labels) {
  true <- factor(true_labels, levels = sort(unique(c(true_labels, pred_labels))))
  pred <- factor(pred_labels, levels = sort(unique(c(true_labels, pred_labels))))
  
  cm <- table(true, pred)
  cat("Confusion Matrix:\n")
  print(cm)
  
  classes <- levels(true)
  metrics <- list()
  
  for (class in classes) {
    tp <- sum(true == class & pred == class)
    fp <- sum(true != class & pred == class)
    fn <- sum(true == class & pred != class)
    tn <- sum(true != class & pred != class)
    
    precision <- if (tp + fp == 0) 0 else tp / (tp + fp)
    recall <- if (tp + fn == 0) 0 else tp / (tp + fn)
    specificity <- if (tn + fp == 0) 0 else tn / (tn + fp)
    f1 <- if (precision + recall == 0) 0 else 2 * precision * recall / (precision + recall)
    accuracy <- (tp + tn) / (tp + tn + fp + fn)
    
    metrics[[paste0("Class_", class)]] <- list(
      Precision = precision,
      Recall = recall,
      Specificity = specificity,
      F1 = f1,
      Accuracy = accuracy
    )
  }
  
  macro_precision <- mean(sapply(metrics, function(x) x$Precision), na.rm = TRUE)
  macro_recall <- mean(sapply(metrics, function(x) x$Recall), na.rm = TRUE)
  macro_specificity <- mean(sapply(metrics, function(x) x$Specificity), na.rm = TRUE)
  macro_f1 <- mean(sapply(metrics, function(x) x$F1), na.rm = TRUE)
  overall_accuracy <- mean(sapply(metrics, function(x) x$Accuracy), na.rm = TRUE)
  
  ari <- mclust::adjustedRandIndex(true_labels, pred_labels)
  
  return(list(
    Confusion_Matrix = cm,
    Macro_Precision = macro_precision,
    Macro_Recall = macro_recall,
    Macro_Specificity = macro_specificity,
    Macro_F1 = macro_f1,
    Overall_Accuracy = overall_accuracy,
    ARI = ari,
    Class_Metrics = metrics
  ))
}

# Function to plot overlapped data
plot_overlapped_data <- function(series_list, labels, title = "PTB Dataset - Original ECG Data", max_series = NULL) {
  cat("Creating overlapped data plot...\n")
  
  # Use all series if max_series is NULL
  if (is.null(max_series)) {
    series_indices <- seq_len(length(series_list))
  } else {
    n_series <- min(max_series, length(series_list))
    series_indices <- sample(seq_len(length(series_list)), n_series)
  }
  
  # Create plotting data
  plot_data <- data.frame()
  
  for (i in series_indices) {
    series_data <- series_list[[i]]
    time_points <- seq(0, 1, length.out = nrow(series_data))
    
    for (j in seq_len(ncol(series_data))) {
      series_df <- data.frame(
        time = time_points,
        value = series_data[, j],
        dimension = paste("Lead", j),
        series_id = i,
        label = labels[i],
        label_name = paste("Class", labels[i]),
        row.names = NULL
      )
      
      plot_data <- rbind(plot_data, series_df)
    }
  }
  
  p <- ggplot(plot_data, aes(x = .data$time, y = .data$value, color = .data$label_name, group = .data$series_id)) +
    geom_line(alpha = 0.6, linewidth = 0.5) +
    facet_wrap(~ .data$dimension, scales = "free_y", ncol = 3) +
    labs(title = title,
         x = "Time", y = "Amplitude", color = "Class") +
    scale_color_manual(values = c("Class 0" = "blue", "Class 1" = "red")) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

# Function to plot clustering results
plot_clustering_results <- function(series_list, true_labels, cluster_assignments, title = "PTB - Clustering Results", max_series = NULL) {
  cat("Creating clustering results plot...\n")
  
  # Use all series if max_series is NULL
  if (is.null(max_series)) {
    series_indices <- seq_len(length(series_list))
  } else {
    n_series <- min(max_series, length(series_list))
    series_indices <- sample(seq_len(length(series_list)), n_series)
  }
  
  # Create plotting data
  plot_data <- data.frame()
  
  for (i in series_indices) {
    series_data <- series_list[[i]]
    time_points <- seq(0, 1, length.out = nrow(series_data))
    
    for (j in seq_len(ncol(series_data))) {
      series_df <- data.frame(
        time = time_points,
        value = series_data[, j],
        dimension = paste("Lead", j),
        series_id = i,
        true_label = true_labels[i],
        cluster = cluster_assignments[i],
        true_label_name = paste("Class", true_labels[i]),
        row.names = NULL
      )
      
      plot_data <- rbind(plot_data, series_df)
    }
  }
  
  p <- ggplot(plot_data, aes(x = .data$time, y = .data$value, color = .data$true_label_name, group = .data$series_id)) +
    geom_line(alpha = 0.6, linewidth = 0.5) +
    facet_wrap(~ .data$dimension, scales = "free_y", ncol = 3) +
    labs(title = title,
         x = "Time", y = "Amplitude", color = "True Class") +
    scale_color_manual(values = c("Class 0" = "blue", "Class 1" = "red")) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

# Main function
main <- function() {
  cat("=== PTB Dataset Analysis with WICMAD ===\n")
  
  # Set data directory
  data_dir <- "data/ptb-xl-a-large-publicly-available-electrocardiography-dataset-1.0.3"
  
  # Load data
  cat("1. Loading PTB dataset...\n")
  data <- load_ptb_data(data_dir, reveal_ratio = 0.15)
  
  # Print data summary
  cat("Training data: ", length(data$train_series), " series\n")
  cat("Label distribution:\n")
  print(table(data$train_labels))
  cat("First few labels:", head(data$train_labels, 10), "\n")
  
  # Plot original data (overlapped) - subsample for visualization
  cat("\n2. Creating original data visualization...\n")
  
  # Create a subsample of 100 curves with 5% anomalies for plotting
  n_plot_samples <- 100
  n_plot_anomalies <- 5  # 5% of 100
  n_plot_normal <- 95    # 95% of 100
  
  # Find normal and abnormal indices
  normal_indices <- which(data$train_labels == 0)
  abnormal_indices <- which(data$train_labels == 1)
  
  # Sample for plotting
  plot_normal_indices <- sample(normal_indices, min(n_plot_normal, length(normal_indices)))
  plot_abnormal_indices <- sample(abnormal_indices, min(n_plot_anomalies, length(abnormal_indices)))
  plot_indices <- c(plot_normal_indices, plot_abnormal_indices)
  
  # Shuffle the indices
  plot_indices <- sample(plot_indices)
  
  # Create subsampled data for plotting
  plot_series <- data$train_series[plot_indices]
  plot_labels <- data$train_labels[plot_indices]
  
  cat("Plotting subsample: 100 curves (95 normal, 5 abnormal)\n")
  
  original_plot <- plot_overlapped_data(plot_series, plot_labels, 
                                       "PTB Dataset - Original ECG Data (Subsample)")
  
  # Save original data plot
  pdf("ptb_original_data.pdf", width = 12, height = 10)
  print(original_plot)
  dev.off()
  cat("Original data plot saved to ptb_original_data.pdf\n")
  
  # Prepare data for WICMAD
  cat("\n3. Preparing data for WICMAD...\n")
  wicmad_data <- prepare_wicmad_data(data$train_series, data$train_labels)
  
  cat("Number of series:", length(wicmad_data$Y), "\n")
  cat("Time points:", length(wicmad_data$t), "\n")
  cat("Series dimensions:", nrow(wicmad_data$Y[[1]]), "x", ncol(wicmad_data$Y[[1]]), "\n")
  
  # Run WICMAD clustering
  cat("\n4. Running WICMAD clustering...\n")
  wicmad_result <- wicmad(
    Y = wicmad_data$Y,
    t = wicmad_data$t,
    n_iter = 10000,
    burn = 3000,
    thin = 1,
    warmup_iters = 500
  )
  
  # Extract cluster assignments
  cluster_assignments <- wicmad_result$Z[nrow(wicmad_result$Z), ]
  
  # Map WICMAD clusters to binary classification (normal vs abnormal)
  # Find the cluster that contains the most normal samples (label 0)
  normal_samples <- which(wicmad_data$labels == 0)
  cluster_counts <- table(cluster_assignments[normal_samples])
  normal_cluster <- as.numeric(names(cluster_counts)[which.max(cluster_counts)])
  
  # Map to binary: normal cluster = 0, all others = 1 (abnormal)
  binary_cluster_assignments <- ifelse(cluster_assignments == normal_cluster, 0, 1)
  
  cat("WICMAD found", length(unique(cluster_assignments)), "clusters\n")
  cat("Normal cluster:", normal_cluster, "\n")
  cat("Binary mapping: cluster", normal_cluster, "-> normal (0), others -> abnormal (1)\n")
  
  # Calculate metrics
  cat("\n5. Calculating clustering metrics...\n")
  metrics <- calculate_clustering_metrics(wicmad_data$labels, binary_cluster_assignments)
  
  cat("Clustering Performance:\n")
  cat("Macro Precision:", round(metrics$Macro_Precision, 4), "\n")
  cat("Macro Recall:", round(metrics$Macro_Recall, 4), "\n")
  cat("Macro Specificity:", round(metrics$Macro_Specificity, 4), "\n")
  cat("Macro F1 Score:", round(metrics$Macro_F1, 4), "\n")
  cat("Overall Accuracy:", round(metrics$Overall_Accuracy, 4), "\n")
  cat("Adjusted Rand Index:", round(metrics$ARI, 4), "\n")
  
  # Plot clustering results - use same subsample for consistency
  cat("\n6. Creating clustering results visualization...\n")
  
  # Use the same subsample indices for clustering results
  plot_cluster_series <- data$train_series[plot_indices]
  plot_cluster_labels <- wicmad_data$labels[plot_indices]
  plot_cluster_assignments <- binary_cluster_assignments[plot_indices]
  
  clustering_plot <- plot_clustering_results(
    plot_cluster_series, 
    plot_cluster_labels, 
    plot_cluster_assignments,
    "PTB - Clustering Results (Subsample)"
  )
  
  # Save clustering results plot
  pdf("ptb_clustering_results.pdf", width = 12, height = 10)
  print(clustering_plot)
  dev.off()
  cat("Clustering results plot saved to ptb_clustering_results.pdf\n")
  
  cat("\n=== Analysis Complete ===\n")
  cat("Generated files:\n")
  cat("- ptb_original_data.pdf: Original ECG data (overlapped)\n")
  cat("- ptb_clustering_results.pdf: Clustered ECG data\n")
  cat("\nFinal Performance:\n")
  cat("Macro Precision:", round(metrics$Macro_Precision, 4), "\n")
  cat("Macro Recall:", round(metrics$Macro_Recall, 4), "\n")
  cat("Macro F1 Score:", round(metrics$Macro_F1, 4), "\n")
  cat("Overall Accuracy:", round(metrics$Overall_Accuracy, 4), "\n")
  cat("Adjusted Rand Index:", round(metrics$ARI, 4), "\n")
}

# Run the analysis
main()
