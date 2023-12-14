
library(dplyr)
library(data.table)
library(ggplot2)
library(ggdendro)
library(cowplot)
library(ggtree)
library(aplot)
library(tidyverse)
library(patchwork)
library(Seurat)
library(BPCells)


human2murine_gene <- function(mouse_human_genes, goi) {
    goi_idx <- mouse_human_genes[mouse_human_genes$Symbol %in% goi, ]$DB.Class.Key
    murine_class <- subset(mouse_human_genes, Common.Organism.Name == "mouse, laboratory")
    murine_name <- subset(murine_class, DB.Class.Key %in% goi_idx)
    murine_goi <- murine_name$Symbol
    return (murine_goi)
}



goi_mgBP <- function(MG_folders, goi) {
    mg_goi <- open_matrix_dir(MG_folders[1])[goi, ]
    for (i in 2:length(MG_folders)) {
            cur_mg <- open_matrix_dir(MG_folders[i])
            mg_goi <- methods::cbind2(mg_goi, cur_mg[goi, ])
    }
    return(mg_goi)
}

refine_metanames <- function(meta) {
        colnames(meta)[colnames(meta) == "orig.ident"] <- "accession"
        colnames(meta)[colnames(meta) == "gse_alias"] <- "experiment"
        colnames(meta)[colnames(meta) == "DISEASE"] <- "disease"
        colnames(meta)[colnames(meta) == "TISSUE"] <- "tissue"
        colnames(meta)[colnames(meta) == "celltype"] <- "detail_celltype"
        colnames(meta)[colnames(meta) == "macro_celltype"] <- "celltype"
        return(meta)
}

# This function calculates the fraction of cells with expression above a certain threshold for each gene of interest (goi) in each combination of accession and cell type. The results are returned as a data frame.
# 
# Args:
#   mg_goi: A matrix with gene expression data. Rows represent genes and columns represent cells.
#   mg_meta: A data frame with metadata for each cell. It should include 'accession' and 'celltype' columns.
#   goi: A vector of genes of interest.
#   threshold: A numeric value indicating the expression threshold. Cells with expression above this threshold are considered as expressing the gene. Default is 0.5.
# 
# Returns:
#   A data frame with the fraction of cells expressing the gene ('pct.exp') for each gene of interest in each combination of accession and cell type. The data frame also includes 'accession', 'celltype', 'gene', and 'num_cells' columns.

fraction_accCT_table <- function(mg_goi, mg_meta, goi, threshold = 0.5) {
  # Check if all genes of interest are in the gene expression data
  if (length(intersect(goi, rownames(mg_goi))) == 0) {
    stop("Check your goi content!\n")
  }
  
  # Get the unique accessions and cell types
  accessions <- unique(mg_meta$accession)
  celltypes <- unique(mg_meta$celltype)
  
  # Initialize a data frame to store the results
  df <- data.frame(matrix(ncol = length(goi), nrow = length(accessions) * length(celltypes)))
  colnames(df) <- goi
  rownames(df) <- paste(rep(accessions, each = length(celltypes)), rep(celltypes, times = length(accessions)), sep = "_")
  
  # Concatenate 'accession' and 'celltype' columns to create a new 'acc_ct' column
  mg_meta$acc_ct <- paste(mg_meta$accession, mg_meta$celltype, sep = ";")
  
  # Traverse each accession
  for (j in 1:length(accessions)) {
    # Subset the metadata for the current accession
    cur_acc <- mg_meta[mg_meta$accession %in% accessions[j], ]
    
    # Get the unique cell types for the current accession
    all_ct <- unique(cur_acc$celltype)
    
    # Traverse each cell type
    for (z in 1:length(all_ct)) {
      # Subset the metadata for the current cell type
      ct_acc <- cur_acc[cur_acc$celltype %in% all_ct[z], ]
      
      # Subset the gene expression data for the cells in the current cell type
      idx_mg <- mg_goi[, intersect(colnames(mg_goi), rownames(ct_acc))]
      
      # Calculate the fraction of cells with expression above the threshold for each gene of interest and store in the data frame
      df[paste(accessions[j], all_ct[z], sep = "_"), ] <- rowSums(idx_mg >= threshold) / ncol(idx_mg)
    }
  }
  
  # Calculate the number of cells for each combination of accession and cell type
  num_cells <- as.data.frame(table(mg_meta$accession, mg_meta$celltype))
  rownames(num_cells) <- paste(num_cells$Var1, num_cells$Var2, sep = "_")
  num_cells <- num_cells[match(rownames(df), rownames(num_cells)), ]
  df$num_cells <- num_cells$Freq
  df$acc_ct <- rownames(df)
  
  # Merge the data frame with the metadata
  metadata <- mg_meta %>%
    group_by(accession, celltype) %>%
    summarise_all(first) %>%
    mutate(acc_ct = paste(accession, celltype, sep = "_"))
  frac_table <- merge(df, metadata, by = "acc_ct", all.x = TRUE)
  frac_table <- frac_table[frac_table$num_cells > 0, ]
  
  # Return the data frame
  return(frac_table)
}


# This function calculates the average expression and the percentage of cells with expression above a certain threshold for each gene of interest (goi) in each cell type group. The results are returned as a data frame.
# 
# Args:
#   mg_goi: A matrix with gene expression data. Rows represent genes and columns represent cells.
#   global_meta: A data frame with metadata for each cell. It should include 'groups' and 'celltype' columns.
#   goi: A vector of genes of interest.
#   threshold: A numeric value indicating the expression threshold. Cells with expression above this threshold are considered as expressing the gene. Default is 0.5.
# 
# Returns:
#   A data frame with the average expression ('avg.exp') and the percentage of cells expressing the gene ('pct.exp') for each gene of interest in each cell type group. The data frame also includes 'groups', 'celltype', and 'gene' columns.

BPCells_global_feature <- function(mg_goi, global_meta, goi, threshold = 0.5) {
  # Concatenate 'groups' and 'celltype' columns to create a new 'groups_ct' column
  global_meta$groups_ct <- paste(global_meta$groups, global_meta$celltype, sep = ";")
  
  # Get the unique cell type groups
  groups_ct <- unique(global_meta$groups_ct)
  
  # Initialize a list to store the results for each cell type group
  ps <- list()
  
  # Traverse each cell type group
  for (j in 1:length(groups_ct)) {
    # Get the indices of the cells in the current cell type group
    cur_idx <- rownames(global_meta[global_meta$groups_ct %in% groups_ct[j], ])
    
    # Subset the gene expression data for the cells in the current cell type group
    idx_mg <- mg_goi[, intersect(colnames(mg_goi), cur_idx)]
    
    # Calculate the average expression for each gene of interest
    gct_goi <- rowMeans(expm1(idx_mg))
    
    # Calculate the percentage of cells with expression above the threshold for each gene of interest
    pct_mg <- rowSums(idx_mg > threshold) / ncol(idx_mg)
    
    # Store the results in the list
    ps[[j]] <- cbind(avg.exp = log1p(gct_goi), pct.exp = pct_mg, group_celltype = groups_ct[j], gene = goi)
  }
  
  # Combine all the results into a single data frame
  gene_cluster <- do.call(rbind, ps)
  
  # Convert the result to a data frame and create new columns
  gene_cluster <- as.data.frame(gene_cluster)
  gene_cluster$groups <- gsub(";.*", "", gene_cluster$group_celltype)
  gene_cluster$celltype <- gsub(".*;", "", gene_cluster$group_celltype)
  gene_cluster$avg.exp <- as.numeric(gene_cluster$avg.exp)
  gene_cluster$pct.exp <- as.numeric(gene_cluster$pct.exp)
  
  # Return the data frame
  return(gene_cluster)
}



# This function creates a dot plot for each gene of interest (goi) showing the average expression and the percentage of cells expressing the gene in each cell type group. The dot plots are saved as PNG files in the specified output directory.
# 
# Args:
#   gene_cluster: A data frame with the average expression ('avg.exp') and the percentage of cells expressing the gene ('pct.exp') for each gene of interest in each cell type group. The data frame should also include 'groups', 'celltype', and 'gene' columns.
#   goi: A vector of genes of interest.
#   output_dir: A string specifying the directory where the dot plots should be saved.
# 
# Returns:
#   NULL. The function saves the dot plots as PNG files in the specified output directory.

feature_dotplot <- function(gene_cluster, goi, output_dir) {
  # Traverse each gene of interest
  for (i in 1:length(goi)) {
    # Subset the data for the current gene of interest
    goi_df <- gene_cluster[gene_cluster$gene == goi[i], ]
    
    # Set the row names to be the 'group_celltype' column
    rownames(goi_df) <- goi_df$group_celltype
    
    # Pivot the data frame to a wider format and replace infinite and NA values with 0
    mat <- goi_df %>% 
      select(-pct.exp) %>%
      pivot_wider(names_from = celltype, values_from = avg.exp) %>%
      mutate_all(~ifelse(is.infinite(.), 500, .)) %>%
      mutate_all(~ifelse(is.na(.), 0, .)) %>%
      data.frame()
    
    # Calculate the distance matrix and perform hierarchical clustering
    mat_dist <- dist((mat[, 4:ncol(mat)]  %>% as.matrix() %>% t()))
    groups_clust <- hclust(mat_dist)
    
    # Get the order of the cell types from the clustering result
    ct_order <- groups_clust$labels[groups_clust$order]
    
    # Convert the clustering result to a dendrogram
    ddgram_col <- as.dendrogram(groups_clust)
    
    # Replace '+' and ' ' with '.' in the 'celltype' column
    goi_df$celltype <- gsub("\\+|\\ ", "\\.", goi_df$celltype)
    
    # Create the dot plot
    dotplot <- goi_df %>% 
      mutate(`% Expressing` = pct.exp,
             celltype = factor(celltype, Levels = ct_order)) %>% 
      filter(avg.exp != 0, `% Expressing` > 0) %>% 
      ggplot(aes(x = celltype, y = groups, color = avg.exp, size = `% Expressing`)) + 
      geom_point() + 
      cowplot::theme_cowplot() +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
            panel.grid.major = element_line(colour = "grey"),
            panel.grid.minor = element_line(colour = "grey"),
            panel.background = element_rect(fill = "white"),) +
      xlab('') +
      ylab('') +
      scale_color_gradientn(colours = viridis::viridis(20), limits = c(0, 5), oob = scales::squish, name = 'log1p (avgExp)')
    
    # Create the ggtree plot
    ggtree_plot_col <- ggtree(ddgram_col) +
      layout_dendrogram() +
      theme(legend.position='none') +
      xlim2(dotplot)
    
    # Combine the ggtree plot and the dot plot
    p <- ggtree_plot_col +
      dotplot + 
      plot_layout(ncol = 1, widths = c(0.1, 0.5), heights = c(1, 5)) +
      plot_annotation(title = paste0("Global Feature Dot Plot: ", goi[i]),
                      theme = theme(plot.title = element_text(size = 20, hjust = 0.5)))
    
    # Save the plot as a PNG file in the output directory
    p_title <- paste0(output_dir, goi[i], "_global_feature.png")
    ggsave(filename = p_title,
           plot = p, width = 12, height = 8, units = "in")
  }
  
  # Return NULL
  return(NULL)
}


# This function creates a Dot plot for each gene of interest (goi) showing the mean expression fraction and the number of cells expressing the gene in each cell type. The results are returned as a ggplot object.
# 
# Args:
#   frac_table: A data frame with the number of cells ('num_cells') and the expression data for each gene of interest in each cell type. The data frame should also include 'celltype' and 'gene' columns.
#   goi: A vector of genes of interest.
# 
# Returns:
#   A ggplot object with the Dot plot.

gene_ct_fraction_plot <- function(frac_table, goi) {
  # Calculate the mean expression fraction for each gene of interest in each cell type
  frac_mean <- frac_table %>%
    mutate(across(goi, ~replace(., is.na(.), 0))) %>%
    group_by(celltype) %>%
    summarise(across(goi, mean, na.rm = TRUE))
  
  # Pivot the data frame to a longer format
  longer_frac_mean <- frac_mean %>%
    pivot_longer(cols = -celltype, names_to = "gene", values_to = "exp_frac")
  
  # Exclude rows with NA 'celltype'
  longer_frac_mean <- longer_frac_mean[!is.na(longer_frac_mean$celltype), ]
  
  # Calculate the total number of cells for each cell type
  exp_cell_table <- frac_table %>%
    group_by(celltype) %>%
    summarise(exp_cells = sum(num_cells, na.rm = TRUE))
  
  # Merge the data frames
  longer_frac_mean <- merge(longer_frac_mean, exp_cell_table, by = "celltype", all.x = TRUE)
  
  # Calculate the total number of cells expressing each gene for each cell type
  longer_frac_mean$exp_cells <- longer_frac_mean$exp_cells * longer_frac_mean$exp_frac
  
  # Create a Dot plot
  p <- ggplot(longer_frac_mean, aes(x = gene, y = celltype)) +
    geom_point(aes(size = exp_cells, color = exp_frac, alpha = 0.8)) +
    scale_color_gradient(low = "white", high = "red", limits = c(0, 0.5)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = element_text(size = 6)) +
    labs(color = "Mean Fraction") +
    coord_flip()
  
  # Return the ggplot object
  return(p)
}


# This function creates a Dot plot for a specific cell type of interest (ctoi) showing the expression fraction and the number of cells expressing each gene of interest (goi) in each accession. The results are returned as a ggplot object.
# 
# Args:
#   frac_table: A data frame with the number of cells ('num_cells') and the expression data for each gene of interest in each cell type. The data frame should also include 'accession', 'celltype', and 'gene' columns.
#   goi: A vector of genes of interest.
#   ctoi: A string specifying the cell type of interest. Default is "Astro".
# 
# Returns:
#   A ggplot object with the fraction Dotplot plot.

acc_ct_dotplot <- function(frac_table, goi, ctoi = "Astro") {
  # Subset the data for the cell type of interest
  ct_frac <- subset(frac_table, celltype == ctoi)
  
  # Pivot the data frame to a longer format and replace NA values with 0
  longer_ct_frac <- ct_frac %>%
    mutate(across(all_of(goi), ~replace(., is.na(.), 0))) %>%
    group_by(accession) %>%
    reframe(across(all_of(goi))) %>%
    pivot_longer(cols = -accession, names_to = "gene", values_to = "exp_frac")
  
  # Calculate the total number of cells for each accession
  exp_cell_table <- ct_frac %>%
    group_by(accession) %>%
    summarise(exp_cells = sum(num_cells, na.rm = TRUE))
  
  # Merge the data frames
  longer_ct_frac <- merge(longer_ct_frac, exp_cell_table, by = "accession", all.x = TRUE)
  
  # Calculate the total number of cells expressing each gene for each accession
  longer_ct_frac$exp_cells <- longer_ct_frac$exp_cells * longer_ct_frac$exp_frac
  
  # Create a Dot plot
  p <- ggplot(longer_ct_frac, aes(x = gene, y = exp_frac)) +
    geom_point(aes(size = exp_cells, color = exp_frac, alpha = 0.2)) +
    scale_color_gradient(low = "white", high = "red", limits = c(0.5, 1)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = element_text(size = 6)) +
    labs(color = "Accession Fraction")
  
  # Return the ggplot object
  return(p)
}


# This function creates a Dot plot for each gene of interest (goi) showing the differential expression fraction and the difference in the number of cells expressing the gene between two diseases in each cell type. The results are returned as a ggplot object.
# 
# Args:
#   frac_table: A data frame with the number of cells ('num_cells'), the disease, and the expression data for each gene of interest in each cell type. The data frame should also include 'accession', 'celltype', and 'gene' columns.
#   goi: A vector of genes of interest.
#   groupby: A string specifying the column to group by. Default is 'celltype'.
# 
# Returns:
#   A ggplot object with the Dot plot.

acc_ct_Diff_frac <- function(frac_table, goi, groupby='celltype') {
  # Get the unique diseases
  diseases <- unique(frac_table$disease[!is.na(frac_table$disease)])
  
  # Subset the data for each disease
  df_disease1 <- frac_table[frac_table$disease == diseases[1], ]
  df_disease2 <- frac_table[frac_table$disease == diseases[2], ]
  
  # Get the cell types shared by the two diseases
  shared_ct <- intersect(df_disease1$celltype, df_disease2$celltype)
  
  # Subset the data for the shared cell types
  df_disease1 <- df_disease1[match(shared_ct, df_disease1$celltype), ]
  df_disease2 <- df_disease2[match(shared_ct, df_disease2$celltype), ]
  
  # Calculate the mean expression fraction for each gene of interest in each cell type for each disease
  df_mean1 <- fraction_mean(df_disease1, goi)
  df_mean2 <- fraction_mean(df_disease2, goi)

  # Calculate the ratio of the mean expression fractions
  df_diff <- df_mean1[, -1] / df_mean2[, -1]
  df_diff$celltype <- shared_ct
  
  # Pivot the data frame to a longer format and replace infinite and NA values with 0
  df_diff_long <- df_diff %>%
    pivot_longer(cols = -celltype, names_to = "gene", values_to = "mean_DEfraction") %>%
    mutate(across(mean_DEfraction, ~replace(., is.na(.)|is.infinite(.), 0)))
  
  # Calculate the total number of cells for each cell type for each disease
  exp_cell_1 <- df_disease1 %>%
    group_by(celltype) %>%
    summarise(exp_cells = sum(num_cells, na.rm = TRUE))
  exp_cell_2 <- df_disease2 %>%
    group_by(celltype) %>%
    summarise(exp_cells = sum(num_cells, na.rm = TRUE))

  # Calculate the difference in the total number of cells between the two diseases
  exp_diff <- exp_cell_1[, 2] - exp_cell_2[, 2]
  exp_diff$celltype <- shared_ct
  
  # Merge the data frames
  df_diff_long <- merge(df_diff_long, exp_diff, by = "celltype", all.x = TRUE)
  
  # Calculate the total number of cells expressing each gene for each cell type
  df_diff_long$exp_cells <- df_diff_long$exp_cells * df_diff_long$mean_DEfraction
  
  # Create a Dot plot
  p <- ggplot(df_diff_long, aes(x = gene, y = celltype)) +
    geom_point(aes(size = abs(exp_cells), color = mean_DEfraction), alpha = 1) +
    scale_color_gradient(low = "white", high = "red", limits = c(0, 10)) +
    scale_size(limits = c(100, 1000)) +  # Adjust the range as needed
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),
          axis.text.y = element_text(size = 12)) +
    labs(color = "Differential Fraction")
  
  # Return the ggplot object
  return(p)
}



# This function creates a bar plot for each gene of interest (goi) showing the expression fraction in each cell type for a specific parameter (e.g., disease). The results are returned as a list of ggplot objects.
# 
# Args:
#   frac_table: A data frame with the number of cells ('num_cells'), the parameter of interest, and the expression data for each gene of interest in each cell type. The data frame should also include 'accession', 'celltype', and 'gene' columns.
#   goi: A vector of genes of interest. Default is 'goi'.
#   threshold: A numeric value indicating the expression threshold. Cells with expression above this threshold are considered as expressing the gene. Default is 0.5.
#   param: A string specifying the parameter of interest. Default is 'disease'.
#   floor: A numeric value indicating the minimum number of cells for a cell type to be included in the plot. Default is 50.
# 
# Returns:
#   A list of ggplot objects with the bar plots.

fraction_Barplot <- function(frac_table, goi = goi,
                          threshold = 0.5,
                          param = "disease",
                          floor = 50
                          ) {
    frac_table$param <- frac_table[, grep(param, colnames(frac_table))]
    sample_names <- names(table(frac_table$param))
    # frac_table$ct_disease <- paste0(frac_table$celltype,
    #               "+", frac_table$disease)
    ps <- list()
    for (i in 1:length(goi)) {
        cur_frac <- frac_table[frac_table$gene == goi[i], ]
        all_sample_table <- c(length(unique(cur_frac$accession[cur_frac$param == sample_names[1]])),
                          length(unique(cur_frac$accession[cur_frac$param == sample_names[2]])))
        cur_frac <- cur_frac[cur_frac$cell_counts >= floor &
                                cur_frac$expressed_ratio > 0, ]
        if (nrow(cur_frac) == 0) {
          cat("Current gene:",
            goi[i],
            "is not found in Fraction Table \n")
          next
        }
        cur_frac$param_ct <- paste0(cur_frac$celltype, "_", cur_frac$param)
        sample_table <- c(length(unique(cur_frac$accession[cur_frac$param == sample_names[1]])),
                          length(unique(cur_frac$accession[cur_frac$param == sample_names[2]])))
        p <- ggplot(cur_frac,
                    aes(x = param_ct, y = expressed_ratio,
                    fill = param)) +
                    geom_boxplot() +
                    # geom_violin(trim = FALSE) +
                    geom_point(
                          position = position_jitter(
                          seed = 1, width = 0.2, height = 0),
                          alpha = 2,
                          size = 1.5) +
                    # geom_boxplot(width = 0.05, fill = "white") +
                    # facet_wrap(~gene) +
                    labs(x = "Genes of Interest",
                        y = "Fraction of Cells Expressing Target Genes",
                        size = 100) +
                theme(axis.text = element_text(size = 12, angle = 90),
                      axis.text.x = element_text(vjust = 0.5, angle = 90),
                      legend.key.size = unit(1.5, "cm"),
                      legend.text = element_text(size = 10),
                      plot.title = element_text(size = 20, hjust = 0.5),
                      panel.grid.major = element_line(colour = "grey"),
                      panel.grid.minor = element_line(colour = "grey"),
                      panel.background = element_rect(fill = "white"),
                      plot.subtitle = element_text(size = 14))
                # theme(panel.background = element_rect(fill = "white"))
        p <- p + labs(title = paste0("Fraction plot of GOI: ",
                                    goi[i]),
                      subtitle = paste0(
                                "Total ",
                                sample_names[1],
                                " Sample Number: ",
                                all_sample_table[1],
                                "\nTotal ",
                                sample_names[1],
                                " Sample Number Pass Cutoff: ",
                                sample_table[1],
                                "\nTotal ",
                                sample_names[2],
                                " Sample Number: ",
                                all_sample_table[2],
                                "\nTotal ",
                                sample_names[2],
                                " Sample Number Pass Cutoff: ",
                                sample_table[2],
                                # length(unique(frac_table$cell_counts)),
                                "\nTotal Cell Counts: ",
                                sum(unique(frac_table$cell_counts)),
                                "\nPass Cutoff Cell Counts: ",
                                round(sum(cur_frac$cell_counts *
                                  cur_frac$expressed_ratio), digits = 0)
                                ,
                                "\nExpression Cutoff: ", threshold,
                                "\nMinimum Cell Number: ", floor
                                ),
                     # tag = paste0("Normal Sample Numbe: ")
                    )
        ps[[i]] <- p
    }
    return(ps)
    cat("Successfully Plot Fraction Plot\n")
}


# This function creates a scatter plot for each gene of interest (goi) showing the correlation of expression fractions between two parameters (e.g., diseases) in each cell type. The results are returned as a list of ggplot objects.
# 
# Args:
#   frac_table: A data frame with the number of cells ('num_cells'), the parameter of interest, and the expression data for each gene of interest in each cell type. The data frame should also include 'accession', 'celltype', and 'gene' columns.
#   goi: A vector of genes of interest.
#   threshold: A numeric value indicating the expression threshold. Cells with expression above this threshold are considered as expressing the gene. Default is 0.
#   param: A string specifying the parameter of interest. Default is 'param'.
#   floor: A numeric value indicating the minimum number of cells for a cell type to be included in the plot. Default is 10.
# 
# Returns:
#   A list of ggplot objects with the scatter plots.

corr_fraction_plot <- function(frac_table, goi,
                          threshold = 0,
                          param = "param",
                          floor = 10
                          ) {
    frac_table$param <- frac_table[, grep(param, colnames(frac_table))]
    sample_names <- names(table(frac_table$param))
    # frac_table$ct_disease <- paste0(frac_table$celltype,
    #               "+", frac_table$disease)
    ps <- list()
    for (i in 1:length(goi)) {
        cur_frac <- frac_table[frac_table$gene == goi[i], ]
        all_sample_table <- c(length(unique(cur_frac$accession[cur_frac$param == sample_names[1]])),
                          length(unique(cur_frac$accession[cur_frac$param == sample_names[2]])))
        cur_frac <- cur_frac[cur_frac$cell_counts >= floor &
                                cur_frac$expressed_ratio > threshold, ]
        if (nrow(cur_frac) == 0) {
          cat("Current gene:",
            goi[i],
            "is not found in Fraction Table \n")
          next
        }
        cur_frac$groups_param <- paste0(cur_frac$groups, "_", cur_frac$param)
        sample_table <- c(length(unique(cur_frac$accession[cur_frac$param == sample_names[1]])),
                          length(unique(cur_frac$accession[cur_frac$param == sample_names[2]])))
        p <- ggplot(cur_frac,
                    aes(x = groups_param, y = expressed_ratio,
                    fill = param)) +
                    geom_boxplot() +
                    # geom_violin(trim = FALSE) +
                    geom_point(
                          position = position_jitter(
                          seed = 1, width = 0.2, height = 0),
                          alpha = 2,
                          size = 1.5) +
                    # geom_boxplot(width = 0.05, fill = "white") +
                    # facet_wrap(~condition) +
                    labs(x = "Genes of Interest",
                        y = "Fraction of Cells Expressing Target Genes",
                        size = 100) +
                theme(axis.text = element_text(size = 12, angle = 45),
                        legend.key.size = unit(1.5, "cm"),
                        legend.text = element_text(size = 10),
                        plot.title = element_text(size = 20, hjust = 0.5),
                        panel.grid.major = element_line(colour = "grey"),
                        panel.grid.minor = element_line(colour = "grey"),
                        panel.background = element_rect(fill = "white"),
                        plot.subtitle = element_text(size = 14))
                # theme(panel.background = element_rect(fill = "white"))
        p <- p + labs(title = paste0("Fraction plot of GOI: ",
                                    goi[i]),
                      subtitle = paste0(
                                "Total ",
                                sample_names[1],
                                " Sample Number: ",
                                all_sample_table[1],
                                "\nTotal ",
                                sample_names[1],
                                " Sample Number Pass Cutoff: ",
                                sample_table[1],
                                "\nTotal ",
                                sample_names[2],
                                " Sample Number: ",
                                all_sample_table[2],
                                "\nTotal ",
                                sample_names[2],
                                " Sample Number Pass Cutoff: ",
                                sample_table[2],
                                # length(unique(frac_table$cell_counts)),
                                "\nTotal Cell Counts: ",
                                sum(unique(frac_table$cell_counts)),
                                "\nPass Cutoff Cell Counts: ",
                                round(sum(cur_frac$cell_counts *
                                  cur_frac$expressed_ratio), digits = 0)
                                ,
                                "\nExpression Cutoff: ", threshold,
                                "\nMinimum Cell Number: ", floor
                                ),
                     # tag = paste0("Normal Sample Numbe: ")
                    )
        ps[[i]] <- p
    }
    return(ps)
    cat("Successful Cor-Fraction Plot\n")
}


# Convert the data frame to a character vector
goi <- as.character(unlist(read.table("ECM_genes.txt", header = FALSE, sep = "\t")))
mouse_human_genes <- read.csv("HOM_MouseHumanSequence.rpt", sep="\t")
murine_goi <- human2murine_gene(mouse_human_genes, goi)
# prepare the count matrix
mg_goi <- counts[intersect(unique(murine_goi), rownames(counts)), ]

meta <- refine_metanames(meta)
meta$celltype <- meta$detail_celltype


frac_table <- fraction_accCT_table(mg_goi, mg_meta = meta, goi = rownames(mg_goi), threshold = 0.5)

# Level0
meta$groups <- meta$disease
gene_cluster <- BPCells_global_feature(mg_goi, global_meta = meta, goi = rownames(mg_goi), threshold = 0.5)
feature_dotplot(gene_cluster, goi = c("Ddr2", "Lamc2"), output_dir)


# Level1
p <- gene_ct_fraction_plot(frac_table, goi = rownames(mg_goi))
p_title <- paste0(output_dir, "global_acc_ct_dotplot.png")
ggsave(filename = p_title,
    plot = p, width = 8, height = 18, units = "in")


# Level2
ctoi <- "Astro"
p <- acc_ct_dotplot(tmp, goi, ctoi = ctoi)
p_title <- paste0(output_dir, ctoi, "_ct_acc_dotplot.png")
ggsave(filename = p_title,
    plot = p, width = 20, height = 8, units = "in")


# Level3
p <- acc_ct_Diff_frac(frac_table, goi, groupby = 'celltype')
p_title <- paste0(output_dir, "def_dotplot.png")
ggsave(filename = p_title,
    plot = p, width = 28, height = 8, units = "in")
