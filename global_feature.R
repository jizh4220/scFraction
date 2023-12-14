# source ~/conda/bin/activate BPCells
# Rscript global_Feature.R goi.txt user_id
args <- commandArgs(trailingOnly = TRUE)
# test if there is at least one argument: if not, return an error
if (length(args) != 2) {
  stop("Two arguments must be supplied \n", call. = FALSE)
}
# img_id <- "LUAD;whole_lung;GSM5702473"

user_id <- args[2]
library(ggplot2)
library(ggdendro)
library(cowplot)
library(ggtree)
library(aplot)
library(tidyverse)
library(patchwork)
library(Seurat)
library(BPCells)
`%!in%` <- Negate(`%in%`)


print("=====================================")
print("=====         Global Gene Explorer         =====")
print("=====================================")


setwd(file.path(user_id))
goi <- suppressWarnings(as.character(read.table(args[1], header = FALSE, sep = ";")))

output_dir <- "global_Feature/"
suppressWarnings(dir.create(output_dir, recursive = TRUE))


fread_counts <- function(path, header = FALSE) {
    library(data.table)
    if (header == TRUE) {
        counts <- as.data.frame(fread(path, header = TRUE))
    } else {
        counts <- as.data.frame(fread(path))
    }
    return(counts)
}

MG_folders <- list.files("~/PulmoScope/global/MGDB",
        pattern = "val",
        recursive = TRUE, full.names = TRUE)
MG_folders <- fs::path_dir(MG_folders)
global_meta <- fread_counts("~/PulmoScope/global/global_meta_fourth_clean.csv")
rownames(global_meta) <- global_meta$idx

# load necessary libraries

BPCells_global_feature <- function(MG_folders, goi, output_dir) {
        mg_goi <- open_matrix_dir(MG_folders[1])[goi, ]
        for (i in 2:length(MG_folders)) {
                cur_mg <- open_matrix_dir(MG_folders[i])
                mg_goi <- methods::cbind2(mg_goi, cur_mg[goi, ])
        }
        # mg_goi <- mg_goi[Matrix::rowSums(mg_goi) > 3, ]
        # # Normalize by reads-per-cell
        # mg_goi <- multiply_cols(mg_goi, 1/Matrix::colSums(mg_goi))
        # # Log normalization
        # mg_goi <- log1p(mg_goi * 10000) # Log normalization
        global_meta$groups_ct <- paste(global_meta$groups, global_meta$celltype, sep = ";")
        groups_ct <- unique(global_meta$groups_ct)
        ps <- list()
        # group_celltype wise traversal
        for (j in 1:length(groups_ct)) {
                cur_idx <- rownames(global_meta[global_meta$groups_ct %in% groups_ct[j], ])
                idx_mg <- mg_goi[, intersect(colnames(mg_goi), cur_idx)]
                gct_goi <- rowMeans(expm1(idx_mg))
                # Calculate the percentage of cells with expression above the threshold
                pct_mg <- rowSums(idx_mg > 0) / ncol(idx_mg) * 100
                # store
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
        # Save the data frame to a CSV file
        write.csv(gene_cluster, paste0(output_dir, "global_features_counts.csv"))
        return(gene_cluster)
}


BPCells_global_feature <- function(gene_cluster, output_dir) {
        goi <- unique(gene_cluster$gene)
        for (i in 1:length(goi)) {
                goi_df <- gene_cluster[gene_cluster$gene == goi[i], ]
                
                rownames(goi_df) <- goi_df$group_celltype
                mat <- goi_df %>% 
                        select(-pct.exp) %>%  # drop unused columns to faciliate widening
                        pivot_wider(names_from = celltype, values_from = avg.exp) %>% 
                        mutate_all(~ifelse(is.infinite(.), 500, .)) %>% 
                        mutate_all(~ifelse(is.na(.), 0, .)) %>%  # replace Inf with largest finite number
                        data.frame()

                mat_dist <- dist((mat[, 4:ncol(mat)]  %>% as.matrix() %>% t()))
                groups_clust <- hclust(mat_dist)
                ct_order <- groups_clust$labels[groups_clust$order]

                ddgram_col <- as.dendrogram(groups_clust)
                # Create the ggtree plot
                ggtree_plot_col <- ggtree(ddgram_col) +
                                # geom_tiplab(size = 2) + 
                                layout_dendrogram() +
                                theme(legend.position='none') +
                                xlim2(dotplot)

                goi_df$celltype <- gsub("\\+|\\ ", "\\.", goi_df$celltype)
                dotplot <- goi_df %>% 
                        mutate(`% Expressing` = pct.exp,
                                celltype = factor(celltype, levels = ct_order)) %>% 
                        filter(avg.exp != 0, `% Expressing` > 0) %>% 
                        
                        ggplot(aes(x = celltype, y = groups, color = avg.exp, size = `% Expressing`)) + 
                        geom_point() + 
                        cowplot::theme_cowplot() + 
                        # theme(axis.line  = element_blank()) +
                        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
                                panel.grid.major = element_line(colour = "grey"),
                                panel.grid.minor = element_line(colour = "grey"),
                                panel.background = element_rect(fill = "white"),) +
                        xlab('') +
                        ylab('') +
                        scale_color_gradientn(colours = viridis::viridis(20), limits = c(0, 5), oob = scales::squish, name = 'log1p (avgExp)')

                p <-    ggtree_plot_col +
                        dotplot + 
                        plot_layout(ncol = 1, widths = c(0.1, 0.5), heights = c(1, 5)) +
                        plot_annotation(title = paste0("Global Feature Dot Plot: ",
                                                goi[i]),
                                        theme = theme(plot.title = element_text(size = 20, hjust = 0.5)))
                # p <- plot_grid(dotplot, nrow = 1, rel_widths = c(0.5,2), align = 'h')
                p_title <- paste0(output_dir, goi[i], "_global_feature.png")
                ggsave(filename = p_title,
                        plot = p, width = 12, height = 8, units = "in")
        }

        # Save the data frame to a CSV file
        write.csv(gene_cluster, "/home/lijia/zhangjiaxuan/PulmoScope/global_feature/test_5m_global_features.csv")
}

gene_cluster <- BPCells_global_feature(MG_folders, goi, output_dir) 
BPCells_global_feature(gene_cluster, output_dir = "global_Feature/")
