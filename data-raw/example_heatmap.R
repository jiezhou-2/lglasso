# Create sample data ===================================================
library(RColorBrewer)
set.seed(43)
data <- matrix(rnorm(500), 50, 10)
colnames(data) <- paste0("Sample_", 1:10)
rownames(data) <- paste0("Gene_", 1:50)

head(data)

# Annotations ===================================================

# create a data frame for column annotation
ann_df <- data.frame(Group = rep(c("Disease", "Control"), c(5, 5)),
                     Lymphocyte_count = rnorm(10))
row.names(ann_df) <- colnames(data)
head(ann_df)

gene_functions_df <- data.frame(gene_functions = rep(c('Oxidative_phosphorylation',
                                                       'Cell_cycle',
                                                       'Immune_regulation',
                                                       'Signal_transduction',
                                                       'Transcription'), rep(10, 5)))
row.names(gene_functions_df) <- rownames(data)

ann_colors <- list(
  gene_functions = c("Oxidative_phosphorylation" = "#F46D43",
                     "Cell_cycle" = "#708238",
                     "Immune_regulation" = "#9E0142",
                     "Signal_transduction" = "beige",
                     "Transcription" = "violet"),
  Group = c("Disease" = "darkgreen",
            "Control" = "blueviolet"),
  Lymphocyte_count = brewer.pal(5, 'PuBu')
)



# Base heatmap ===================================================
heat_plot <- pheatmap(data,
                      col = brewer.pal(10, 'RdYlGn'), # choose a colour scale for your data
                      cluster_rows = T, cluster_cols = T, # set to FALSE if you want to remove the dendograms
                      clustering_distance_cols = 'euclidean',
                      clustering_distance_rows = 'euclidean',
                      clustering_method = 'ward.D',
                      annotation_row = gene_functions_df, # row (gene) annotations
                      annotation_col = ann_df, # column (sample) annotations
                      annotation_colors = ann_colors, # colours for your annotations
                      annotation_names_row = F,
                      annotation_names_col = F,
                      fontsize_row = 10,          # row label font size
                      fontsize_col = 7,          # column label font size
                      angle_col = 45, # sample names at an angle
                      legend_breaks = c(-2, 0, 2), # legend customisation
                      legend_labels = c("Low", "Medium", "High"), # legend customisation
                      show_colnames = T, show_rownames = F, # displaying column and row names
                      main = "Super heatmap with annotations") # a title for our heatmap
