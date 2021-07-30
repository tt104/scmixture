library(plotly)
library(htmlwidgets)

# Open as cluster means as data.frame
cluster_data <- read.csv(snakemake@input[[1]], header=TRUE)

# Find the mean gene expression for each gene across all clusters
mean_values <- colMeans(cluster_data)
# Get the gene names from the columns
gene_names <- colnames(cluster_data)
# Create list of column numbers representing genes
gene_ids <- c(1:ncol(cluster_data))

create_plot <- function(cluster_num) {
  # Get the gene counts for the selected cluster
  cluster_values <- unname(cluster_data[cluster_num,])
  # Plot the mean gene expression across clusters
  fig <- plot_ly(cluster_data, x = ~gene_names, y = ~mean_values, type = 'bar', 
                 name = 'Mean Expression of All Clusters',
                 hoverlabel = list(namelength = -1),
                 marker = list(color = 'grey'))
  # Plot the gene expression of the selected cluster
  fig <- fig %>% add_markers(y = ~cluster_values, 
                             name = paste('Cluster', cluster_num, 'Mean Expression'),
                             marker = list(
                               color = ~cluster_values,
                               colorscale = 'Jet',
                               line = list(color = 'black', width = 1)))
  # Name the graph and hide the x-axis
  fig <- fig %>% layout(title = paste("Average Gene Expression of All Clusters Against Cluster", cluster_num),
                        xaxis = list(title = "Gene", zeroline = FALSE, showline = FALSE, 
                                     showticklabels = FALSE, showgrid = FALSE),
                        yaxis = list(title = "Gene Expression"),
                        plot_bgcolor='rgb(245, 245, 245)',
                        hovermode = "x unified")
  # Save the graph to HTML
  htmlwidgets::saveWidget(fig, file.path(gsub(" ", "", paste(snakemake@output[[1]], "/Cluster", cluster_num, "Means.html"))))
}

# Create the output directory if it does not already exist
dir.create(file.path(snakemake@output[[1]]), showWarnings = FALSE)

# Make a graph for each cluster
print(paste('Creating mean expression graphs for', nrow(cluster_data), 'clusters'))
for (cluster in 1:nrow(cluster_data)) {
  create_plot(cluster)
}
