library(gprofiler2)
library(htmlwidgets)

# Iterate through up and down genes
for (i in 1:length(snakemake@input[])) {
  genelist <- read.csv(snakemake@input[[i]], header=FALSE)[,1]
  
  if(length(genelist)>0)
  {
    gostres <- gost(query = genelist ,organism = snakemake@config[["organism"]])
    if(!is.null(gostres))
    {
      p <- gostplot(gostres, interactive = TRUE)
      htmlwidgets::saveWidget(p,snakemake@output[[i]])
    }
    else
    {
      system(paste("echo No enrichment > ",snakemake@output[[i]],sep=''))
    }
  }
}
