### Codes for performing similarity analyses between pairwise cell groups 
# Here is an example showing the similarity between our LMC subclusters from the brachial and lumbar segments 

library(batchelor)
# load the integrated data including all LMC neurons from the brachial and lumbar segments 
load("e135_brachial_lumbar_merged_LMC_withSubclusters.RData")
# Get the cell labels given by seperate analysis of each segment
cell.levels.brachial <- c('cl1','cl2','cl3','cl4','rl1','rl2','rl3','rl4','rl5','cm1','cm2','cm3','rm1','rm2','rm3','rm4')
cell.levels.lumbar <- c("c/rl1", "c/rl2", "c/rl3", "c/rl4", "c/rl5", "c/rl6","cm1",   "cm2" ,  "cm3","rm" )
# Get the joint UMAP space
dr.umap = Seurat::combined@reductions$umap@cell.embeddings
dr.brachial = dr.umap[combined$conditions == "Brachial", ]
dr.lumbar = dr.umap[combined$conditions == "Lumbar", ]
# Identifying MNN and computing average overlap ratio using different numbers of neighbors 
k.mnn.all = c(10, 15, 20)
sim.avg <- 0
for (k.mnn in k.mnn.all) {
  out <- findMutualNN(dr.brachial, dr.lumbar, k1=k.mnn)
  mn = data.frame(first = out$first, second = out$second)
  a = combined$subclusters[combined$conditions == "Brachial"]
  a <- factor(a, levels = cell.levels.brachial)
  names(a) = 1:length(a)
  mn$first.name = a[mn$first]
  b = combined$subclusters[combined$conditions == "Lumbar"]
  b <- factor(b, levels = cell.levels.lumbar)
  names(b) = 1:length(b)
  mn$second.name = b[mn$second]
  
  mn$first.name <- factor(mn$first.name, levels = cell.levels.brachial)
  mn$second.name <- factor(mn$second.name, levels = cell.levels.lumbar)
  sim <- c()
  for (i in cell.levels.brachial) {
    sim <- rbind(sim, as.numeric(table(mn$second.name[mn$first.name == i])))
  }
  rownames(sim) <- cell.levels.brachial
  colnames(sim) <- cell.levels.lumbar
  # Computing an average overlap ratio 
  sim.avg <- sim.avg + sim/as.numeric(table(a))/2/k.mnn
}
sim.avg <- sim.avg/3

## Visualize the similarity between cell groups using heatmap
sim.scale = sim.avg
library(ComplexHeatmap)
library(colorspace)
color.heatmap = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = "RdYlBu")))(100)
ht = Heatmap(sim.scale,col =color.heatmap, cluster_rows = FALSE, cluster_columns = FALSE,show_heatmap_legend = T, row_names_side = "left",row_title = "Brachial",column_title = "Similarity", name = NULL,
              cell_fun = function(j, i, x, y, width, height, fill) {
                if(sim.scale[i, j] > 0.2)
                  grid.text(sprintf("%.2f", sim.scale[i, j]), x, y, gp = gpar(fontsize = 8))
              })
draw(ht)





