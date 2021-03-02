library(monocle)
## load data from seurat object
names(pbmc)
result=pbmc[[]]
result=result[result$seurat_clusters %in% c(0,2),] #c("Naive CD4 T", "Memory CD4 T")
head(result)
result$seurat_clusters=gsub(0,"NaiveT",result$seurat_clusters)
result$seurat_clusters=gsub(2,"MemoryT",result$seurat_clusters)
head(result)
dim(result)
count_data=pbmc[['RNA']]@counts

## build data frame for monocle
sample_sheet=result
gene_annotation=data.frame(rownames(count_data))
rownames(gene_annotation)=rownames(count_data)
colnames(gene_annotation)='gene_short_name'
expr_matrix=count_data[,rownames(result)]
## Store Data in a CellDataSet Object
pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)

## monocle object
HSMM <- newCellDataSet(expr_matrix,
                       phenoData = pd,
                       featureData = fd,
                       expressionFamily=negbinomial.size())
#Size factors help us normalize for differences in mRNA recovered across cells, and 
#"dispersion" values will help us perform differential expression analysis later.
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

## gene selection
HSMM <- detectGenes(HSMM, min_expr = 0.1)
print(head(fData(HSMM)))

HSMM <- detectGenes(HSMM, min_expr = 0.1)
print(head(fData(HSMM)))
expressed_genes <- row.names(subset(fData(HSMM),
                                    num_cells_expressed >= 10))
print(head(pData(HSMM)))

## cell selection
# valid_cells <- row.names(subset(pData(HSMM),
#                                 Cells.in.Well == 1 &
#                                   Control == FALSE &
#                                   Clump == FALSE &
#                                   Debris == FALSE &
#                                   Mapped.Fragments > 1000000))
valid_cells <- row.names(subset(pData(HSMM)))
HSMM <- HSMM[,valid_cells]

# Log-transform each value in the expression matrix.
L <- log(exprs(HSMM[expressed_genes,]))

#install.packages("reshape")
library(reshape)
# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))

# Plot the distribution of the standardized gene expression values.
qplot(value, geom = "density", data = melted_dens_df) +
  stat_function(fun = dnorm, size = 0.5, color = 'red') +
  xlab("Standardized log(FPKM)") +
  ylab("Density")

##Clustering cells without marker genes
disp_table <- dispersionTable(HSMM)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
plot_ordering_genes(HSMM)

# HSMM@auxClusteringData[["tSNE"]]$variance_explained <- NULL
plot_pc_variance_explained(HSMM, return_all = F) # norm_method='log'

HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 6,
                        reduction_method = 'tSNE', verbose = T)
HSMM <- clusterCells(HSMM, num_clusters = 2)
plot_cell_clusters(HSMM, 1, 2, color = "CellType",
                   markers = c("IL7R", "GZMK"))

plot_cell_clusters(HSMM, 1, 2, color = "seurat_clusters")

##Clustering cells using markers

## please pratice by yourself

## Trajectory step 1: choose genes that define a cell's progress
## differential genes
## Time-consuming, ~5 minutes
diff_test_res <- differentialGeneTest(HSMM[expressed_genes,],
                                      fullModelFormulaStr = "~seurat_clusters")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))

HSMM_myo <- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM_myo)

## Trajectory step 2: reduce data dimensionality
HSMM_myo <- reduceDimension(HSMM_myo, max_components = 2,
                            method = 'DDRTree')

## Trajectory step 3: order cells along the trajectory
HSMM_myo <- orderCells(HSMM_myo)
tpc=plot_cell_trajectory(HSMM_myo, color_by = "seurat_clusters")
tps=plot_cell_trajectory(HSMM_myo, color_by = "State")
tpt=plot_cell_trajectory(HSMM_myo, color_by = "Pseudotime")
plot_cell_trajectory(HSMM_myo, color_by = "State") +
  facet_wrap(~State, nrow = 1)

HSMM_expressed_genes <-  row.names(subset(fData(HSMM_myo),
                                          num_cells_expressed >= 10))
HSMM_filtered <- HSMM_myo[HSMM_expressed_genes,]
head(ordering_genes)
my_genes <- row.names(subset(fData(HSMM_filtered),
                             gene_short_name %in% c("HES4",    "ISG15",   "TNFRSF4", "MIB2")))
cds_subset <- HSMM_filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "seurat_clusters")

## Analyzing Branches in Single-Cell Trajectories
## about ~5 mins
BEAM_res <- BEAM(HSMM_myo, branch_point = 1, cores = 15)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
plot_genes_branched_heatmap(HSMM_myo[row.names(subset(BEAM_res,
                                                  qval < 1e-3)),],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)

