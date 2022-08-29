#20220218 nanaliu cytofworkflow from biotree
library(cytofWorkflow)
library(HDCytoData)
library(CATALYST)
library(readxl)
library(ggplot2)
library(reshape2)


fs <- Bodenmiller_BCR_XL_flowSet()
fs
ehub <- ExperimentHub()
ehub <- query(ehub,'HDCytoData')

url <- "http://imlspenticton.uzh.ch/robinson_lab/cytofWorkflow" 
# 这个cytof的panel的抗体信息表格：
panel <- "PBMC8_panel_v3.xlsx"
download.file(file.path(url, panel), destfile = panel, mode = "wb")
panel <- read_excel(panel)
head(data.frame(panel))



url <- "http://imlspenticton.uzh.ch/robinson_lab/cytofWorkflow"
md <- "PBMC8_metadata.xlsx"
download.file(file.path(url, md), destfile = md, mode = "wb")
md <- read_excel(md)
head(data.frame(md))
table(md[,3:4])

# spot check that all panel columns are in the flowSet object
all(panel$fcs_colname %in% colnames(fs))
# 有了样本的表型信息，panel的抗体信息，以及表达量矩阵，就可以构建对象：
# specify levels for conditions & sample IDs to assure desired ordering
md$condition <- factor(md$condition,
                       levels = c('Ref','BCRXL'))
md$sample_id <- factor(md$sample_id,
                       levels = md$sample_id[order(md$condition)])

sce <- prepData(fs,panel,md,
                features = panel$fcs_colname)
p <- plotExprs(sce,color_by = 'condition')
p$facet$params$ncol <- 6
p


if(!file.exists(file = 'input_for_cytofWorkflow.Rdata')){
  library(readxl)
  url <- "http://imlspenticton.uzh.ch/robinson_lab/cytofWorkflow"
  md <- "PBMC8_metadata.xlsx"
  download.file(file.path(url, md), destfile = md, mode = "wb")
  md <- read_excel(md)
  head(data.frame(md))
  table(md[,3:4])}
  # 样本的表型信息

dim(exprs(fs[1]))
exprs(fs[1])[1:6,1:5]
# unable to find an inherited method for function ‘exprs’ 
#for signature ‘"flowSet"’


panel
md
sce <- prepData(fs, panel, md, features = panel$fcs_colname)

pro='basic_qc'

p <- plotExprs(sce, color_by = "condition")
p$facet$params$ncol <- 4
p

n_cells(sce)
plotCounts(sce,group_by = 'sample_id',
           color_by = 'condition')
ggsave2(filename = paste0(pro,'_plotCounts.pdf'))
#save where, -getwd(). and the saved name and format

pbMDS(sce,color_by = 'condition',label_by = 'sample_id')

plotExprHeatmap(sce, scale = "last",
                hm_pal = rev(hcl.colors(10, "YlGnBu")))

plotNRS(sce, features = "type", color_by = "condition")



set.seed(1234)
sce <- cluster(sce,features = 'type',
               xdim = 10,ydim = 10,
               maxK = 20,seed = 1234)
plotExprHeatmap(sce,features = 'type',
                by='cluster_id',k='meta20',
                row_clust = F,
                bars = T,perc = T)

plotClusterExprs(sce,k='meta20',features = 'type')

plotMultiHeatmap(sce,
                 hm1 = 'type',k='meta20',
                 row_anno = F,bars = T,perc = T)

sce <- runDR(sce, "TSNE", cells = 1e3, features = "type")

p1 <- plotDR(sce, "TSNE", color_by = "meta20") + 
  theme(legend.position = "none")
p1

plotDR(sce, "TSNE", color_by = "meta20", 
       facet_by = "sample_id")
plotDR(sce, "TSNE", color_by = "meta20",
       facet_by = "condition")

sce <- runDR(sce, "UMAP", cells = 1e3, features = "type")

p2 <- plotDR(sce, "UMAP", color_by = "meta20")
lgd <- get_legend(p2)
p2 <- p2 + theme(legend.position = "none")
plot_grid(p1, p2, lgd, 
          nrow = 1, rel_widths = c(5, 5, 2))


plotDR(sce, "UMAP", color_by = "meta20",
       facet_by = "sample_id")

plotDR(sce, "UMAP", color_by = "meta20",
       facet_by = "condition")

plotCodes(sce, k = "meta20")

plotMultiHeatmap(sce, 
                 hm1 = "type",  
                 k = "som100", 
                 m = "meta20", 
                 row_anno = FALSE, 
                 col_anno = FALSE, 
                 bars = TRUE, 
                 perc = TRUE)
sce@metadata

plot_grid(labels = c("A", "B"),
          plotDR(sce, "UMAP", color_by = "meta20"),
          plotDR(sce, "UMAP", color_by = "meta8"))

plotAbundances(sce, k = "meta20", by = "sample_id")
plotAbundances(sce, k = "meta20", by = "cluster_id", 
               shape_by = "patient_id")


p=plotExprHeatmap(sce, features = "type", 
                  by = "cluster_id", 
                  k = "meta20", 
                  row_clust = F,
                  bars = TRUE, 
                  perc = TRUE)

p

plotAbundances(sce, k = "meta20", 
               by = "sample_id")

plotAbundances(sce, k = "meta20", 
               by = "cluster_id", 
               shape_by = "patient_id")




ns <- table(cluster_id = cluster_ids(sce, 'meta20'), 
            sample_id = sample_ids(sce))
fq <- prop.table(ns, 2) * 100
df <- as.data.frame(fq)

head(df)
dat=dcast(df,cluster_id~sample_id)
ei <- metadata(sce)$experiment_info
ei
dat$p=apply(dat,1,function(x){
  t.test(as.numeric(x[-1])~ei$condition)$p.value
})
dat$p

(da_formula1 <- createFormula(ei, 
                              cols_fixed = "condition", 
                              cols_random = "sample_id"))
(da_formula2 <- createFormula(ei, 
                              cols_fixed = "condition", 
                              cols_random = c("sample_id", "patient_id")))

contrast <- createContrast(c(0, 1))
da_res1 <- diffcyt(sce, 
                   formula = da_formula1, 
                   contrast = contrast,
                   analysis_type = "DA", 
                   method_DA = "diffcyt-DA-GLMM",
                   clustering_to_use = "meta20", 
                   verbose = FALSE)
da_res2 <- diffcyt(sce, 
                   formula = da_formula2, 
                   contrast = contrast,
                   analysis_type = "DA", 
                   method_DA = "diffcyt-DA-GLMM",
                   clustering_to_use = "meta20", 
                   verbose = FALSE)

names(da_res1)
rowData(da_res1$res)
FDR_cutoff=0.05
table(rowData(da_res1$res)$p_adj < FDR_cutoff)
table(rowData(da_res2$res)$p_adj < FDR_cutoff)
topTable(da_res2, show_props = TRUE, 
         format_vals = TRUE, digits = 2)
plotDiffHeatmap(sce, rowData(da_res2$res), 
                all = TRUE, fdr = FDR_cutoff)

df=topTable(da_res2, 
            show_props = TRUE, 
            format_vals = TRUE, 
            digits = 2,
            show_all_cols=T)
df=as.data.frame(df)
df

