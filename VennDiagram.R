## Salihoglu ##
###################### v e n n D i a g r a m ###################################
library(ggplot2)
library(ggVennDiagram)

MERS <-read.table("/MERS_GSE56192/ballgown/ref_only/usingData/gene_results_sig.tsv")
MERS = MERS[-1,]
colnames(MERS) <- c("id", "feature", "fc","pval","qval","gene_name","de")
SARS <-read.table("/SARS_GSE56192/ballgown/ref_only/using_data/gene_results_sig.tsv")
SARS = SARS[-1,]
colnames(SARS) <- c("id", "feature", "fc","pval","qval","gene_name","de")
SARSCOV2 <-read.table("/SARS-COV2_GSE147507/GSE147507-gene_results_sig.tsv")
SARSCOV2 = SARSCOV2[-1,]
colnames(SARSCOV2) <- c("id", "feature", "fc","pval","qval","gene_name","de")
PBMC <-read.table("/PBMC_PRJCA002326/rs_de/de/ballgown/using_data/gene_results_sig.tsv")
PBMC = PBMC[-1,]
colnames(PBMC) <- c("id", "feature", "fc","pval","qval","gene_name","de")

x <- list( MERS_CoV = MERS$gene_name,
           SARS_CoV = SARS$gene_name,
           Lung_SARS_CoV2 = SARSCOV2$gene_name,
           Blood_SARS_CoV2 = PBMC$gene_name
)
venn <- ggVennDiagram(x,label_geom = "text",label_color = "Black", label_size = 4, edge_size = 1, label = "count")+ scale_color_brewer(palette = "Set1")+scale_fill_distiller(palette = "Blues", direction = 1)