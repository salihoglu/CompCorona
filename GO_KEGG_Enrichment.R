## Salihoglu ##


library(clusterProfiler)
library(dplyr)
library(org.Hs.eg.db)
library(ggridges)
library(enrichplot)
library(ComplexUpset)
library(ComplexHeatmap)
library(ggplot2)
library(UpSetR)

###################### SARSCOV2 ########################################
SARSCOV2 <-read.table("/SARS-COV2_GSE147507/GSE147507-gene_results_sig.tsv")
SARSCOV2 = SARSCOV2[-1,]
colnames(SARSCOV2) <- c("id", "feature", "fc","pval","qval","gene_name","de")

SARSCOV2$fc <- as.numeric(SARSCOV2$fc)

SARSCOV2$diffexpressed <- "Not sig"
SARSCOV2$diffexpressed[SARSCOV2$fc > 1] <- "Upregulated"
SARSCOV2$diffexpressed[SARSCOV2$fc < 1] <- "Downregulated"
genes_entrez<- bitr(SARSCOV2$gene_name, fromType="SYMBOL",
                    toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
colnames(genes_entrez) <- c("gene_name","ENTREZID")
update_SARSCOV2 <- merge(genes_entrez, SARSCOV2, by="gene_name")

original_gene_list <- update_SARSCOV2$fc
names(original_gene_list) <- update_SARSCOV2$ENTREZID
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

g <- scale_colour_gradient(
  low = "#00BB8E",
  high = "#6F00BB",
  space = "Lab",
  guide = "colourbar",
  aesthetics = "colour"
)
s <- labs(colour = "logFC") 

#________________________________________________________________ GO enrichment

go_enr<- enrichGO(gene         = names(gene_list),
                  OrgDb         = org.Hs.eg.db,
                  keyType       = "ENTREZID",
                  ont           = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  minGSSize = 1,
                  maxGSSize = 1500,
                  readable      = TRUE)
go_enr_df <- as.data.frame(go_enr)

go2 <- pairwise_termsim(go_enr)

gse <- gseGO(geneList= gene_list, 
             ont ="ALL", 
             keyType = "ENTREZID", 
             pvalueCutoff = 0.05, 
             OrgDb = 'org.Hs.eg.db',
             scoreType = "pos",
             verbose      = FALSE
             
)

gse_df <- as.data.frame(gse)
#________________________________________________________________ KEGG enrichment

kegg_enr<- enrichKEGG(gene = update_SARSCOV2$ENTREZID, 
                      organism= "hsa",
                      pAdjustMethod = "BH",
                      minGSSize = 1,
                      maxGSSize = 600,
                      pvalueCutoff = 0.05
)
kegg_enr <- setReadable(kegg_enr, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
kegg_enr_df <- data.frame(kegg_enr)


gse_kegg <- gseKEGG(geneList=gene_list, 
                    organism= "hsa",
                    pAdjustMethod = "none",
                    minGSSize = 1,
                    maxGSSize = 600,
                    keyType       = "ncbi-geneid",
                    pvalueCutoff = 0.05
)
gse_kegg_df <- as.data.frame(gse_kegg)

#__________________________________________________________________PLOTS

ridgeplot(gse_kegg, showCategory = 15) + labs(x = "SARS-CoV2 Enrichment distribution")+geom_density_ridges_gradient() +
  scale_fill_viridis_c()

treeplot(go2, showCategory = 50,color = "p.adjust")

emapplot(go2, showCategory = 20)+ viridis::scale_fill_viridis("pAdjust")+ggtitle("GO Enrichment for SARS-CoV2")

upsetplot(go_enr, 10)

dotplot(gse, showCategory=45, split=".sign",font.size = 6) + facet_grid(.~.sign)+g+ggtitle("GO Enrichment for SARS-COV2")

dotplot(gse_kegg, showCategory=45, split=".sign",font.size = 6) + facet_grid(.~.sign)+g+ggtitle("KEGG Enrichment for SARS-COV2")

barplot(kegg_enr, showCategory = 20)+ viridis::scale_fill_viridis()

cnetplot(kegg_enr, showCategory = 15, foldChange =gene_list,  colorEdge = TRUE)+ g  + s

cnetplot(go_enr, categorySize="pvalue",
         colorEdge = FALSE,
         cex_category = 2,
         showCategory = 10, 
         cex_gene = 3 ,
         layout = "kk",
         fixed = TRUE,
         foldChange = gene_list, 
         vertex.label.font=4) + g  + s

heatplot(go_enr, foldChange=gene_list, showCategory=20)
#___________________________________________________________________________________
SARS_gene <- SARS$gene_name
MERS_gene <- MERS$gene_name
HCOV_gene <- HCOV$gene_name
Blood_SARSCOV2_gene <- PBMC$gene_name
Lung_SARSCOV2_gene <- SARSCOV2$gene_name
s <- combine(SARS_gene, MERS_gene, Blood_SARSCOV2_gene, Lung_SARSCOV2_gene)

lt <- list( Sars_CoV = SARS_gene,
            Mers_CoV = MERS_gene,
            Blood_Sars_CoV2 = Blood_SARSCOV2_gene,
            Lung_Sars_CoV2 = Lung_SARSCOV2_gene
)

testr <- as.data.frame(lt)
tests <- as.data.frame(SARSCOV2$fc)
colnames(tests) <- "SARS_CoV2_FC"
rownames(tests) <- SARSCOV2$gene_name
testrs <- merge(testr, tests, by=0)

Infections = c('Sars_CoV', 'Mers_CoV', 'HcoV_229E', 'Blood_Sars_CoV2', 'Lung_Sars_CoV2')      
lt <- list_to_matrix(lt)
merge(SARSCOV2$fc, lt, by=0)
m1 = make_comb_mat(lt) # the default mode is `distinct`
m2 = make_comb_mat(lt, mode = "intersect")
m3 = make_comb_mat(lt, mode = "union")
m = make_comb_mat(lt)
set_name(m)
comb_name(m)
set_size(m)
comb_size(m)
comb_degree(m)
t(m)
m = make_comb_mat(lt)
UpSet(m)
UpSet(m,  comb_order = order(comb_size(m)))
UpSet(m2,comb_order = order(comb_size(m2)))

################################################
ss = set_size(m2)
cs = comb_size(m2)
ht = UpSet(m2, 
           set_order = order(ss),
           comb_order = order(comb_degree(m2), -cs),
           top_annotation = HeatmapAnnotation(
             "Gene Intersections" = anno_barplot(cs, 
                                                 ylim = c(0, max(cs)*1.1),
                                                 border = FALSE, 
                                                 gp = gpar(fill = "black"), 
                                                 height = unit(9, "cm")
             ), 
             annotation_name_side = "left", 
             annotation_name_rot = 90),
           left_annotation = rowAnnotation(
             "Gene Size" = anno_barplot(-ss, 
                                        baseline = 0,
                                        axis_param = list(
                                          at = c(0, -1000,  -2000, -3000,-4000),
                                          labels = c(0,  1000,  2000, 3000,4000),
                                          labels_rot = 0),
                                        border = FALSE, 
                                        gp = gpar(fill = "black"), 
                                        width = unit(4, "cm")
             ),
             set_name = anno_text(set_name(m2), 
                                  location = 0.5, 
                                  just = "center",
                                  width = max_text_width(set_name(m2)) + unit(4, "mm"))
           ), 
           right_annotation = NULL,
           show_row_names = FALSE)
ht = draw(ht)
od = column_order(ht)
decorate_annotation("Gene Intersections", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = c("left", "bottom"), 
            gp = gpar(fontsize = 6, col = "#404040"), rot = 45)
})
###########################################################################
set_size(8, 8)

set.seed(0)   # keep the same jitter for identical plots

upset(
  lt,
  Infections,
  annotations = list(
    # 1st method - passing list:
    'Length'=list(
      aes=aes(x=intersection, y=length),
      # provide a list if you wish to add several geoms
      geom=geom_boxplot(na.rm=TRUE)
    ),
    # 2nd method - using ggplot
    'Rating'=(
      # note that aes(x=intersection) is supplied by default and can be skipped
      ggplot(mapping=aes(y=rating))
      # checkout ggbeeswarm::geom_quasirandom for better results!
      + geom_jitter(aes(color=log10(votes)), na.rm=TRUE)
      + geom_violin(alpha=0.5, na.rm=TRUE)
    ),
    # 3rd method - using `upset_annotate` shorthand
    'Budget'=upset_annotate('budget', geom_boxplot(na.rm=TRUE))
  ),
  min_size=10,
  width_ratio=0.1
)
############################################################################

colnames(SARS_gene) <- "genes"
colnames(MERS_gene) <- "genes"
colnames(Blood_SARSCOV2_gene) <- "genes"
colnames(Lung_SARSCOV2_gene) <- "genes"
set1 <- rbind(SARS_gene, MERS_gene, HCOV_gene, Blood_SARSCOV2_gene, Lung_SARSCOV2_gene)
set1 <- distinct(set1)
set1 <- as.data.frame(set1)
set2 <- combine(set1,s)
a <- as.list(s)


m2 = make_comb_mat(s, mode = "intersect")
upset(
  m2,
  base_annotations=list(
    'Intersection size'=intersection_size(counts=FALSE)
  ),
  min_size=100,
  width_ratio=0.1,
  height_ratio=0.7,
  sort_sets='ascending',
  themes=upset_modify_themes(
    list(
      'intersections_matrix'=theme(text=element_text(size=8)),
      'overall_sizes'=theme(axis.text.x=element_text(angle=90))
    )
  ),
  annotations=list(
    'Rating'=list(
      aes=aes(x=intersection, y=rating),
      geom=list(
        geom_violin(width=1.1, alpha=0.5),
        stat_summary(
          fun.y=mean,
          fun.ymin=function(x) { quantile(x, 0.25) },
          fun.ymax=function(x) { quantile(x, 0.75) },
          colour='grey'
        )
      )
    )
  ),
  queries=list(
    upset_query(
      intersect=c('MERS_gene', 'SARS_gene'),
      color='red',
      fill='red',
      only_components=c('intersections_matrix', 'Intersection size')
    ),
    upset_query(
      set='Drama',
      fill='blue'
    ),
    upset_query(
      intersect=c('SARS_gene', 'HCOV_gene'),
      fill='yellow',
      only_components=c('Rating')
    )
  )
)

##########################################################################
write.csv(gse_df,"GO_KEGG_Enrichment/SARS-COV2_GSE147507/gseGO_df_SARScov2.csv")
write.csv(gse_kegg_df,"GO_KEGG_Enrichment/SARS-COV2_GSE147507/gseKEGG_df_SARScov2.csv")
write.csv(kegg_enr_df,"GO_KEGG_Enrichment/SARS-COV2_GSE147507/kegg_enr_df_SARScov2.csv")
write.csv(go_enr_df,"GO_KEGG_Enrichment/SARS-COV2_GSE147507/go_enr_df_SARScov2.csv")

############################################################################