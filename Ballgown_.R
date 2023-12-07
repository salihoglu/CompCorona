#Jason Walker, jason.walker[AT]wustl.edu
#Malachi Griffith, mgriffit[AT]wustl.edu
#Obi Griffith, obigriffith[AT]wustl.edu
#The Genome McDonnell Institute, Washington University School of Medicine

#R tutorial for Informatics for RNA-sequence Analysis workshops

#Tutorial_Part1

library(ballgown)
library(genefilter)
library(devtools)
library(ggplot2)
library(gplots)
library(GenomicRanges)

# Load phenotype data from a file we saved in the current working directory
dir()
pheno_data = read.csv("CONT_vs_MERS.csv") #CHANGE WITH YOUR FILE NAME !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Load ballgown data structure and save it to a variable "bg"
bg = ballgown(samples=as.vector(pheno_data$path), pData=pheno_data)

# Display a description of this object
bg

# Load all attributes including gene name
bg_table2 = texpr(bg, 'all')
bg_gene_names2 = unique(bg_table2[, 9:10])

# Save the ballgown object to a file for later use
save(bg, file='bg.rda')
load("bg.rda")
bg_table2 = texpr(bg, 'all')
# Perform differential expression (DE) analysis with no filtering
results_transcripts2 = stattest(bg, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")
results_genes2 = stattest(bg, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
results_genes2 = merge(results_genes2, bg_gene_names2, by.x=c("id"), by.y=c("gene_id"))
results_genes2[,"de"] = log2(results_genes2[,"fc"])

# Save a tab delimited file for both the transcript and gene results
write.table(results_transcripts2, "transcript_results2.tsv", sep="\t", quote=FALSE, row.names = FALSE)
write.table(results_genes2, "gene_results2.tsv", sep="\t", quote=FALSE, row.names = FALSE)

# Filter low-abundance genes. Here we remove all transcripts with a variance across the samples of less than one
bg_filt2 = subset (bg,"rowVars(texpr(bg)) > 1", genomesubset=TRUE)

# Load all attributes including gene name
bg_filt_table2 = texpr(bg_filt2 , 'all')
bg_filt_gene_names2 = unique(bg_filt_table2[, 9:10])

# Perform DE analysis now using the filtered data
results_filt_transcripts2 = stattest(bg_filt2, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")
results_filt_genes2 = stattest(bg_filt2, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
results_filt_genes2 = merge(results_filt_genes2, bg_filt_gene_names2, by.x=c("id"), by.y=c("gene_id"))
results_filt_genes2[,"de"] = log2(results_filt_genes2[,"fc"])

# Output the filtered list of genes and transcripts and save to tab delimited files
write.table(results_filt_transcripts2, "transcript_results_filtered2.tsv", sep="\t", quote=FALSE, row.names = FALSE)
write.table(results_filt_genes2, "gene_results_filtered2.tsv", sep="\t", quote=FALSE, row.names = FALSE)

# Identify the significant genes with p-value < 0.01
sig_transcripts2 = subset(results_filt_transcripts2, results_filt_transcripts2$pval<0.05)
sig_transcripts2 = sig_transcripts2[order(sig_transcripts2[,"pval"], decreasing=FALSE),]
sig_genes2 = subset(results_filt_genes2, results_filt_genes2$pval<0.05)
sig_genes2 = sig_genes2[order(sig_genes2[,"pval"], decreasing=FALSE),]
sig_genes2

# Output the signifant gene results to a pair of tab delimited files
write.table(sig_transcripts2, "transcript_results_sig2.tsv", sep="\t", quote=FALSE, row.names = FALSE)
write.table(sig_genes2, "gene_results_sig_06_03.tsv", sep="\t", quote=FALSE, row.names = FALSE)


#### Write a simple table of differentially expressed transcripts to an output file
#Each should be significant with a log2 fold-change >= 2
sigpi2 = which(results_filt_genes2[,"pval"]<0.05)
sigp2 = results_filt_genes2[sigpi2,]
write.table(sigp2, file="pVal05_genes_coming_short.txt", sep="\t", row.names=FALSE, quote=FALSE)
sigde2 = which(abs(sigp2[,"de"]) >= 0) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
sig_tn_de2 = sigp2[sigde2,]

#Order the output by or p-value and then break ties using fold-change
o = order(sig_tn_de2[,"qval"], -abs(sig_tn_de2[,"de"]), decreasing=FALSE)

output2 = sig_tn_de2[o,c("gene_name","id","fc","pval","qval","de")]
write.table(output2, file="fc-Sig_DE_genes_coming_short.txt", sep="\t", row.names=FALSE, quote=FALSE)

#View output
output

#You can open the file "SigDE.txt" in Excel, Calc, etc.
#It should have been written to the current working directory that you set at the beginning of the R tutorial
dir()

library(ComplexHeatmap)
library(RColorBrewer)
library(viridis)

annotation_col = data.frame(
  Class = factor(c("Mock","Mock","Mock","Mock","Mock","Mock","Mock","Mock","Mock","Mock","Mock","Mock","Mock","Mock","Mock","Mock",
                   "MERS","MERS","MERS","MERS","MERS","MERS","MERS","MERS","MERS","MERS","MERS","MERS","MERS","MERS","MERS","MERS"))
)

ann_colors = list(
  Class = c(MERS ="#E43B3B" , Mock ="#3B7FE4" ))
pheatmap(sig_genes2, treeheight_row = 0, annotation_col = annotation_col,annotation_colors = ann_colors, name = "value", color = viridis(10))



save(bg_filt, results_filt_transcripts, results_filt_genes, sig_transcripts, sig_genes, output, file='bg_filt_results.rda')

#To exit R type:
quit(save="yes")
