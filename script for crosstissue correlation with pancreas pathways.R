setwd('G:/My Drive/lab files/pancreas targeting/HPR and PCD analyses')
load('G:/My Drive/lab files/sex-difference myokine study/github/raw files/GTEx NA included env.RData')
library(WGCNA)
library(RColorBrewer)

library(pheatmap)
library(MetBrewer)
library(reshape2)
library(ggplot2)
library(dplyr)
library(forcats)
library(enrichR)
working_dataset=GTEx_subfiltered
row.names(working_dataset) = working_dataset$gene_tissue
working_dataset$gene_tissue=NULL
working_dataset = as.data.frame(t(working_dataset))
test1 = working_dataset[,grepl('ITIH5', colnames(working_dataset))]
colnames(test1)


sex_table = read.delim('G:/My Drive/Datasets/Human/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt')
sex_table$GTEx_ID = gsub('GTEX-', '', sex_table$SUBJID)
sex_table$sexMF = ifelse(sex_table$SEX==1, 'M', 'F')
table(sex_table$sexMF)
new_trts = sex_table[sex_table$GTEx_ID %in% row.names(working_dataset),]
table(new_trts$sexMF)
#  F   M 
#327 653

males = new_trts[new_trts$sexMF=='M',]
females = new_trts[!new_trts$sexMF=='M',]

#Subset two datasets for expression based on sex
working_datasetM = working_dataset[row.names(working_dataset) %in% males$GTEx_ID,]
working_datasetF = working_dataset[row.names(working_dataset) %in% females$GTEx_ID,]

tissue1 <- working_dataset[,grepl('Adipose - Subcutaneous', colnames(working_dataset)) | grepl('Adipose - Visceral (Omentum)', colnames(working_dataset), fixed=T) | grepl('Brain - Hypothalamus', colnames(working_dataset), fixed=T) | grepl('Colon - Transverse', colnames(working_dataset), fixed=T) | grepl('Spleen', colnames(working_dataset), fixed=T) | grepl('Small Intestine - Terminal Ileum', colnames(working_dataset), fixed=T) | grepl('Artery - Coronary', colnames(working_dataset), fixed=T) | grepl('Stomach', colnames(working_dataset), fixed=T) | grepl('Thyroid', colnames(working_dataset), fixed=T) | grepl('Pancreas', colnames(working_dataset), fixed=T) | grepl('Muscle - Skeletal', colnames(working_dataset), fixed=T) | grepl('Pituitary', colnames(working_dataset), fixed=T) | grepl('Liver', colnames(working_dataset), fixed=T) | grepl('Kidney - Cortex', colnames(working_dataset), fixed=T) | grepl('Heart - Left Ventricle', colnames(working_dataset), fixed=T) | grepl('Colon - Sigmoid', colnames(working_dataset), fixed=T) | grepl('Adrenal Gland', colnames(working_dataset), fixed=T) |  grepl('Artery - Aorta', colnames(working_dataset), fixed=T),]

########################################################################################
#see where HPR is coming and correlation with insulin

deg_table = read.delim('G:/My Drive/lab files/pancreas targeting/GTEx islet enrichments/datasets/uniprot-human-genes and goterms mapping.tab')

path_set = deg_table[grepl('insulin secretion', deg_table$Gene.ontology..biological.process.),]
ins_genes = path_set$Gene.names...primary..
ins_genes = paste0(ins_genes, '_', 'Pancreas')

tissue_set = colnames(working_dataset[,grepl('ITIH5', colnames(working_dataset))])
tissue_set = gsub('ITIH5_', '', tissue_set, fixed = T)
hpr_set = paste0('HPR_', tissue_set)

origin = tissue1[,colnames(tissue1) %in% hpr_set]
target = tissue1[,colnames(tissue1) %in% ins_genes]


inter_cors = bicorAndPvalue(target, target, use = 'p')
new_heat1 = inter_cors$bicor
pvals = ifelse(inter_cors$p<0.01, '8', '')
pdf(file = 'Internal Cor structure - Insulin Secretion Genes.pdf')
breaksList = seq(-1, 1, by = .1)
pheatmap(new_heat1, fontsize_number = 3, display_numbers = pvals, show_rownames = F, show_colnames = F,  number_color = "black", color = colorRampPalette(rev(brewer.pal(n = 7, name = "PiYG")))(length(breaksList)), breaks = breaksList, main='Insulin Secretion Genes', fontsize_row = 2, fontsize_col = 5, cluster_rows = T, cluster_cols = T)
dev.off()



full_cors = bicorAndPvalue(origin, target, use = 'p')
cor_table = melt(full_cors$bicor)
head(cor_table)
colnames(cor_table) = c('HPR_location', 'gene_tissue', 'bicor')
new_p = melt(full_cors$p)

cor_table$pvalue = new_p$value
cor_table = na.omit(cor_table)
qest = qvalue(cor_table$pvalue)
cor_table$qvalue = qest$qvalues
cor_table$gene_symbol = gsub("\\_.*","",cor_table$gene_tissue)
cor_table$tissue = gsub(".*_","",cor_table$HPR_location)
cor_table = cor_table[!is.na(cor_table$tissue),]
cor_table = cor_table[order(cor_table$pvalue, decreasing = F),]
cc1 = cor_table
cc1$logp = -log10(cc1$pvalue)
head(cc1)
pdf(file = 'Tissue-specific HPR correlations with pancreas insulin secretion genes.pdf')
ggplot(cc1, aes(x=fct_reorder(tissue, logp), y=logp, fill=tissue)) + geom_violin()+ geom_boxplot(width=0.2) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + geom_hline(yintercept=median(cc1$logp), color='firebrick3') + xlab('') + ylab('-log10(pvalue) HPR-Insulin Secretion Gene Correlation')+ ggtitle('Tissue-specific HPR correlations with pancreas insulin secretion genes') + theme(legend.position = "none")
dev.off()






#################################################
#run HPR enrichments
gene_1 = 'HPR'
candid_gene = gene_1
origin_tissue = 'Liver'
gene_tissue1 = paste0(gene_1, '_', origin_tissue)
tissue1 <- working_dataset[,grepl('Adipose - Subcutaneous', colnames(working_dataset)) | grepl('Adipose - Visceral (Omentum)', colnames(working_dataset), fixed=T) | grepl('Brain - Hypothalamus', colnames(working_dataset), fixed=T) | grepl('Colon - Transverse', colnames(working_dataset), fixed=T) | grepl('Spleen', colnames(working_dataset), fixed=T) | grepl('Small Intestine - Terminal Ileum', colnames(working_dataset), fixed=T) | grepl('Artery - Coronary', colnames(working_dataset), fixed=T) | grepl('Stomach', colnames(working_dataset), fixed=T) | grepl('Thyroid', colnames(working_dataset), fixed=T) | grepl('Pancreas', colnames(working_dataset), fixed=T) | grepl('Muscle - Skeletal', colnames(working_dataset), fixed=T) | grepl('Pituitary', colnames(working_dataset), fixed=T) | grepl('Liver', colnames(working_dataset), fixed=T) | grepl('Kidney - Cortex', colnames(working_dataset), fixed=T) | grepl('Heart - Left Ventricle', colnames(working_dataset), fixed=T) | grepl('Colon - Sigmoid', colnames(working_dataset), fixed=T) | grepl('Adrenal Gland', colnames(working_dataset), fixed=T) |  grepl('Artery - Aorta', colnames(working_dataset), fixed=T),]
tissue1 = as.data.frame(tissue1)

#tissue1 = tissue1[!row.names(tissue1) %in% males$GTEx_ID,]
origin = tissue1[,colnames(tissue1)==gene_tissue1]
target = tissue1[,!colnames(tissue1)==gene_tissue1]

full_cors = bicorAndPvalue(origin, target, use = 'p')
cor_table = melt(full_cors$bicor)
cor_table$Var1=NULL
colnames(cor_table) = c('gene_tissue', 'bicor')
new_p = melt(full_cors$p)

cor_table$pvalue = new_p$value[match(cor_table$gene_tissue, new_p$Var2)]
cor_table = na.omit(cor_table)
qest = qvalue(cor_table$pvalue)
cor_table$qvalue = qest$qvalues
cor_table$gene_symbol = gsub("\\_.*","",cor_table$gene_tissue)
cor_table$tissue = gsub(".*_","",cor_table$gene_tissue)
cor_table = cor_table[!is.na(cor_table$tissue),]
cor_table = cor_table[order(cor_table$pvalue, decreasing = F),]

res1 = cor_table[cor_table$pvalue<0.01,]
res1 = na.omit(res1)
write.csv(res1, file = paste0('significant tissue enrichments ', gene_tissue1, '.csv'), row.names = F)

sig_table = cor_table[cor_table$qvalue<0.1,]
sig_table$qcat =ifelse(sig_table$qvalue<0.01, 'q<0.01', 'q<0.1')
sig_table$qcat =ifelse(sig_table$qvalue<0.0001, 'q<0.0001', paste0(sig_table$qcat))
table(sig_table$tissue[sig_table$qcat=='q<0.01'])

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=25, face="bold")
  )

binned_sig_prots= sig_table %>%
  dplyr::group_by(qcat, tissue) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))

tissue_freqs = binned_sig_prots %>% dplyr::group_by(qcat) %>% dplyr::summarise(sum(n))

binned_sig_prots$tot_count = tissue_freqs$`sum(n)`[match(binned_sig_prots$qcat, tissue_freqs$qcat)]
binned_sig_prots$qcat = paste0(binned_sig_prots$qcat, ' ', binned_sig_prots$tot_count, ' genes')
col_scheme = rev(met.brewer('Austria', length(unique(sig_table$tissue))))
names(col_scheme) = unique(sig_table$tissue)
pdf(file = paste0('crosstissue gene enrichments ', gene_tissue1, '.pdf'))
ggplot(binned_sig_prots, aes(x = "", y = freq, fill =tissue)) + 
  geom_bar(stat = "identity", width = 1, position = position_fill()) + blank_theme + theme(plot.title=element_text(size=25, face="bold")) +
  theme(axis.text.x=element_blank())+ scale_fill_manual(values=col_scheme) +
  coord_polar(theta = "y") + 
  facet_wrap( ~ qcat)
dev.off()

########perform when removing origin
sig_table1 = sig_table[!sig_table$tissue==origin_tissue,]
binned_sig_prots= sig_table1 %>%
  dplyr::group_by(qcat, tissue) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))

tissue_freqs = binned_sig_prots %>% dplyr::group_by(qcat) %>% dplyr::summarise(sum(n))

binned_sig_prots$tot_count = tissue_freqs$`sum(n)`[match(binned_sig_prots$qcat, tissue_freqs$qcat)]
binned_sig_prots$qcat = paste0(binned_sig_prots$qcat, ' ', binned_sig_prots$tot_count, ' genes')
col_scheme = rev(met.brewer('Austria', length(unique(sig_table$tissue))))
names(col_scheme) = unique(sig_table$tissue)
pdf(file = paste0('crosstissue (origin removed) gene enrichments ', gene_tissue1, '.pdf'))
ggplot(binned_sig_prots, aes(x = "", y = freq, fill =tissue)) + 
  geom_bar(stat = "identity", width = 1, position = position_fill()) + blank_theme + theme(plot.title=element_text(size=25, face="bold")) +
  theme(axis.text.x=element_blank())+ scale_fill_manual(values=col_scheme) +
  coord_polar(theta = "y") + 
  facet_wrap( ~ qcat)
dev.off()

binned_sig_prots= sig_table %>%
  dplyr::group_by(qcat, tissue) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))


tissue_list = binned_sig_prots[binned_sig_prots$qcat=='q<0.01',]
tissue_list = tissue_list[order(tissue_list$n, decreasing = T),]

select_tissue = tissue_list$tissue[1]
pp1 = res1[res1$qvalue<0.1,]
pp1 = pp1[pp1$tissue %in% select_tissue,]
pp1_length = ifelse(length(row.names(pp1)) > 500, as.numeric(500), as.numeric(length(row.names(pp1))))
pp2 = pp1[1:pp1_length,]
gg1 = pp2$gene_symbol

setEnrichrSite("Enrichr")
dbs <- listEnrichrDbs()
dbs1 <- c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "Reactome_2022", "DSigDB")

enriched <- enrichr(gg1, dbs1)
names(enriched[1])
pdf(file = paste0(select_tissue, ' gene correlations with ', candid_gene, ' ', names(enriched[1]), '.pdf'))
plotEnrich(enriched[[1]], showTerms = 15, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()
pdf(file = paste0(select_tissue, ' gene correlations with ', candid_gene, ' ', names(enriched[2]), '.pdf'))
plotEnrich(enriched[[2]], showTerms = 15, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()
pdf(file = paste0(select_tissue, ' gene correlations with ', candid_gene, ' ', names(enriched[3]), '.pdf'))
plotEnrich(enriched[[3]], showTerms = 15, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()

select_tissue = tissue_list$tissue[2]
pp1 = res1[res1$qvalue<0.1,]
pp1 = pp1[pp1$tissue %in% select_tissue,]
pp1_length = ifelse(length(row.names(pp1)) > 500, as.numeric(500), as.numeric(length(row.names(pp1))))
pp2 = pp1[1:pp1_length,]
gg1 = pp2$gene_symbol

setEnrichrSite("Enrichr")
dbs <- listEnrichrDbs()
dbs1 <- c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "Reactome_2022", "DSigDB")

enriched <- enrichr(gg1, dbs1)
names(enriched[1])
pdf(file = paste0(select_tissue, ' gene correlations with ', candid_gene, ' ', names(enriched[1]), '.pdf'))
plotEnrich(enriched[[1]], showTerms = 15, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()
pdf(file = paste0(select_tissue, ' gene correlations with ', candid_gene, ' ', names(enriched[2]), '.pdf'))
plotEnrich(enriched[[2]], showTerms = 15, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()
pdf(file = paste0(select_tissue, ' gene correlations with ', candid_gene, ' ', names(enriched[3]), '.pdf'))
plotEnrich(enriched[[3]], showTerms = 15, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()


select_tissue = tissue_list$tissue[3]
pp1 = res1[res1$qvalue<0.1,]
pp1 = pp1[pp1$tissue %in% select_tissue,]
pp1_length = ifelse(length(row.names(pp1)) > 500, as.numeric(500), as.numeric(length(row.names(pp1))))
pp2 = pp1[1:pp1_length,]
gg1 = pp2$gene_symbol

setEnrichrSite("Enrichr")
dbs <- listEnrichrDbs()
dbs1 <- c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "Reactome_2022", "DSigDB")

enriched <- enrichr(gg1, dbs1)
names(enriched[1])
pdf(file = paste0(select_tissue, ' gene correlations with ', candid_gene, ' ', names(enriched[1]), '.pdf'))
plotEnrich(enriched[[1]], showTerms = 15, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()
pdf(file = paste0(select_tissue, ' gene correlations with ', candid_gene, ' ', names(enriched[2]), '.pdf'))
plotEnrich(enriched[[2]], showTerms = 15, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()
pdf(file = paste0(select_tissue, ' gene correlations with ', candid_gene, ' ', names(enriched[3]), '.pdf'))
plotEnrich(enriched[[3]], showTerms = 15, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()

select_tissue = 'Pancreas'
pp1 = res1[res1$qvalue<0.1,]
pp1 = pp1[pp1$tissue %in% select_tissue,]
pp1_length = ifelse(length(row.names(pp1)) > 500, as.numeric(500), as.numeric(length(row.names(pp1))))
pp2 = pp1[1:pp1_length,]
gg1 = pp2$gene_symbol

setEnrichrSite("Enrichr")
dbs <- listEnrichrDbs()
dbs1 <- c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "Reactome_2022", "DSigDB")

enriched <- enrichr(gg1, dbs1)
names(enriched[1])
pdf(file = paste0(select_tissue, ' gene correlations with ', candid_gene, ' ', names(enriched[1]), '.pdf'))
plotEnrich(enriched[[1]], showTerms = 15, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()
pdf(file = paste0(select_tissue, ' gene correlations with ', candid_gene, ' ', names(enriched[2]), '.pdf'))
plotEnrich(enriched[[2]], showTerms = 15, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()
pdf(file = paste0(select_tissue, ' gene correlations with ', candid_gene, ' ', names(enriched[3]), '.pdf'))
plotEnrich(enriched[[3]], showTerms = 15, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()

