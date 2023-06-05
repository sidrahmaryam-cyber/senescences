library(Seurat)
library(ggplot2)
setwd("~/sidrah/senescence")
S = readRDS('S_prepared.h5ad.to.RDS')

# remove UNK 
Idents(S) = S@meta.data$cluster_name
s = subset(S, idents = c("19: Unk_2", "18: Unk_1"), invert= TRUE )

### cell aggregation 

s@meta.data$cluster_name <- gsub("^.{0,3}", "", s@meta.data$cluster_name)

s$aggregated = 'Proximal_Tubule'
s$aggregated[s$cluster_name== 'TAL_1'] = 'Loop_of_Henle'
s$aggregated[s$cluster_name== 'TAL_2'] = 'Loop_of_Henle'
s$aggregated[s$cluster_name== 'DCT_1'] = 'Distal_Tubule'
s$aggregated[s$cluster_name== ' DCT_2'] = 'Distal_Tubule'
s$aggregated[s$cluster_name== ' IC-A'] = 'Connecting_segment_&_collecting_duct'
s$aggregated[s$cluster_name== ' IC-B'] = 'Connecting_segment_&_collecting_duct'

s$aggregated[s$cluster_name== 'PC-CD/CNT'] = 'Connecting_segment_&_collecting_duct'
s$aggregated[s$cluster_name== 'vSMC/Fib'] = 'Vascular_and_Interstitial'
s$aggregated[s$cluster_name== 'EC'] = 'Endothelial_cells'
s$aggregated[s$cluster_name== ' Mac'] = 'Leukocytes'

s$aggregated[s$cluster_name== ' Mast-Cell'] = 'Leukocytes'
s$aggregated[s$cluster_name== ' T-Cell'] = 'Leukocytes'
s$aggregated[s$cluster_name== ' B-Cell'] = 'Leukocytes'
s$aggregated[s$cluster_name== ' Pod'] = 'Podocytes'

table(s$aggregated, s$cluster_name)
table(s$aggregated)

### plotting

### plots

pdf(file = "umap_treatment.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 7)
DimPlot(s, reduction = 'UMAP', group.by = 'treatment')+theme(text = element_text(size = 15, face = 'bold'))
dev.off()

pdf(file = "umap_cluster.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5)
DimPlot(s, reduction = 'UMAP', group.by = 'cluster_name')+theme(text = element_text(size = 15, face = 'bold'))
dev.off()

pdf(file = "umap_aggregated.pdf",   # The directory you want to save the file in
    width = 9, # The width of the plot in inches
    height = 6)
DimPlot(s, reduction = 'UMAP', group.by = 'aggregated')+theme(text = element_text(size = 15, face = 'bold'))
dev.off()

#################################### Dot plot

#### removing leukocytes

non_immun = subset(s, subset = aggregated== 'Leukocytes', invert= TRUE)

sen_gen = c('CDKN1A',
            'CDKN2A',
            'CDKN2B',
            'BCL2',
            'CCL2',
            'CXCL8',
            'IL1A',
            'IL1B',
            'IL6',
            'MIF',
            'NFKB1',
            'PDGFA',
            'TP53')

pdf(file = "dotplot2_NON_immune.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 4)
DotPlot(non_immun, features = sen_gen, group.by = 'treatment', cols = c("red", "blue"))+theme(text = element_text(size = 15, face = 'bold'))+ ggtitle('Senescence genes')
dev.off()

#### adpkd gene plot
adkpd =c('MYC','PKD1', 'PKD2')
pdf(file = "adpkd_genes_non_immune.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 4)
DotPlot(non_immun, features = adkpd, group.by = 'treatment', cols = c("red", "blue"))+theme(text = element_text(size = 15, face = 'bold'))+ ggtitle('ADPKD genes')
dev.off()



################################# feature plots


pdf(file = "CDKN1A_feature.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 4)
FeaturePlot(s, features = 'CDKN1A', split.by = 'treatment')
dev.off()


pdf(file = "CDKN2A_feature.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 4)
FeaturePlot(s, features = 'CDKN2A', split.by = 'treatment')
dev.off()


pdf(file = "CDKN2B_feature.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 4)
FeaturePlot(s, features = 'CDKN2B', split.by = 'treatment')
dev.off()

################################ Dotplot

c = c('CDKN1A'
      , 'CDKN2A' , 'CDKN2B', 'CCL2', 'CXCL8')

pdf(file = "DotPlot_condition_genes.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 10)
DotPlot(s, features = c, group.by = 'aggregated',split.by = 'treatment', cols = c("red", "blue"))+theme(text = element_text(size = 15, face = 'bold'))
dev.off()
 


#######
library(ggplot2)



heal = subset(s, subset = treatment == 'control')
cond = subset(s, subset= treatment == 'ADKPD')


c = 'CXCL8'
#, 'CDKN2A' , 'CDKN2B', 'CCL2', 'CXCL8')

p1 = DotPlot(heal, features = c,group.by = 'aggregated',) + NoLegend() + ggtitle('Control')+theme(text = element_text(size = 15, face = 'bold'))
p2 = DotPlot(cond, features = c, group.by = 'aggregated') + ggtitle('ADKPD')+ theme(axis.title.y = element_blank(),axis.ticks=element_blank(),
                                                           axis.text.y = element_blank())+theme(text = element_text(size = 15, face = 'bold'))




pdf(file = "CXCL8.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10)
p1 + p2
dev.off()
