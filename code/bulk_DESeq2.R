setwd("~/Documents/Projects/Combi-Seq-analysis/code")

library('BiocParallel')
library('DESeq2')

data = read.csv('../results/bulk/count_ploy.csv',
                sep=',', header=TRUE, row.names=1)
meta = read.csv('../results/bulk/meta_poly.csv',
                sep=',', header=TRUE, row.names=1)
meta$Drug = relevel(as.factor(meta$Drug), 'DMSO_DMSO')

dds=DESeqDataSetFromMatrix(countData = data,
                           colData = meta,
                           design = ~ Drug)
vsd <- vst(dds, blind=TRUE)
write.csv(assay(vsd), '../results/bulk/vst_poly.csv')
dds=DESeq(dds,parallel = T)
write.csv(results(dds, contrast = c('Drug', 'medium_medium', 'DMSO_DMSO')), '../results/bulk/DE-analysis/poly/medium_medium.csv')
write.csv(results(dds, contrast = c('Drug', 'Trametinib_DMSO', 'DMSO_DMSO')), '../results/bulk/DE-analysis/poly/Trametinib_DMSO.csv')
write.csv(results(dds, contrast = c('Drug', 'YM155_DMSO', 'DMSO_DMSO')), '../results/bulk/DE-analysis/poly/YM155_DMSO.csv')
write.csv(results(dds, contrast = c('Drug', 'Trametinib_YM155', 'DMSO_DMSO')), '../results/bulk/DE-analysis/poly/Trametinib_YM155.csv')

data = read.csv('../results/bulk/count_combi.csv',
                sep=',', header=TRUE, row.names=1)
meta = read.csv('../results/bulk/meta_combi.csv',
                sep=',', header=TRUE, row.names=1)
meta$Drug = relevel(as.factor(meta$Drug), 'DMSO_DMSO')

dds=DESeqDataSetFromMatrix(countData = data,
                           colData = meta,
                           design = ~ Drug)
vsd <- vst(dds, blind=TRUE)
write.csv(assay(vsd), '../results/bulk/vst_combi.csv')
dds=DESeq(dds,parallel = T)

write.csv(results(dds, contrast = c('Drug', 'Trametinib_DMSO', 'DMSO_DMSO')), '../results/bulk/DE-analysis/combi/Trametinib_DMSO.csv')
write.csv(results(dds, contrast = c('Drug', 'YM155_DMSO', 'DMSO_DMSO')), '../results/bulk/DE-analysis/combi/YM155_DMSO.csv')
write.csv(results(dds, contrast = c('Drug', 'Trametinib_YM155', 'DMSO_DMSO')), '../results/bulk/DE-analysis/combi/Trametinib_YM155.csv')
write.csv(results(dds, contrast = c('Drug', 'YM155_Trametinib', 'DMSO_DMSO')), '../results/bulk/DE-analysis/combi/YM155_Trametinib.csv')