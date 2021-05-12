library(biomaRt)

genelist = rownames(read.csv('../data/raw/190517_4x4_alignments/DMSO_1_BIOReadsPerGene.out.tab',
                   sep='\t', row.names = 1))
species='hsapiens_gene_ensembl'
from='ensembl_gene_id'
to='hgnc_symbol'

mart = useMart("ensembl", dataset = species)

genelist=getBM(values=genelist,attributes = c(from,to), 
               filters = from,mart = mart)
write.csv(genelist, '../results/hgnc_symbol.csv')
