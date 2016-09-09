library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicFeatures)
library(AnnotationDbi)
library(ggplot2)
library(tibble)
library(readr)

# Load or generate ENSEMBL txdb -------------------------------------------
if (!file.exists('ensembl.txdb')) {
  txdb <- makeTxDbFromBiomart(dataset="hsapiens_gene_ensembl")
  saveDb(txdb, 'ensembl.txdb')
} else {
  txdb <- loadDb('ensembl.txdb')
}


# Get transcript, gene and CDS ranges -------------------------------------
tr <- transcripts(txdb, columns=c('tx_id', 'tx_name', 'gene_id'))
seqlevelsStyle(tr) <- 'UCSC'
tr <- keepStandardChromosomes(tr)

gn <- genes(txdb)
seqlevelsStyle(gn) <- 'UCSC'
gn <- keepStandardChromosomes(gn)

cd <- cds(txdb, columns=c('cds_id', 'tx_name', 'gene_id'))
seqlevelsStyle(cd) <- 'UCSC'
cd <- keepStandardChromosomes(cd)
#TODO: Add 2kb upstream and 200bp downstream of genes


# Visualize length distr. differences -------------------------------------
x = aggregate(data.frame(med=width(tr), stringsAsFactors = F),
              list(gene=unlist(tr$gene_id)),
              median)
x$gene.length <- width(gn)[match(x$gene, gn$gene_id)]

ggplot(x, aes(gene.length, med)) +
  geom_point(color=densCols(x$gene.length, x$med),
             size=0.8) +
  theme_bw()

ggplot(x, aes(gene.length, med)) +
  geom_point(color=densCols(log10(x$gene.length), log10(x$med)),
             size=0.8) +
  scale_x_log10() +
  scale_y_log10() + theme_bw()

summary(width(gn))
summary(width(tr))
summary(width(cd))


# Get sequences -----------------------------------------------------------
#sq.dnastr <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gn)
#colSums(alphabetFrequency(sq.dnastr)) #ensure there is only ACGTN
sq.gn <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gn, as.character=T)
sq.gn <- chartr('ACGTN', '12345', sq.gn)
write_tsv(tibble(sequence=sq.gn,
                 gene_id=gn$gene_id), 'sequences_gene.tsv')

sq.tr <- getSeq(BSgenome.Hsapiens.UCSC.hg38, tr, as.character=T)
sq.tr <- chartr('ACGTN', '12345', sq.tr)
write_tsv(tibble(sequence=sq.tr,
                 tx_id=tr$tx_name,
                 gene_id=unlist(tr$gene_id)), 'sequences_transcript.tsv')


# sq.cd <- getSeq(BSgenome.Hsapiens.UCSC.hg38, cd, as.character=T)
# sq.cd <- chartr('ACGTN', '12345', sq.cd)
# write_tsv(tibble(sequence=sq.cd,
#                  chr=seqnames(cd),
#                  start=start(cd),
#                  end=end(cd),
#                  tx_id=unlist(cd$tx_name),
#                  gene_id=unlist(cd$gene_id)), 'sequences_cds.tsv')


gc()