library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicFeatures)
library(AnnotationDbi)
library(ggplot2)
library(tibble)
library(readr)
library(biomaRt)

set.seed(42)

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


# Get gene information ----------------------------------------------------
mart <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl",mart)
mart.df.gn <- getBM(attributes = c('gene_biotype', 'ensembl_gene_id'),
                 filters = 'ensembl_gene_id',
                 values = gn$gene_id,
                 mart = mart)
mart.df.gn <- mart.df.gn[match(gn$gene_id, mart.df.gn$ensembl_gene_id),]
rownames(mart.df.gn) <- NULL

unique(mart.df.gn$gene_biotype)

qplot(width(gn[!grepl('pseudo', mart.df.gn$gene_biotype)]), bins=100, color=I('white')) +
  scale_x_log10() +
  theme_bw()

qplot(width(gn[mart.df.gn$gene_biotype == 'protein_coding']), color=I('white'), bins=100) +
  scale_x_log10() +
  theme_bw()

summary(width(gn[mart.df.gn$gene_biotype == 'protein_coding']))



# Transcript information --------------------------------------------------
mart.df.tr <- getBM(attributes = c('gene_biotype', 'ensembl_gene_id', 'ensembl_transcript_id'),
                 filters = 'ensembl_transcript_id',
                 values = tr$tx_name,
                 mart = mart)
mart.df.tr <- mart.df.tr[match(gn$gene_id, mart.df.tr$ensembl_gene_id),]
rownames(mart.df.tr) <- NULL

unique(mart.df.tr$gene_biotype)

qplot(width(tr[!grepl('pseudo', mart.df.tr$gene_biotype)]), bins=100, color=I('white')) +
  scale_x_log10() +
  theme_bw()

qplot(width(tr[mart.df.tr$gene_biotype == 'protein_coding']), color=I('white'), bins=100) +
  scale_x_log10() +
  theme_bw()

summary(width(tr[mart.df.tr$gene_biotype == 'protein_coding']))


# Get sequences -----------------------------------------------------------
sq.gn <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gn, as.character=T)
write_tsv(tibble(sequence=sq.gn), 'sequences_gene.seq', col_names=F)
write_tsv(tibble(ids=gn$gene_id), 'sequences_gene.id.tsv', col_names=F)

sq.tr <- getSeq(BSgenome.Hsapiens.UCSC.hg38, tr, as.character=T)
write_tsv(tibble(sequence=sq.tr), 'sequences_transcript.seq', col_names=F)
write_tsv(tibble(tx_id=tr$tx_name,
                 gene_id=unlist(tr$gene_id)), 'sequences_transcript.id.tsv')

# sq.cd <- getSeq(BSgenome.Hsapiens.UCSC.hg38, cd, as.character=T)

tr.sub <- tr[(width(tr) > 1000 & width(tr) < 10000)]
tr.sub.genes <- unique(unlist(tr.sub$gene_id))
tr.sub.validation.genes <- tr.sub.genes[sample(seq_along(tr.sub.genes), length(tr.sub.genes)*0.1)]
tr.sub.validation <- tr.sub[unlist(tr.sub$gene_id) %in% tr.sub.validation.genes]
tr.sub.training <- tr.sub[!(unlist(tr.sub$gene_id) %in% tr.sub.validation.genes)]

sq.tr.sub.training <- getSeq(BSgenome.Hsapiens.UCSC.hg38, tr.sub.training, as.character=T)
sq.tr.sub.validation <- getSeq(BSgenome.Hsapiens.UCSC.hg38, tr.sub.validation, as.character=T)

write_tsv(tibble(sequence=sq.tr.sub.training), 'sequences_transcript_sub_training.seq', col_names=F)
write_tsv(tibble(tx_id=tr.sub.training$tx_name,
                 gene_id=unlist(tr.sub.training$gene_id)), 'sequences_transcript_sub_training.id.tsv')

write_tsv(tibble(sequence=sq.tr.sub.validation), 'sequences_transcript_sub_validation.seq', col_names=F)
write_tsv(tibble(tx_id=tr.sub.validation$tx_name,
                 gene_id=unlist(tr.sub.validation$gene_id)), 'sequences_transcript_sub_validation.id.tsv')


gc()