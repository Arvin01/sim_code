library(dplyr)
tissuedat <- read.csv("../../../data/gtex_tissue_gene_reads/heart.csv")
exon <- read.table("../../../data/gencode.v19.genes.patched_contigs_exons.txt", header = TRUE)

## Associate genes with chromosome position -----------------------------------------
gene_names <- rep(NA, length = length(exon$Gene))
for (index in 1:length(gene_names)) {
    gene_names[index] <- sub("([[:alnum:]]+\\.[[:digit:]]+)\\_[[:digit:]]+", "\\1",
                             exon$Gene[index])
    if (index%%10000 == 0) {
        cat("Iteration", index, "\n")
    }
}
## exon$Gene_trunc <- gene_names
## write.csv("exon", file = "exon.csv", row.names = FALSE)
## match_out <- match(gene_names, tissuedat$Name)
## all(sort(unique(match_out)) == 1:length(unique(match_out)))

match_out <- match(tissuedat$Name, gene_names)
tissuedat$CHR <- exon$CHR[match_out]
rm(exon)
## ----------------------------------------------------------------------------------

## Select top expressed 10000 genes for analysis -----------------------------------
expr <- select(tissuedat, -c(Name, Description, CHR))
gene_levels <- rowSums(log(expr + 1))
genes_order <- order(gene_levels, decreasing = TRUE)
subtissue <- tissuedat[genes_order[1:10000], ]

rm(tissuedat)
rm(expr)
## ----------------------------------------------------------------------------------


## get subject genders
pheno <- read.table("../../../data/GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt",
                    header = TRUE)
