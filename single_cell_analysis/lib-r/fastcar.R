library(Seurat)
library(FastCAR)
library(Matrix)
library(qlcMatrix)
library(R.utils)


options(error = function() {traceback(3); q(1)})
options(stringsAsFactors = FALSE)


writeMtx <- function(dir, counts) {
    matrixPath <- file.path(dir, "matrix.mtx")
    genesPath <- file.path(dir, "features.tsv")
    barcodesPath <- file.path(dir, "barcodes.tsv")

    writeMM(counts, matrixPath)
    # easier to load if the genes file has at least columns even though seurat objects
    # don't have yet explicit geneIds/geneSyms data, we just duplicate whatever the matrix has now
    write.table(as.data.frame(cbind(
        rownames(counts),
        rownames(counts),
        rep("Gene Expression", nrow(counts))
    )), file=genesPath, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
    write(colnames(counts), file = barcodesPath)

    gzip(matrixPath)
    gzip(genesPath)
    gzip(barcodesPath)
}


renameGenes <- function(genes) {
    genes <- gsub("^GRCh38_+", "", genes)
    genes <- gsub("^SARS-CoV-2i_", "SARS-CoV-2-", genes)
    genes <- gsub("SARS-CoV-2-antisense", "Antisense", genes)
    return(genes)
}


process <- function(cell.dir, full.dir, output, cutoff, prob) {
    if (endsWith(cell.dir, ".h5")) {
        cell.mtx <- Read10X_h5(cell.dir)
        full.mtx <- Read10X_h5(full.dir)
    } else {
        cell.mtx <- read.cell.matrix(cell.dir)
        full.mtx <- read.full.matrix(full.dir)
    }
    backgrnd <- determine.background.to.remove(full.mtx, cell.mtx, cutoff, prob)

    to.corr <- sort(backgrnd[backgrnd > 0], decreasing = TRUE)
    names(to.corr) <- renameGenes(names(to.corr))
    report <- file(file.path(output, "fastcar.txt"), "w")
    write(sprintf("Params: cutoff=%d, prob=%f", cutoff, prob), file = report)
    for (gene in names(to.corr)) {
        write(sprintf("%s: %d", gene, to.corr[gene]), file = report)
    }
    close(report)

    new.mtx <- remove.background(cell.mtx, backgrnd)
    writeMtx(output, new.mtx)
}


args <- unlist(commandArgs(trailingOnly = TRUE))
process(
    args[1],
    args[2],
    args[3],
    as.numeric(args[4]),
    as.numeric(args[5])
)
