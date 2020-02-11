
YEASTGENOME_SACCER3_RELEASE <- "https://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_R64-2-1_20150113.tgz"
YEASTGENOME_SACCER3_GFF <- "S288C_reference_genome_R64-2-1_20150113/saccharomyces_cerevisiae_R64-2-1_20150113.gff"

.download_sacCer3_gff <- function(){
  dir <- tempdir()
  file <- paste0(dir,"\\RELEASE")
  ans <- try(curl::curl_download(YEASTGENOME_SACCER3_RELEASE,file),
             silent = TRUE)
  if(is(ans,"try-error")){
    stop("Download of release tar unsuccessful.")
  }
  file_untared <- untar(file, YEASTGENOME_SACCER3_GFF, exdir = dir)
  if(file_untared != 0){
    stop("Extraction of release tar unsuccessful.")
  }
  gff <- rtracklayer::import.gff3(paste0(dir,"/",YEASTGENOME_SACCER3_GFF))
  gff <- gff[is.na(mcols(gff)$orf_classification) |
               mcols(gff)$orf_classification != "Dubious"]
  gff <- gff[as.vector(strand(gff)) %in% c("+","-")]
  gff <- GRanges(seqnames = seqnames(gff),
                 ranges = ranges(gff),
                 strand = strand(gff),
                 mcols(gff)[,c("type","ID","Name","Parent","gene")])
  gff
}

.assemble_tx <- function(types = c("rRNA_gene","tRNA_gene","snRNA_gene",
                                   "ncRNA_gene","snoRNA_gene",
                                   "transposable_element_gene")){
  gff <- .download_sacCer3_gff()
  # get protein encoding transcripts
  genes <- mcols(gff)$type == "gene"
  children <- lengths(mcols(gff)$Parent) != 0L
  m <- match(gff[children]$Parent,gff[genes]$ID)
  f <- unlist(!is.na(m)) & gff[children]$type == "mRNA"
  mRNA <- f
  m <- match(gff[children]$Parent,gff[children][mRNA]$ID)
  f <- unlist(!is.na(m)) & gff[children]$type == "CDS"
  CDS <- f
  tx <- split(gff[children][CDS],factor(unlist(gff[children][CDS]$Parent),
                                        unique(unlist(gff[children][CDS]$Parent))))
  df <- DataFrame(gene_id = gff[genes]$ID,
                  gene_name = gff[genes]$gene,
                  transcript_name = unlist(gff[children][mRNA]$ID))
  df <- df[match(names(tx),df$transcript_name),]
  mcols(tx) <- df
  # get non coding transcripts
  RNA_tx <- gff[mcols(gff)$type %in% types]
  gene_name <- RNA_tx$gene
  gene_name[!is.na(mcols(RNA_tx)$Name)] <- RNA_tx[!is.na(mcols(RNA_tx)$Name)]$Name
  gene_name[!is.na(mcols(RNA_tx)$gene)] <- RNA_tx[!is.na(mcols(RNA_tx)$gene)]$gene
  df <- DataFrame(gene_id = RNA_tx$ID,
                  gene_name = gene_name,
                  transcript_name = RNA_tx$Name)
  RNA_tx <- split(RNA_tx, RNA_tx$Name)
  df <- df[match(names(RNA_tx),df$transcript_name),]
  mcols(RNA_tx) <- df
  # combine
  c(tx,RNA_tx)
}

import.RMBase <- function(bs, tx, organism, genome, type){
  files <- downloadRMBaseFiles(organism, genome, type)
  grl <- getRMBaseData(files, GenomeInfoDb::seqlevels(tx))
  gr <- unlist(GenomicRanges::GRangesList(grl))
  gr <- shiftToTranscriptCoordinates(gr, tx)
  mcols(gr)$mod_id <- seq_along(gr)
  colnames(mcols(gr)) <- gsub("mod_type","mod",colnames(mcols(gr)))
  # Cleanup base
  seq <- getSeq(bs,tx)
  seq <- relist(unlist(unlist(seq)),
                IRanges::PartitioningByWidth(sum(nchar(seq))))
  seq_rna <- as(seq,"RNAStringSet")
  gr <- Modstrings::removeIncompatibleModifications(gr, seq_rna)
  colnames(mcols(gr)) <- gsub("^mod$","mod_type",colnames(mcols(gr)))
  #
  makeEpiTxDbfromGRanges(gr)
}

import_from_tRNAdb <- function(bs, tx){
  seq <- getSeq(bs,tx)
  seq <- relist(unlist(unlist(seq)),
                IRanges::PartitioningByWidth(sum(nchar(seq))))
  seq_rna <- as(seq,"RNAStringSet")
  makeEpiTxDbFromtRNAdb("Saccharomyces cerevisiae", tx, seq_rna)
}

start.import <- function(tx){
  etdb <- import.RMBase(BSgenome.Scerevisiae.UCSC.sacCer3, tx, "yeast",
                        "sacCer3",
                        listAvailableModFromRMBase("yeast", "sacCer3"))
  db <- dbConnect(SQLite(), "inst/extdata/EpiTxDb.Sc.sacCer3.RMBase.sqlite")
  sqliteCopyDatabase(etdb$conn, db)
  dbDisconnect(etdb$conn)
  dbDisconnect(db)
  etdb <- import_from_tRNAdb(BSgenome.Scerevisiae.UCSC.sacCer3, tx)
  db <- dbConnect(SQLite(), "inst/extdata/EpiTxDb.Sc.sacCer3.tRNAdb.sqlite")
  sqliteCopyDatabase(etdb$conn, db)
  dbDisconnect(etdb$conn)
  dbDisconnect(db)
  return(TRUE)
}

library(BSgenome.Scerevisiae.UCSC.sacCer3)
library(EpiTxDb)
library(RSQLite)
library(GenomicRanges)

tx <- .assemble_tx()
seqlevels(tx) <- gsub("chrmt","chrM",seqlevels(tx))

start.import(tx)
