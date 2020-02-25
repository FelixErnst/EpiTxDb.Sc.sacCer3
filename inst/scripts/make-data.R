# get annotation data
library(AnnotationHub)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
library(EpiTxDb)
library(RSQLite)
library(GenomicRanges)

# get annotation hub
ah <- AnnotationHub()

# get ensembl annotation
edb <- query(ah, c("EnsDb","Saccharomyces cerevisiae", "99"))[[1]]
seqlevelsStyle(edb) <- "UCSC"

# get transcript annotations from ensemble db
assemble_tx <- function(edb, genome){
    # get exons
    tx <- exonsBy(edb,"tx")
    tx_id <- IRanges::CharacterList(Map(rep,names(tx),lengths(tx)))
    mcols(tx, level="within")[,"tx_id"] <- tx_id
    genome(tx) <- genome
    #
    tx[lengths(tx) != 0L]
}

tx <- assemble_tx(edb, "sacCer3")

################################################################################
# functions for import
################################################################################

import.RMBase <- function(bs, organism, genome, type){
    metadata <- data.frame(name = c("Data source","Organism","Genome",
                                    "Coordinates"),
                           value = c("RMBase v2.0","Saccharomyces cerevisiae",
                                     "sacCer3","per Genome"))
    #
    seq <- getSeq(bs)
    seq_rna <- as(seq,"RNAStringSet")
    #
    makeEpiTxDbFromRMBase(organism, genome, type, sequences = seq_rna,
                          metadata = metadata)
}

import_from_tRNAdb <- function(bs, tx){
  metadata <- data.frame(name = c("Data source","Organism","Genome",
                                  "Coordinates"),
                         value = c("tRNAdb","Saccharomyces cerevisiae",
                                   "sacCer3","per Genome"))
    #
    seq <- getSeq(bs,tx)
    seq <- relist(unlist(unlist(seq)),
                  IRanges::PartitioningByWidth(sum(nchar(seq))))
    seq_rna <- as(seq,"RNAStringSet")
    makeEpiTxDbFromtRNAdb("Saccharomyces cerevisiae", seq_rna,
                          metadata = metadata)
}

start.import <- function(bs, tx){
    etdb <- import.RMBase(bs, "yeast", "sacCer3",
                          listAvailableModFromRMBase("yeast", "sacCer3"))
    db <- dbConnect(SQLite(), "hub/EpiTxDb.Sc.sacCer3.RMBase.sqlite")
    sqliteCopyDatabase(etdb$conn, db)
    dbDisconnect(etdb$conn)
    dbDisconnect(db)
  
    etdb <- import_from_tRNAdb(bs, tx)
    db <- dbConnect(SQLite(), "hub/EpiTxDb.Sc.sacCer3.tRNAdb.sqlite")
    sqliteCopyDatabase(etdb$conn, db)
    dbDisconnect(etdb$conn)
    dbDisconnect(db)
    return(TRUE)
}

start.import(BSgenome.Scerevisiae.UCSC.sacCer3, tx)
