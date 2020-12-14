# get annotation data
library(AnnotationHub)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
library(org.Sc.sgd.db)
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
                                   "sacCer3","per Transcript"))
    #
    seq <- getSeq(bs,tx)
    seq <- relist(unlist(unlist(seq)),
                  IRanges::PartitioningByWidth(sum(nchar(seq))))
    seq_rna <- as(seq,"RNAStringSet")
    makeEpiTxDbFromtRNAdb("Saccharomyces cerevisiae", seq_rna,
                          metadata = metadata)
}

import_from_snoRNAdb <- function(snoRNAdb, orgdb){
  metadata <- data.frame(name = c("Data source","Organism","Genome",
                                  "Coordinates"),
                         value = c("snoRNAdb","Saccharomyces cerevisiae","sacCer3",
                                   "per Transcript"))
  # Modifications
  mod_id <- seq_len(nrow(snoRNAdb))
  mod_name <- paste0(snoRNAdb$mod,"_",snoRNAdb$position)
  mod_type <- snoRNAdb$mod
  mod_start <- snoRNAdb$position
  mod_end <- snoRNAdb$position
  
  transcripts <- select(orgdb,as.character(snoRNAdb$seqnames),
                        c("REFSEQ","SGD","ENTREZID"),"GENENAME")
  
  modifications <- data.frame(mod_id = mod_id,
                              mod_name = mod_name,
                              mod_type = mod_type,
                              mod_start = mod_start,
                              mod_end = mod_end,
                              mod_strand = "+",
                              sn_id = as.integer(transcripts$ENTREZID),
                              sn_name = transcripts$REFSEQ,
                              stringsAsFactors = FALSE)
  
  # Reactions
  gene_fbl <- select(orgdb,keys = "NOP1",
                     columns = c("GENENAME","ENSEMBL","ENTREZID"),
                     keytype = "GENENAME")
  
  gene_dkc <- select(orgdb,keys = "CBF5",
                     columns = c("GENENAME","ENSEMBL","ENTREZID"),
                     keytype = "GENENAME")
  
  gene_pus7p <- select(orgdb,keys = "PUS7",
                       columns = c("GENENAME","ENSEMBL","ENTREZID"),
                       keytype = "GENENAME")
  
  gene_pus1p <- select(orgdb,keys = "PUS1",
                       columns = c("GENENAME","ENSEMBL","ENTREZID"),
                       keytype = "GENENAME")
  
  rx_rank <- 1L
  mod_type <- snoRNAdb$mod
  genename <- character(length(mod_type))
  ensembl <- character(length(mod_type))
  ensembltrans <- character(length(mod_type))
  entrezid <- character(length(mod_type))
  
  genename[mod_type == "Y" & snoRNAdb$modifier != "Pus7p"] <- gene_dkc$GENENAME
  genename[mod_type %in% c("Am","Gm","Cm","Um")] <- gene_fbl$GENENAME
  genename[snoRNAdb$modifier == "Pus7p"] <- gene_pus7p$GENENAME
  genename[snoRNAdb$modifier == "Pus1p"] <- gene_pus1p$GENENAME
  
  ensembl[mod_type == "Y" & snoRNAdb$modifier != "Pus7p"] <- gene_dkc$ENSEMBL
  ensembl[mod_type %in% c("Am","Gm","Cm","Um")] <- gene_fbl$ENSEMBL
  ensembl[snoRNAdb$modifier == "Pus7p"] <- gene_pus7p$ENSEMBL
  ensembl[snoRNAdb$modifier == "Pus1p"] <- gene_pus1p$ENSEMBL
  
  entrezid[mod_type == "Y" & snoRNAdb$modifier != "Pus7p"] <- gene_dkc$ENTREZID
  entrezid[mod_type %in% c("Am","Gm","Cm","Um")] <- gene_fbl$ENTREZID
  entrezid[snoRNAdb$modifier == "Pus7p"] <- gene_pus7p$ENTREZID
  entrezid[snoRNAdb$modifier == "Pus1p"] <- gene_pus1p$ENTREZID
  
  reactions <- data.frame(mod_id = mod_id,
                          rx_genename = genename,
                          rx_rank = rx_rank,
                          rx_ensembl = ensembl,
                          rx_ensembltrans = ensembltrans,
                          rx_entrezid = entrezid,
                          stringsAsFactors = FALSE)
  
  # Specifiers
  specifier_genename <- snoRNAdb$modifier
  specifier_genename <- toupper(specifier_genename)
  specifier_f <- !(specifier_genename %in% c("unknown",""))
  specifier_type <- rep("snoRNA",length(specifier_genename))
  specifier_type[specifier_genename==""] <- NA
  specifier_type[!grepl("^SNR",specifier_genename)] <- "protein"
  specifier_genename <- strsplit(as.character(specifier_genename),",")[specifier_f]
  specifier_genename[specifier_genename == "PUS7P"] <- "PUS7"
  specifier_genename[specifier_genename == "PUS1P"] <- "PUS1"
  specifier_genename[specifier_genename == "SPB1P"] <- "SPB1"
  specifier_lengths <- lengths(specifier_genename)
  specifier_mod_id <- unlist(Map(rep,mod_id[specifier_f],specifier_lengths))
  specifier_type <- unlist(Map(rep,specifier_type[specifier_f],specifier_lengths))
  specifier_entrezid <- mapIds(orgdb,unlist(specifier_genename),"ENTREZID",
                               "GENENAME")
  specifier_ensembl <- mapIds(orgdb,unlist(specifier_genename),"ENSEMBL",
                              "GENENAME")
  
  specifiers <- data.frame(mod_id = specifier_mod_id,
                           spec_type = specifier_type,
                           spec_genename = unlist(specifier_genename),
                           spec_ensembl = specifier_ensembl,
                           spec_entrezid = specifier_entrezid,
                           stringsAsFactors = FALSE)
  # References
  references <- data.frame(mod_id = mod_id,
                           ref_type = "PMID",
                           ref = as.character(snoRNAdb$pmid))
  
  # clean up
  rm_mod_id <- modifications[modifications$mod_type %in% c("m1Y","acp3Y"), "mod_id"]
  specifiers <- specifiers[!(specifiers$mod_id %in% rm_mod_id),]
  reactions <- reactions[!(reactions$mod_id %in% rm_mod_id),]
  makeEpiTxDb(modifications, reactions, specifiers, references,
              metadata = metadata)
}

start.import <- function(bs, tx, orgdb){
    etdb <- import_from_snoRNAdb(read.csv("D:/AWS/EpiTxDb.Sc.sacCer3/yeast_snoRNA.csv",sep = ";",fileEncoding = "UTF-8-BOM"),
                                 orgdb)
    db <- dbConnect(SQLite(), "hub/EpiTxDb.Sc.sacCer3.snoRNAdb.sqlite")
    sqliteCopyDatabase(etdb$conn, db)
    dbDisconnect(etdb$conn)
    dbDisconnect(db)
  
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

start.import(BSgenome.Scerevisiae.UCSC.sacCer3, tx, org.Sc.sgd.db)
