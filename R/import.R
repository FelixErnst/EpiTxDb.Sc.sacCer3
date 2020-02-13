# 
# YEASTGENOME_SACCER3_RELEASE <- "https://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_R64-2-1_20150113.tgz"
# YEASTGENOME_SACCER3_GFF <- "S288C_reference_genome_R64-2-1_20150113/saccharomyces_cerevisiae_R64-2-1_20150113.gff"
# 
# .download_sacCer3_gff <- function(){
#   dir <- tempdir()
#   file <- paste0(dir,"\\RELEASE")
#   ans <- try(curl::curl_download(YEASTGENOME_SACCER3_RELEASE,file),
#              silent = TRUE)
#   if(is(ans,"try-error")){
#     stop("Download of release tar unsuccessful.")
#   }
#   file_untared <- untar(file, YEASTGENOME_SACCER3_GFF, exdir = dir)
#   if(file_untared != 0){
#     stop("Extraction of release tar unsuccessful.")
#   }
#   gff <- rtracklayer::import.gff3(paste0(dir,"/",YEASTGENOME_SACCER3_GFF))
#   gff <- gff[is.na(mcols(gff)$orf_classification) |
#                mcols(gff)$orf_classification != "Dubious"]
#   gff <- gff[as.vector(strand(gff)) %in% c("+","-")]
#   gff <- GRanges(seqnames = seqnames(gff),
#                  ranges = ranges(gff),
#                  strand = strand(gff),
#                  mcols(gff)[,c("type","ID","Name","Parent","gene")])
#   gff
# }
# 
# .assemble_tx <- function(types = c("rRNA_gene","tRNA_gene","snRNA_gene",
#                                    "ncRNA_gene","snoRNA_gene",
#                                    "transposable_element_gene")){
#   gff <- .download_sacCer3_gff()
#   # get protein encoding transcripts
#   genes <- mcols(gff)$type == "gene"
# 
#   hits_mRNA <- findMatches(gff$Parent,gff[genes]$ID)
#   hits_mRNA <- hits_mRNA[queryHits(hits_mRNA) %in% which(gff$type == "mRNA")]
# 
#   hits_CDS <- findMatches(gff$Parent,gff[queryHits(hits_mRNA)]$ID)
#   hits_CDS <- hits_CDS[queryHits(hits_CDS) %in% which(gff$type == "CDS")]
# 
#   tx <- do.call(pc,list(split(gff[genes],subjectHits(hits_mRNA)),
#                         split(gff[queryHits(hits_mRNA)],subjectHits(hits_mRNA)),
#                         split(gff[queryHits(hits_CDS)],subjectHits(hits_mRNA)[subjectHits(hits_CDS)])))
#   names(tx) <- unlist(gff[queryHits(hits_mRNA)]$ID)
#   mcols <- mcols(tx,level="within")
#   stopifnot(all(names(tx) == unlist(mcols[mcols[,"type"] == "mRNA"][,"ID"])))
# 
#   df <- DataFrame(gene_id = gff[genes]$ID,
#                   gene_name = gff[genes]$gene,
#                   transcript_name = unlist(gff[queryHits(hits_mRNA)]$ID))
#   mcols(tx) <- df
#   stopifnot(all(names(tx) == df$transcript_name))
# 
#   # get non coding transcripts
#   non_coding <- mcols(gff)$type %in% types
# 
#   hits_nc <- findMatches(gff$Parent,gff[non_coding]$ID)
#   hits_nc <- hits_nc[queryHits(hits_nc) %in% which(gff$type %in% c("noncoding_exon","intron","CDS","plus_1_translational_frameshift"))]
# 
#   nc_tx <- gff[non_coding]
#   gene_name <- nc_tx$ID
#   gene_name[!is.na(mcols(nc_tx)$Name)] <- nc_tx[!is.na(mcols(nc_tx)$Name)]$Name
#   gene_name[!is.na(mcols(nc_tx)$gene)] <- nc_tx[!is.na(mcols(nc_tx)$gene)]$gene
#   transcript_name <- nc_tx$Name
#   df <- DataFrame(gene_id = nc_tx$ID,
#                   gene_name = gene_name,
#                   transcript_name = nc_tx$Name)
# 
#   nc_tx <- split(nc_tx,seq_along(nc_tx))
#   names(nc_tx) <- transcript_name
#   nc_tx[unique(subjectHits(hits_nc))] <-
#     pc(nc_tx[unique(subjectHits(hits_nc))],
#        split(gff[queryHits(hits_nc)],subjectHits(hits_nc)))
# 
#   mcols(nc_tx) <- df
#   stopifnot(all(names(nc_tx) == df$transcript_name))
#   # combine
#   c(tx,nc_tx)
# }
# 
# import.RMBase <- function(bs, tx, organism, genome, type){
#   seq <- getSeq(bs,tx)
#   seq <- relist(unlist(unlist(seq)),
#                 IRanges::PartitioningByWidth(sum(nchar(seq))))
#   seq_rna <- as(seq,"RNAStringSet")
#   #
#   makeEpiTxDbfromRMBase(organism, genome, type, tx, seq_rna)
# }
# 
# import_from_tRNAdb <- function(bs, tx){
#   seq <- getSeq(bs,tx)
#   seq <- relist(unlist(unlist(seq)),
#                 IRanges::PartitioningByWidth(sum(nchar(seq))))
#   seq_rna <- as(seq,"RNAStringSet")
#   makeEpiTxDbFromtRNAdb("Saccharomyces cerevisiae", tx, seq_rna)
# }
# 
# start.import <- function(bs, tx){
#   tx_sub <- tx[mcols(tx,level="within")[,"type"] %in% c("CDS","noncoding_exon")]
#   tx_sub <- tx_sub[lengths(tx_sub) != 0L]
#   etdb <- import.RMBase(bs, tx_sub, "yeast",
#                         "sacCer3",
#                         listAvailableModFromRMBase("yeast", "sacCer3"))
#   db <- dbConnect(SQLite(), "inst/extdata/EpiTxDb.Sc.sacCer3.RMBase.sqlite")
#   sqliteCopyDatabase(etdb$conn, db)
#   dbDisconnect(etdb$conn)
#   dbDisconnect(db)
#   
#   etdb <- import_from_tRNAdb(bs, tx_sub)
#   db <- dbConnect(SQLite(), "inst/extdata/EpiTxDb.Sc.sacCer3.tRNAdb.sqlite")
#   sqliteCopyDatabase(etdb$conn, db)
#   dbDisconnect(etdb$conn)
#   dbDisconnect(db)
#   return(TRUE)
# }
# 
# library(BSgenome.Scerevisiae.UCSC.sacCer3)
# library(EpiTxDb)
# library(RSQLite)
# library(GenomicRanges)
# 
# tx <- .assemble_tx()
# seqlevels(tx) <- gsub("chrmt","chrM",seqlevels(tx))
# 
# start.import(BSgenome.Scerevisiae.UCSC.sacCer3, tx)
