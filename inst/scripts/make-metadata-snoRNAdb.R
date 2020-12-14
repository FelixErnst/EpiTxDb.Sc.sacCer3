
# Base data for all data sets --------------------------------------------------
library(S4Vectors)

df_Base <- DataFrame(
  BiocVersion = "3.12",
  SourceVersion = NA,
  Coordinate_1_based = TRUE,
  Maintainer = "Felix G.M. Ernst <felix.gm.ernst@outlook.com>"
)

snoRNAdbURL <- "https://www-snorna.biotoul.fr/"

df <- rbind(
  cbind(df_Base,
        DataFrame(Title = "EpiTxDb snoRNAdb for Saccharomyces cerevisisae sacCer3", 
                  Description = paste0(
                    "Information from the snoRNAdb was downloaded, updated ",
                    "with data from the Fournier LAB at UMASS and stored as ",
                    "EpiTxDb database. The ",
                    "information provided match sacCer3 release sequences."), 
                  SourceType = "BED",
                  SourceUrl = snoRNAdbURL,
                  DataProvider = "RMBase v2.0",
                  RDataClass = "SQLiteFile", 
                  DispatchClass = "SQLiteFile",
                  RDataPath = "EpiTxDb.Sc.sacCer3/EpiTxDb.Sc.sacCer3.snoRNAdb.sqlite")),
  cbind(df_Base,
        DataFrame(Title = "Sequences of snoRNA targets of Saccharomyces cerevisisae sacCer3", 
                  Description = paste0(
                    "Fasta file for snoRNA targets based on genomic sequences ",
                    "for Saccharomyces cerevisisae sacCer3."),
                  SourceType = "FASTA",
                  SourceUrl = "https://www.ncbi.nlm.nih.gov/gene",
                  DataProvider = "NCBI",
                  RDataClass = "FaFile", 
                  DispatchClass = "FaFile",
                  RDataPath = "EpiTxDb.Sc.sacCer3/snoRNA.targets.sacCer3.fa"))
)

df$Species <- "Saccharomyces cerevisiae S288C"
df$TaxonomyId <- "559292"
df$SourceVersion <- Sys.time()
df$Genome <- "sacCer3"
df$Tags <- "EpiTxDb:sacCer3:Modification:Epitranscriptomics"

write.csv(df, file = "inst/extdata/metadata-snoRNAdb.csv", row.names = FALSE)
