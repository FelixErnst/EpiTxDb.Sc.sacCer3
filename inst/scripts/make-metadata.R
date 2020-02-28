
# Base data for all data sets --------------------------------------------------
library(S4Vectors)

df_Base <- DataFrame(
  BiocVersion = "3.11",
  SourceVersion = NA,
  Coordinate_1_based = TRUE,
  Maintainer = "Felix G.M. Ernst <felix.gm.ernst@outlook.com>"
)

RMBaseURL <- "http://rna.sysu.edu.cn/rmbase/"
tRNAdbURL <- "http://trna.bioinf.uni-leipzig.de/"

df <- rbind(
  cbind(df_Base,
        DataFrame(Title = "EpiTxDb RMBase v2.0 for Saccharomyces cerevisiae sacCer3", 
                  Description = paste0(
                    "Information from the RMBase v2.0 was downloaded and ",
                    "imported as EpiTxDb database. All valid modification ",
                    "types for Saccharomyces cerevisiae/sacCer3 were used and ",
                    "checked the nucleotide sequence."), 
                  SourceType = "BED",
                  SourceUrl = RMBaseURL,
                  DataProvider = "RMBase v2.0",
                  RDataClass = "SQLiteFile", 
                  DispatchClass = "SQLiteFile",
                  RDataPath = "EpiTxDb.Sc.sacCer3/EpiTxDb.Sc.sacCer3.RMBase.sqlite")),
  cbind(df_Base,
        DataFrame(Title = "EpiTxDb tRNAdb for Saccharomyces cerevisiae sacCer3", 
                  Description = paste0(
                    "tRNAdb entries for RNA database were imported using ",
                    "tRNAdbImport and the modification information extracted. ",
                    "tRNA sequences were matched against transcripts, ",
                    "results combined and stored as EpiTxDb database."),
                  SourceType = "XML",
                  SourceUrl = tRNAdbURL,
                  DataProvider = "tRNAdb",
                  RDataClass = "SQLiteFile", 
                  DispatchClass = "SQLiteFile",
                  RDataPath = "EpiTxDb.Sc.sacCer3/EpiTxDb.Sc.sacCer3.tRNAdb.sqlite"))
)

df$Species <- "Saccharomyces cerevisiae S288C"
df$TaxonomyId <- "559292"
df$SourceVersion <- Sys.time()
df$Genome <- "sacCer3"
df$Tags <- "EpiTxDb:sacCer3:Modification:Epitranscriptomics"

write.csv(df, file = "inst/extdata/metadata.csv", row.names = FALSE)
