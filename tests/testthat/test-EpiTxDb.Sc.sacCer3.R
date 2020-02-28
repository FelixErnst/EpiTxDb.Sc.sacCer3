context("EpiTxDb.Sc.sacCer3")
test_that("EpiTxDb.Sc.sacCer3:",{
    etdb <- EpiTxDb.Sc.sacCer3.tRNAdb()
    expect_s4_class(etdb,"EpiTxDb")
})
