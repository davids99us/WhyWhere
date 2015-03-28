test_that("wrapper runs on example code", {
  o=whywhere(Bradypus_Pres,Bradypus_files)
  expect_more_than(o$result[1]$AUC,0.6)
  plot.model(o)
  plot.membership(o)
}) 