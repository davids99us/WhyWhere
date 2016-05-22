 
test_that("wrapper runs on example code", {
  o=ww(Bradypus_Pres,Bradypus_files)
  expect_more_than(o$result[1]$AUC,0.6)
  plot.ww(o)
  plot.dseg(o)
}) 