

test_that("DataModel", {
  data(array, package = "BattleField")
  model <- DataModel(array)
  expect_s4_class(model,"BattleField")
  
})
