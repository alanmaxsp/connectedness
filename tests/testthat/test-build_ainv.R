test_that("build_Ainv devuelve Ainv sparse y F con dimensiones correctas", {
  ped <- data.frame(
    animal = c("A", "B", "C", "D"),
    sire   = c("0", "0", "A", "A"),
    dam    = c("0", "0", "B", "B"),
    stringsAsFactors = FALSE
  )
  
  ren <- connectedness::renum_pedigree(ped, verbose = FALSE)
  out <- connectedness::build_Ainv(ren)
  
  expect_type(out, "list")
  expect_true(all(c("Ainv", "F") %in% names(out)))
  expect_s4_class(out$Ainv, "dgCMatrix")
  expect_equal(nrow(out$Ainv), nrow(ren))
  expect_equal(ncol(out$Ainv), nrow(ren))
  expect_length(out$F, nrow(ren))
})

test_that("build_Ainv falla si faltan columnas requeridas", {
  bad <- data.frame(new_id = 1:3, new_sire = c(0,0,1))
  expect_error(
    connectedness::build_Ainv(bad),
    "must have columns: new_id, new_sire, new_dam"
  )
})