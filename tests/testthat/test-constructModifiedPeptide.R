library(MCMS)

context("constructModifiedPeptide")

test_that("constructModifiedPeptide deals with unmodified peptides", {
	p = constructModifiedPeptide("_AGDSDSWK_", 1)
	expect_equal(p, "")
})

test_that("constructModifiedPeptide deals with one or several phospho only", {
	p = constructModifiedPeptide("_AGDS(ph)DSWK_", 1)
	expect_equal(p, "S_4")
	p = constructModifiedPeptide("_AGDSDS(ph)WK_", 1)
	expect_equal(p, "S_6")
	p = constructModifiedPeptide("_AGDS(ph)DS(ph)WK_", 1)
	expect_equal(p, "S_4;S_6")
})

test_that("constructModifiedPeptide deals with ox only", {
	p = constructModifiedPeptide("_AAM(ox)VYMFQYDSTHGK_", 1)
	expect_equal(p, "M_3")
	p = constructModifiedPeptide("_AAMVYM(ox)FQYDSTHGK_", 1)
	expect_equal(p, "M_6")
	p = constructModifiedPeptide("_AAM(ox)VYM(ox)FQYDSTHGK_", 1)
	expect_equal(p, "M_3;M_6")
})

test_that("constructModifiedPeptide deals with acetylation only", {
	p = constructModifiedPeptide("_(ac)AAMVYMFQYDSTHGK_", 1)
	expect_equal(p, "a")
})


test_that("constructModifiedPeptide deals with both phospho, acetylation and ox", {
	p = constructModifiedPeptide("_(ac)AAEM(ox)CY(ph)RK_", 1)
	expect_equal(p, "a;M_4;Y_6")
})


test_that("constructModifiedPeptide deals with both phospho and ox also at start and end positions", {
	p = constructModifiedPeptide("_M(ox)CY(ph)RKS(ph)_", 1)
	expect_equal(p, "M_1;Y_3;S_6")
	p = constructModifiedPeptide("_T(ph)CY(ph)RKM(ox)_", 1)
	expect_equal(p, "T_1;Y_3;M_6")
	# And with ac?
	p = constructModifiedPeptide("_(ac)M(ox)CY(ph)RKS(ph)_", 1)
	expect_equal(p, "a;M_1;Y_3;S_6")
	p = constructModifiedPeptide("_(ac)T(ph)CY(ph)RKM(ox)_", 1)
	expect_equal(p, "a;T_1;Y_3;M_6")
})


test_that("constructModifiedPeptide handles delimiter", {
	p = constructModifiedPeptide("_(ac)M(ox)CY(ph)RKS(ph)_", 1, "+")
	expect_equal(p, "a+M_1+Y_3+S_6")
})
test_that("constructModifiedPeptide handles different positions", {
	p = constructModifiedPeptide("_(ac)M(ox)CY(ph)RKS(ph)_", 10)
	expect_equal(p, "a;M_10;Y_12;S_15")
	p = constructModifiedPeptide("_(ac)M(ox)CY(ph)RKS(ph)_", 12)
	expect_equal(p, "a;M_12;Y_14;S_17")
	# Also with delim
	p = constructModifiedPeptide("_(ac)M(ox)CY(ph)RKS(ph)_", 12, "+")
	expect_equal(p, "a+M_12+Y_14+S_17")
})


test_that("constructModifiedPeptide deals with vectors", {
	p = constructModifiedPeptide(c("_(ac)M(ox)CY(ph)RKS(ph)_", "_AGDS(ph)DS(ph)WK_"),
							 c(1, 1))
	expect_identical(p, c("a;M_1;Y_3;S_6", "S_4;S_6"))

	p = constructModifiedPeptide(c("_(ac)M(ox)CY(ph)RKS(ph)_", "_AGDS(ph)DS(ph)WK_"),
								 c(3, 10))
	expect_identical(p, c("a;M_3;Y_5;S_8", "S_13;S_15"))

	p = constructModifiedPeptide(c("_(ac)M(ox)CY(ph)RKS(ph)_", "_AGDS(ph)DS(ph)WK_"),
								 c(3, 10), "-")
	expect_identical(p, c("a-M_3-Y_5-S_8", "S_13-S_15"))
})


test_that("constructModifiedPeptide recyles positions", {
	p = constructModifiedPeptide(c("_(ac)M(ox)CY(ph)RKS(ph)_", "_AGDS(ph)DS(ph)WK_"),
								 1)
	expect_identical(p, c("a;M_1;Y_3;S_6", "S_4;S_6"))
	p = constructModifiedPeptide(c("_(ac)M(ox)CY(ph)RKS(ph)_", "_AGDS(ph)DS(ph)WK_"),
								 3)
	expect_identical(p, c("a;M_3;Y_5;S_8", "S_6;S_8"))
	p = expect_warning(constructModifiedPeptide(c("_(ac)M(ox)CY(ph)RKS(ph)_", "_AGDS(ph)DS(ph)WK_", "_AGDS(ph)DS(ph)WK_"),
								 c(1, 3)), regexp = "Recycling")
	expect_identical(p, c("a;M_1;Y_3;S_6", "S_6;S_8", "S_4;S_6"))
})

test_that("constructModifiedPeptide works with no position given", {
	p = constructModifiedPeptide(c("_(ac)M(ox)CY(ph)RKS(ph)_", "_AGDS(ph)DS(ph)WK_"))
	expect_identical(p, c("a;M_1;Y_3;S_6", "S_4;S_6"))
})

test_that("constructModifiedPeptide crashes with invalid positions", {
	expect_error(constructModifiedPeptide("_(ac)M(ox)CY(ph)RKS(ph)_", 1:2), "Incompatible sizes")
	expect_error(constructModifiedPeptide(c("_(ac)M(ox)CY(ph)RKS(ph)_", "_AGDS(ph)DS(ph)WK_"), 1:10), "Incompatible sizes")
})
