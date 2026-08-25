harrell_ref <- function(d) {
  survival::concordance(
    survival::Surv(d$time, d$status) ~ d$risk,
    reverse = TRUE
  )
}

pooled_from <- function(pc) (pc$S_ee + pc$S_ec) / (pc$W_ee + pc$W_ec)

test_that("event-censored ties are counted as comparable", {
  # A dies at 100, B censored at 100. B was alive when last seen, so we
  # know A's outcome was worse: the pair IS comparable.
  d <- make_tied_pair_data()
  pc <- pair_counts(d$time, d$status, d$risk, weights_harrell())
  expect_equal(pc$N_ee + pc$N_ec, 4)
})

test_that("event-event ties are NOT comparable", {
  # Two deaths at the same time: we cannot say who came first.
  d <- list(time = c(100, 100, 300), status = c(1, 1, 1), risk = c(0.9, 0.2, 0.1))
  pc <- pair_counts(d$time, d$status, d$risk, weights_harrell())
  expect_equal(pc$N_ee + pc$N_ec, 2)  # (A,C) and (B,C) only
})

test_that("pooled decomposition matches survival::concordance, continuous times", {
  d <- make_test_data(300, ties = FALSE, seed = 1)
  pc <- pair_counts(d$time, d$status, d$risk, weights_harrell())
  expect_equal(pooled_from(pc), harrell_ref(d)$concordance, tolerance = 1e-12)
})

test_that("pooled decomposition matches survival::concordance, TIED times", {
  # Regression: `time > time[i]` dropped 1,062 event-censored tied pairs and
  # biased the estimate by +0.0007 on this fixture.
  d <- make_test_data(300, ties = TRUE, seed = 1)
  pc <- pair_counts(d$time, d$status, d$risk, weights_harrell())
  expect_equal(pooled_from(pc), harrell_ref(d)$concordance, tolerance = 1e-12)
})

test_that("pooled decomposition matches survival::concordance with tied risks", {
  d <- make_test_data(300, ties = FALSE, seed = 2)
  d$risk <- as.numeric(cut(d$risk, 5))
  pc <- pair_counts(d$time, d$status, d$risk, weights_harrell())
  expect_equal(pooled_from(pc), harrell_ref(d)$concordance, tolerance = 1e-12)
})

test_that("pooled decomposition matches survival::concordance with both tied", {
  d <- make_test_data(300, ties = TRUE, seed = 3)
  d$risk <- as.numeric(cut(d$risk, 5))
  pc <- pair_counts(d$time, d$status, d$risk, weights_harrell())
  expect_equal(pooled_from(pc), harrell_ref(d)$concordance, tolerance = 1e-12)
})

test_that("comparable pair count matches survival's", {
  d <- make_test_data(300, ties = TRUE, seed = 1)
  pc <- pair_counts(d$time, d$status, d$risk, weights_harrell())
  cnt <- harrell_ref(d)$count
  expect_equal(
    pc$N_ee + pc$N_ec,
    unname(sum(cnt[c("concordant", "discordant", "tied.x")]))
  )
})

test_that("Uno matches survival timewt = n/G2 exactly when times are distinct", {
  d <- make_test_data(400, ties = FALSE, seed = 5)
  pc <- pair_counts(d$time, d$status, d$risk, weights_uno())
  ref <- survival::concordance(
    survival::Surv(d$time, d$status) ~ d$risk,
    reverse = TRUE, timewt = "n/G2"
  )$concordance
  expect_equal(pooled_from(pc), ref, tolerance = 1e-10)
})

test_that("Uno under TIED times agrees with survival only approximately", {
  # This package implements Uno's published PAIRWISE estimator: every pair
  # is weighted by 1 / Ghat(T_i)^2, as in Uno et al. (2011).
  # survival::concordance uses a counting-process form that applies its
  # weight per event TIME instead. The two coincide exactly when event
  # times are distinct (previous test, agreement at 1e-10) and diverge
  # once times are tied, because tie-clumping amplifies the difference
  # between per-pair and per-time weighting.
  #
  # The magnitude of that gap depends on the fixture's tie pattern, not
  # just on whether ties exist at all. `make_test_data(ties = TRUE)`
  # rounds to whole-unit day-scale times, matching real data. On this
  # fixture family (n = 300, censor_rate = 0.4, seeds 1-5) the *absolute*
  # gap against survival ranges from 8.8e-5 to 3.8e-4. `expect_equal()`'s
  # `tolerance` is relative-scaled (waldo::compare divides by the target's
  # magnitude, here ~0.7), so the number that actually governs this test
  # is the *relative* gap, which ranges from 1.3e-4 to 5.3e-4 over the
  # same seeds -- comfortably inside the 1e-3 tolerance, worst case
  # ~1.9x margin (seed 2). Real survival::lung data (185 unique times
  # among 227 ph.ecog-complete subjects, max tie group 3) measures an
  # absolute gap of 1.3e-4 / relative gap of 2.1e-4. (An earlier draft of
  # this fixture rounded the raw continuous scale directly, which piled
  # 300 subjects onto ~31 distinct times with ~10% at time zero -- a tie
  # pattern no real dataset exhibits -- and inflated the relative gap
  # past 1e-3 on most seeds. That was a fixture defect, not a property of
  # Uno's pairwise estimator, and was fixed rather than compensated for
  # with a looser tolerance.)
  #
  # Per-PAIR weighting is what makes the decomposition possible at all: a
  # per-event-time weight cannot be partitioned into event-event and
  # event-censored pair sets. The convention is therefore forced by the
  # package's purpose, not a defect.
  #
  # Harrell's weighting is unaffected and still matches exactly under
  # ties -- that is asserted separately above at 1e-12 (measured
  # agreement: 0.00e+00 on every tie pattern tested).
  d <- make_test_data(300, ties = TRUE, seed = 5)
  pc <- pair_counts(d$time, d$status, d$risk, weights_uno())
  ref <- survival::concordance(
    survival::Surv(d$time, d$status) ~ d$risk,
    reverse = TRUE, timewt = "n/G2"
  )$concordance
  expect_equal(pooled_from(pc), ref, tolerance = 1e-3)
})

test_that("the weighted identity holds for arbitrary custom weights", {
  d <- make_test_data(200, seed = 7)
  w <- weights_custom(function(t, G) 1 + t, name = "linear")
  pc <- pair_counts(d$time, d$status, d$risk, w)
  C_ee <- pc$S_ee / pc$W_ee
  C_ec <- pc$S_ec / pc$W_ec
  identity_rhs <- (pc$W_ee * C_ee + pc$W_ec * C_ec) / (pc$W_ee + pc$W_ec)
  expect_equal(pooled_from(pc), identity_rhs, tolerance = 1e-12)
})

test_that("negating risk gives 1 - C", {
  d <- make_test_data(200, seed = 9)
  a <- pooled_from(pair_counts(d$time, d$status,  d$risk, weights_harrell()))
  b <- pooled_from(pair_counts(d$time, d$status, -d$risk, weights_harrell()))
  expect_equal(a + b, 1, tolerance = 1e-12)
})

test_that("a negative weight is an error", {
  d <- make_test_data(50, seed = 11)
  w <- weights_custom(function(t, G) -1, name = "bad")
  expect_error(
    pair_counts(d$time, d$status, d$risk, w),
    "negative or non-finite"
  )
})

test_that("a non-finite weight is an error", {
  d <- make_test_data(50, seed = 11)
  w <- weights_custom(function(t, G) Inf, name = "bad")
  expect_error(
    pair_counts(d$time, d$status, d$risk, w),
    "negative or non-finite"
  )
})
