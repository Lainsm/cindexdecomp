new_cindex_weights <- function(name, fn) {
  structure(list(name = name, fn = fn), class = "cindex_weights")
}

check_tau <- function(tau) {
  if (!is.numeric(tau) || length(tau) != 1L) {
    stop("`tau` must be a single numeric value.", call. = FALSE)
  }
  if (!is.finite(tau) || tau <= 0) {
    stop("`tau` must be a positive, finite number.", call. = FALSE)
  }
  invisible(TRUE)
}

fmt_tau <- function(tau) format(tau, digits = 4, trim = TRUE)

#' Pair weightings for concordance decomposition
#'
#' @description
#' A concordance index of the form
#' `C = sum(w_ij * c_ij) / sum(w_ij)` is fully determined by its pair
#' weighting `w_ij`. These constructors supply that weighting to
#' [decompose_cindex()]. `weights_harrell()` gives Harrell's C,
#' `weights_uno()` gives Uno's inverse-probability-of-censoring-weighted
#' C, and `weights_custom()` accepts any user-supplied rule.
#'
#' @param tau Truncation time. For `weights_uno()`, pairs whose earlier
#'   event time exceeds `tau` receive weight zero; `NULL` (default) means
#'   no truncation. For `weights_truncated()`, `tau` is required.
#' @param fn A function of two arguments `(t, G)` returning a single
#'   non-negative finite number. `t` is the event time of the earlier
#'   member of the pair; `G` is a function returning the Kaplan-Meier
#'   censoring survival probability at a given time.
#' @param name A length-one character label used in printed output and
#'   plot annotations.
#'
#' @return An object of class `cindex_weights`.
#'
#' @examples
#' weights_harrell()
#' weights_uno()
#' weights_truncated(tau = 365)
#' weights_custom(function(t, G) 1 / G(t), name = "IPCW (power 1)")
#'
#' @name cindex_weights
NULL

#' @rdname cindex_weights
#' @export
weights_harrell <- function() {
  new_cindex_weights("Harrell", function(t, G) 1)
}

#' @rdname cindex_weights
#' @export
weights_uno <- function(tau = NULL) {
  if (!is.null(tau)) check_tau(tau)
  nm <- if (is.null(tau)) {
    "Uno (IPCW)"
  } else {
    paste0("Uno (IPCW, tau = ", fmt_tau(tau), ")")
  }
  new_cindex_weights(nm, function(t, G) {
    if (!is.null(tau) && t > tau) return(0)
    1 / G(t)^2
  })
}

#' @rdname cindex_weights
#' @export
weights_truncated <- function(tau) {
  check_tau(tau)
  new_cindex_weights(
    paste0("Truncated (tau = ", fmt_tau(tau), ")"),
    function(t, G) as.numeric(t <= tau)
  )
}

#' @rdname cindex_weights
#' @export
weights_custom <- function(fn, name = "Custom") {
  if (!is.function(fn)) {
    stop("`fn` must be a function of two arguments, (t, G).", call. = FALSE)
  }
  if (length(formals(fn)) < 2L) {
    stop("`fn` must accept two arguments: (t, G).", call. = FALSE)
  }
  if (!is.character(name) || length(name) != 1L) {
    stop("`name` must be a single character string.", call. = FALSE)
  }
  new_cindex_weights(name, fn)
}

#' @export
print.cindex_weights <- function(x, ...) {
  cat("<cindex_weights>", x$name, "\n")
  invisible(x)
}
