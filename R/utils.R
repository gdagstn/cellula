#' @noRd
.randomName <- function(){
  
  adjectives = c("rambunctious", "pertinacious", "holistic", "adamantine", 
                 "palyndromic", "hypertonic", "enantiomeric", "sybilline",
                 "brazen", "affine", "empyrean", "adjoining", "trustworthy",
                 "paralogous", "draconian", "obsequious", "nondescript", 
                 "diophantine", "hypogean", "preterintentional", "tawdry",
                 "auspicious", "emphatic", "sardonic", "continental")
  
  birds = c("parakeet", "greylag", "eagle", "stork", "hornbill", "pelican", 
            "kestrel", "flamingo", "magpie", "raven", "kingfisher", "owl", 
            "pigeon", "moorhen", "grebe", "mallard", "hummingbird", "nightingale",
            "gadwall", "goose", "swan", "sparrow", "finch", "vulture", "albatros")
  
  paste0(sample(adjectives, 1), "_", sample(birds, 1))
}

#' @noRd
.bluem <- function(txt) {
  left = "\033[38;5;27m"
  right = "\033[49m\033[39m"
  paste0(left, txt, right)
}

#' @noRd
.redm  <- function(txt) {
  left = "\033[38;5;1m"
  right = "\033[49m\033[39m"
  paste0(left, txt, right)
}


#' @noRd
.reorderNumericLevels <- function(f, rev = FALSE) {
  conv = suppressWarnings(as.numeric(levels(f)))
  if(!any(is.na(conv))) {
    fl = as.character(sort(as.numeric(levels(f))))
    if(rev) fl = fl[rev(seq_along(fl))]
    levels(f) = fl
    return(f)
  } else {
    return(f)
  }
}

#' Rescale
#'
#' Lifted from the \code{scales} package
#' @noRd

.rescalen <- function(x, to = c(0, 1), from = range(x, na.rm = TRUE, finite = TRUE)) {
  if (.zero_range(from) || .zero_range(to)) {
    return(ifelse(is.na(x), NA, mean(to)))
  }
  (x - from[1]) / diff(from) * diff(to) + to[1]
}

#' Zero range
#' 
#' Lifted from the \code{scales} package
#' @noRd

.zero_range <- function(x, tol = 1000 * .Machine$double.eps) {
  if (length(x) == 1) {
    return(TRUE)
  }
  if (length(x) != 2) 
    stop("x must be length 1 or 2")
  if (any(is.na(x))) {
    return(NA)
  }
  if (x[1] == x[2]) {
    return(TRUE)
  }
  if (all(is.infinite(x))) {
    return(FALSE)
  }
  m <- min(abs(x))
  if (m == 0) {
    return(FALSE)
  }
  abs((x[1] - x[2])/m) < tol
}

  
