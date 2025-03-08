#' @noRd
.randomName <- function(){
  
  adjectives <- c("rambunctious", "pertinacious", "holistic", "adamantine", 
                 "palyndromic", "hypertonic", "enantiomeric", "sybilline",
                 "brazen", "affine", "empyrean", "adjoining", "trustworthy",
                 "paralogous", "draconian", "obsequious", "nondescript", 
                 "diophantine", "hypogean", "preterintentional", "tawdry",
                 "auspicious", "emphatic", "sardonic", "continental", "sacrosanct",
				 "platitudinous", "shambolic", "sycophantic", "impervious")
  
  birds <- c("parakeet", "greylag", "eagle", "stork", "hornbill", "pelican", 
            "kestrel", "flamingo", "magpie", "raven", "kingfisher", "owl", 
            "pigeon", "moorhen", "grebe", "mallard", "hummingbird", "nightingale",
            "gadwall", "goose", "swan", "sparrow", "finch", "vulture", "albatros",
			"tit", "jay", "woodpecker", "quetzal")
  
  paste0(sample(adjectives, 1), "_", sample(birds, 1))
}

#' @noRd
.bluem <- function(txt) {
  left <- "\033[38;5;27m"
  right <- "\033[49m\033[39m"
  paste0(left, txt, right)
}

#' @noRd
.redm  <- function(txt) {
  left <- "\033[38;5;1m"
  right <- "\033[49m\033[39m"
  paste0(left, txt, right)
}


#' @noRd
.reorderNumericLevels <- function(f, rev = FALSE) {
  conv <- suppressWarnings(as.numeric(levels(f)))
  if(!any(is.na(conv))) {
    fl <- as.character(sort(as.numeric(levels(f))))
    if(rev) fl <- fl[rev(seq_along(fl))]
    levels(f) <- fl
    return(f)
  } else {
    return(f)
  }
}


#' Check functional dependencies
#' @noRd
checkFunctionDependencies <- function(depdf) {
  
  dependency_check = sapply(depdf$package, function(package) 
    !requireNamespace(package, quietly = TRUE))
  
  depdf$command = apply(depdf, 1, function(x) {
    if(x[2] == "BioC") {
      paste0("BiocManager::install(", "\"", x[1], "\")" )
    } else if(x[2] == "CRAN") {
      paste0("install.packages(\"", x[1], "\")" )
    } else if(grepl("github", x[2])) {
      repo = unlist(strsplit(x[2], split = ":"))[2]
      paste0("remotes::install_github(\"", repo, "/", x[1], "\")" )
    }
    })
  
  if(any(dependency_check)){
    message("Required packages not installed: ", .redm(sum(dependency_check)))
    message("Install the missing packages as follows:")
    invisible(sapply(depdf$command[dependency_check], function(x) {
      message(.bluem(paste0("`", x, "`")))
    }))
  }
  return(any(dependency_check))
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


#' Smooth polylines
#'
#' Uses Gaussian kernel smoothing to smooth polylines
#'
#' @param poly a data frame containing ordered coordinates with polygon vertices
#' @param smoothness numeric, the extent of kernel smoothing. Higher means
#'     rounder shapes. Default is 3.
#' @param min_points numeric, the minimum number of vertices to smooth.
#'     Default is 8.
#' @param n_dense numeric, the number of points to add to the polygon for more
#'     smoothing. Default is 10.
#'
#' @return A data frame containing ordered coordinates for polyline vertices.
#'
#' @details This is a refactoring of \code{\link{[smoothr](smooth_ksmooth)}} to isolate 
#'    the necessary code and avoid heavy GDAL-based dependencies. The code has 
#'    been simplified as well. Internal use only.
#'
#'
#' @author Matthew Strimas-Mackey, modified by Giuseppe D'Agostino
#'
#' @importFrom stats ksmooth
#'
#' @noRd

.smoothPolyline <- function(poly, smoothness = 3, min_points = 8, n_dense = 10) {
  if (nrow(poly) < min_points) {
    poly_sm <- poly
  } else {
    poly_coords <- as.matrix(poly[,c(1,2)])
    d_poly <- sqrt(rowSums(diff(as.matrix(poly_coords))^2))
    bandwidth <- mean(d_poly) * smoothness
    dense_poly <- .addPoints(poly_coords, steps = n_dense)
    npt <- nrow(dense_poly)
    d_dense <- sqrt(rowSums(diff(as.matrix(dense_poly))^2))
    d_x <- c(0, cumsum(d_dense))
    poly_sm <- NULL
    for (i in seq_len(ncol(dense_poly))) {
      ks <- ksmooth(d_x, dense_poly[, i], n.points = length(d_x),
                    kernel = "normal", bandwidth = bandwidth)
      poly_sm <- cbind(poly_sm, ks[["y"]])
    }
    colnames(poly_sm) <- c("x", "y")
    poly_sm <- as.data.frame(poly_sm)
  }
  return(poly_sm)
}

#' Add points
#'
#' Increases the number of points in a polyline while maintaining the shape
#'
#' @param poly a data frame containing ordered coordinates with polyline vertices
#' @param steps numeric, the number of points that should be added between
#'    each point. Default is 5
#'
#' @return A data frame containing densified coordinates for polyline vertices 
#'
#' @details Internal use only.
#'
#' @author Giuseppe D'Agostino
#' 
#' @noRd

.addPoints <- function(poly, steps = 5) {
  colnames(poly) <- NULL
  polygon_coords <- as.matrix(poly[,c(1,2)])
  new_coords <- poly[1,]
  for(i in seq_len(nrow(poly)-1)){
    new_xy <- cbind(seq(polygon_coords[i,1],
                       polygon_coords[i+1,1], length.out = steps+1),
                   seq(polygon_coords[i,2],
                       polygon_coords[i+1,2], length.out = steps+1))
    tmp_coords <- polygon_coords[i,]
    new_coords <- rbind(new_coords, tmp_coords, new_xy)
  }
  new_coords <- new_coords[2:nrow(new_coords),]
  new_coords <- new_coords[!duplicated(new_coords),]
  rownames(new_coords) <- NULL
  colnames(new_coords) <- c("x", "y")
  new_coords
}