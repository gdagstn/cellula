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

#' Qualitative color palette
#' 
#' Preset qualitative color palettes from qualpalr
#' 
#' @param n numeric, the number of colors. Must be between 2 and 30 inclusive.
#' @returns a color palette with n colors 
#' @details
#'     These palettes were generated using the \code{qualpalr} package by 
#'     Johan Larsson, and are optimized for colors that are maximally distant
#'     in a perceptual color space. They are listed here to avoid a dependency
#'     but it is important to credit the original author of the code that 
#'     generated these palettes.

.cpal_qual = function(n){
  
  if(n < 2) stop("n must be at least 2")
  if(n > 30) stop("n cannot be greater than 30")
  
  list(c("#9FCB6A", "#8367CC"),
       c("#C96C69", "#6E9DCE", "#9FCB6A"),
       c("#6DCC76", "#CA69C2", "#89ACCE", "#CD9F7C"),
       c("#73CA6F", "#D37DAD", "#C6DBE8", "#6C7DCC", "#D0A373"),
       c("#73CA6F", "#706EC4", "#DDCAE6", "#CFAB7B", "#C96D94", "#6CB4C7"),
       c("#BFCD6D", "#7669CA", "#E2C9E7", "#C5916F", "#C5EAE0", "#C66FA0", 
         "#7CA5C6"),
       c("#B66CCC", "#6DCC76", "#CB9B6A", "#E2C9E7", "#7590CF", "#76C3CA", 
         "#ECEAC5", "#CA7488"),
       c("#BFCD6D", "#7669CA", "#C96C69", "#6CB4C7", "#B4B2E5", "#C1E6D0", 
         "#CA69C2", "#E2B3C2", "#D0A373"),
       c("#C3D37B", "#8367CC", "#C96C69", "#B8DAE7", "#D5BCE5", "#7395C3", 
         "#CB70B2", "#6ECEA4", "#C89E6D", "#E3D6CB"),
       c("#A8C56F", "#7669CA", "#CF707B", "#C9DDE4", "#CB9B6A", "#E8BACA", 
         "#CA69C2", "#7B9DBA", "#71CBA7", "#E8DAC1", "#B9A6DF"),
       c("#CA69C2", "#6DCC76", "#C4D2EB", "#CB9B6A", "#7669CA", "#CF707B", 
         "#E5C9C6", "#6FC9C1", "#7395C3", "#CACA6E", "#DCB2E4", "#D0DFCA"),
       c("#BFCD6D", "#A76DC9", "#CD7F6D", "#73A6C8", "#6CC884", "#E6BEBE", 
         "#9ED7D3", "#6C7DCC", "#C9CDE8", "#DFA9E0", "#DBDCC8", "#B4A282", 
         "#C96D94"),
       c("#9FCB6A", "#8367CC", "#C96C69", "#75B1C6", "#DCC5C2", "#C89E6D", 
         "#CA69C2", "#C9CDE8", "#7ABB9B", "#BCE8E7", "#7590CF", "#E6E5BE", 
         "#BC9EDA", "#DD9BB9"),
       c("#CACA6E", "#8367CC", "#C96C69", "#6FC9C1", "#DDCDE3", "#CB9B6A", 
         "#75BD75", "#E6CFBA", "#C66FA0", "#C39BDE", "#DFA9B5", "#CADEE2", 
         "#7089C6", "#82A3B2", "#D1E9BF"),
       c("#73CA6F", "#A76DC9", "#CB9B6A", "#73A6C8", "#DDCDE3", "#C96C69", 
         "#CACA6E", "#6C7DCC", "#D3E8C7", "#70C7B2", "#E9D1C4", "#C49CA4", 
         "#B8DAE7", "#C66FA0", "#AAA9E1", "#D9A2DA"),
       c("#92CC71", "#B66CCC", "#6CB4C7", "#C96C69", "#E3D0C4", "#C89E6D", 
         "#706EC4", "#BDC7DE", "#E7C6DE", "#D3E8C7", "#B39FD8", "#CA80A4", 
         "#7ABB9B", "#CCC681", "#7395C3", "#B2E4E2", "#DDA29F"),
       c("#73CA6F", "#B66CCC", "#CB9B6A", "#75B1C6", "#E7C9D0", "#C96C69", 
         "#CACA6E", "#7590CF", "#CADEE2", "#75CFB8", "#D1E9BF", "#7669CA", 
         "#B3BEE2", "#E3CEB1", "#A894BC", "#CC77A5", "#E6B7E6", "#D49A9F"),
       c("#B0CE74", "#7669CA", "#C96C69", "#E8C7EA", "#6FC9C1", "#E8CBC4", 
         "#C66FA0", "#CADEE2", "#82A3B2", "#A894BC", "#ECEAC5", "#D1B87A", 
         "#7590CF", "#B99482", "#DFA2B2", "#B5BEE1", "#72C17E", "#B66CCC", 
         "#BCE9BD"),
       c("#73CA6F", "#AF75CE", "#C5916F", "#ABC0D4", "#ECEAC5", "#C49CA4", 
         "#C66FA0", "#706EC4", "#DDCDE3", "#7395C3", "#ACA1DC", "#E5C8BB", 
         "#C96C69", "#73C8CD", "#DD9FD7", "#C1AE76", "#C4E6E3", "#C3D37B", 
         "#BCE9BD", "#7ABB9B"),
       c("#6DCC76", "#B66CCC", "#C5916F", "#82A3B2", "#ECEAC5", "#DDCAD6", 
         "#C9DDE4", "#7089C6", "#CEBE70", "#B9E7C8", "#DD9FD7", "#DBA0A5", 
         "#C96C69", "#E0C2AC", "#AB98D9", "#B3BEE2", "#7669CA", "#7BCCD2", 
         "#C66FA0", "#B6D482", "#8CB596"),
       c("#6DCC76", "#CA69C2", "#6E9DCE", "#C5916F", "#DDCAE6", "#E5DACA", 
         "#CADEE2", "#ACA1DC", "#D1B87A", "#71BDD0", "#C96C69", "#77C0A6", 
         "#6C7DCC", "#D6E0A4", "#D9A2DA", "#9A6ECD", "#BBE4C4", "#CA979A", 
         "#AEC1DF", "#EAC7CB", "#97B284", "#C96D94"),
       c("#CA69C2", "#73CA6F", "#A5C1E1", "#C5916F", "#ECEAC5", "#6C7DCC", 
         "#DDCDE3", "#C49CA4", "#B9E7C8", "#CACA6E", "#7BCCD2", "#E8CBC4", 
         "#C96C69", "#A78ABA", "#CADEE2", "#7ABB9B", "#C0AF8B", "#DFA9E0", 
         "#C96D94", "#AAA9E1", "#6E9DCE", "#B6C69E", "#8367CC"),
       c("#CA69C2", "#6DCC76", "#6E9DCE", "#C5916F", "#E2C9E7", "#B2E4E2", 
         "#E8DEC7", "#8367CC", "#C3D37B", "#C96C69", "#C1E7BC", "#77C0A6", 
         "#71BDD0", "#BAC8D9", "#CA979A", "#D1B87A", "#D89FC7", "#8385B6", 
         "#B4B2E5", "#B78DD8", "#97B284", "#C96D94", "#EAC7CB", "#C7AF9D"),
       c("#B66CCC", "#73CA6F", "#DFB3A6", "#75B1C6", "#DDCAE6", "#C8E6D6", 
         "#C96C69", "#8589CF", "#E2D9A7", "#D393A4", "#B4A282", "#7ABB9B", 
         "#A9B783", "#D699D0", "#C6DBE8", "#E3D6CB", "#B9A6DF", "#B3BEE2", 
         "#8DD4D4", "#8367CC", "#BBE0A7", "#779DCC", "#C66FA0", "#E8C1CC", 
         "#BE8A76"),
       c("#CA69C2", "#73CA6F", "#6E9DCE", "#CB9B6A", "#BBE1E3", "#EAC7CB", 
         "#8367CC", "#EAE4C7", "#CEBE70", "#D1C9DE", "#C77694", "#75CFB8", 
         "#CD7F6D", "#6CB4C7", "#9A8ABE", "#C1DEA1", "#C49CA4", "#D1BBA7", 
         "#E2B6D9", "#ABC0D4", "#6C7DCC", "#ABAFDD", "#B9E7C8", "#8CB596", 
         "#C894D9", "#A9B783"),
       c("#CA69C2", "#73CA6F", "#7CA5C6", "#C5916F", "#DDCAD6", "#ECEAC5", 
         "#7669CA", "#C3E7E3", "#BFCD6D", "#C96C69", "#E5C1B4", "#75CFB8", 
         "#BCE9BD", "#ABAFDD", "#A789CE", "#C96D94", "#C6D2E2", "#D5B3E2", 
         "#7089C6", "#8BC4D0", "#DDAABC", "#DAC599", "#8CB596", "#BD8BB0", 
         "#D09491", "#B6C69E", "#B4A282"),
       c("#6DCC76", "#B66CCC", "#C5916F", "#A8C7DE", "#E3D6CB", "#6C7DCC", 
         "#DABAE8", "#CACA6E", "#C96D94", "#B9E7C8", "#C96C69", "#C0AF8B", 
         "#73C8CD", "#E8BACA", "#97B284", "#AAA9E1", "#DA95CB", "#DFB3A6", 
         "#D5929A", "#8367CC", "#7395C3", "#A78ABA", "#71CBA7", "#E2E8C4", 
         "#C2E8E9", "#CCCAE1", "#82A3B2", "#B1DA99"),
       c("#9FCB6A", "#8367CC", "#C96C69", "#B8DAE7", "#EBC4D6", "#7395C3", 
         "#D0A373", "#ECEAC5", "#CA69C2", "#8CBEA9", "#BFB993", "#CFB5E5", 
         "#C96D94", "#D998C8", "#D9CBC1", "#B187D0", "#C2C9E7", "#C6E9DA", 
         "#D5929A", "#6CC884", "#6C7DCC", "#76C3CA", "#BDE0AE", "#A09DC4", 
         "#DFB3A6", "#BE8A76", "#97B284", "#82A3B2", "#CACA6E"),
       c("#6DCC76", "#CA69C2", "#6E9DCE", "#C5916F", "#BCE8E7", "#E5CAD4", 
         "#DFE3B4", "#8367CC", "#C96C69", "#DCB2E4", "#A9B783", "#6CB4C7", 
         "#AC8ED6", "#C96D94", "#CFC9E4", "#DEBBAD", "#7C86D0", "#8CB596", 
         "#A4ACDE", "#AEC5D0", "#81CBBE", "#BD8BB0", "#B4A282", "#D09491", 
         "#C0D8C3", "#CACA6E", "#E5DACA", "#DAC599", "#B1DA99", "#DDAABC"
       ))[[n-1]]
}

#' Qualitative color palette - protanopy
#' 
#' Preset qualitative color palettes from qualpalr
#' 
#' @param n numeric, the number of colors. Must be between 2 and 30 inclusive.
#' @returns a color palette with n colors, optimized for protanopy  
#' @details
#'     These palettes were generated using the \code{qualpalr} package by 
#'     Johan Larsson, and are optimized for colors that are maximally distant
#'     in a perceptual color space. They are listed here to avoid a dependency
#'     but it is important to credit the original author of the code that 
#'     generated these palettes.
#'     
#'     These palettes simulate protanopy (color vision deficiency) at a severity
#'     of 0.5.

.cpal_qual_protan <- function(n){
  
  if(n < 2) stop("n must be at least 2")
  if(n > 30) stop("n cannot be greater than 30")
  
  list(
  c("#D4C66B", "#6670CF"),
  c("#D4E4E7", "#6170CC", "#977568"),
  c("#D4C66B", "#6170CC", "#CDDCE6", "#9A797B"),
  c("#D7C56C", "#6170CC", "#977568", "#C6BEE7", "#A9C5B9"),
  c("#D1BB6E", "#6170CC", "#977568", "#98B7D0", "#D8E4D5", "#A38EAC"
  ),
  c("#D7C56C", "#6170CC", "#C7B6BC", "#A5856C", "#A7C1A6", "#91ABD2", 
    "#927894"),
  c("#D7C56C", "#6170CC", "#D7E4E9", "#A5856C", "#91ABD2", "#927894", 
    "#A2B39B", "#C0AFB4"),
  c("#6170CC", "#ADBE74", "#977568", "#BDC5EB", "#C6B4B1", "#97859D", 
    "#A9C5B9", "#F0E8C5", "#B89C69"),
  c("#EBE7C1", "#6170CC", "#A5856C", "#D3CDE9", "#A9C5B9", "#A6B274", 
    "#978791", "#CABCB7", "#A29DDA", "#90AFC8"),
  c("#F0E8C5", "#6170CC", "#9E8071", "#CCC7EA", "#C5AF70", "#C4B7B2", 
    "#977F94", "#A8C3A3", "#D4E4E7", "#92A1B2", "#A094D0"),
  c("#ADBE74", "#6170CC", "#977568", "#CCC4E0", "#927894", "#F0E8C5", 
    "#BFAEAA", "#B89C69", "#A9C1AD", "#D7E4E9", "#929FDC", "#99AFBE"
  ),
  c("#6170CC", "#E9E0AC", "#977568", "#D2CADE", "#A2B7A6", "#B89C69", 
    "#907AA0", "#A4A3DF", "#98B7D0", "#CFC0B6", "#DBE6E0", "#AE9AA4", 
    "#AFBC6D"),
  c("#C4C168", "#6170CC", "#D2CADE", "#977568", "#A1C0C1", "#907AA0", 
    "#CFC0B6", "#A4A3DF", "#A8BE99", "#B49973", "#AE9AA4", "#EBE7C1", 
    "#DBE6E0", "#92A1B2"),
  c("#D7C56C", "#6170CC", "#D7E4E9", "#977568", "#B4A4C2", "#CCC2B9", 
    "#B19B77", "#A2BDB2", "#A7BB83", "#927894", "#EBE7C1", "#8C9BDD", 
    "#AF999A", "#92A1B2", "#BCC6E8"),
  c("#6170CC", "#E9E0AC", "#9A797B", "#CFDAE8", "#A2B39B", "#92A1B2", 
    "#907AA0", "#C2AB79", "#D5C7BD", "#8C9BDD", "#B9AAB5", "#ADBE74", 
    "#C6BEE7", "#D8E4D5", "#A9998B", "#A1C0C1"),
  c("#DFD083", "#6170CC", "#C8D4E9", "#9A797B", "#A2B39B", "#907AA0", 
    "#92A1B2", "#D0C4C3", "#8C9BDD", "#DBE6E0", "#C2AF95", "#A1C0C1", 
    "#BAB0D3", "#EBE3C6", "#AE9AA4", "#A5856C", "#ADBE74"),
  c("#6670CF", "#DBCA75", "#CFDAE5", "#9A797B", "#B4A8D8", "#A7AF96", 
    "#907AA0", "#98B7D0", "#CFC0B6", "#7F95C5", "#EBE3BD", "#BBA575", 
    "#AE9AA4", "#A3B9B6", "#CCC0D0", "#D2DDCE", "#A5856C", "#ADBE74"
  ),
  c("#DBCA75", "#6670CF", "#CFDAE8", "#977568", "#A8C3A3", "#977F94", 
    "#92A1B2", "#CCBCA3", "#AE9C9A", "#A292C7", "#AE946E", "#EBE7C1", 
    "#A1C0C1", "#BAACC0", "#D6CBCB", "#AFB9E5", "#A6B274", "#D2E0D4", 
    "#7B92D0"),
  c("#D7C56C", "#6170CC", "#CCD8E8", "#997D88", "#A2B39B", "#819CCF", 
    "#DAE1D1", "#B6A3B9", "#B89C69", "#99AFBE", "#C8BCA7", "#CFC4C9", 
    "#C3C0E9", "#907AA0", "#9E8071", "#B4A09E", "#E8E0B3", "#A6B274", 
    "#9E95D7", "#B2CBC5"),
  c("#C4C168", "#6670CF", "#B0C8D8", "#977568", "#D6E1C8", "#B3A5BB", 
    "#CCBCA3", "#907AA0", "#A0AEE1", "#A2BDB2", "#D2C5C7", "#A6BE90", 
    "#B89C69", "#D1CCEB", "#B09990", "#9F838F", "#92A1B2", "#DAE4E5", 
    "#E9E0AC", "#A094D0", "#758BC8"),
  c("#DBCA75", "#6670CF", "#CFDAE5", "#9A797B", "#A8BE99", "#92A4C5", 
    "#E1DAC9", "#B6A3B9", "#B89F6C", "#927894", "#B4A09E", "#C9C8EC", 
    "#ABA2DD", "#CFC4C9", "#B4C6BE", "#A5856C", "#C1B89D", "#9EBCCB", 
    "#D4D9A6", "#AFBC6D", "#8F7EBE", "#758BC8"),
  c("#D4C66B", "#6170CC", "#ACC5CC", "#997D88", "#BAB1E1", "#BDA695", 
    "#AFC69E", "#819CCF", "#D6CCC6", "#907AA0", "#E2D5AF", "#9E8071", 
    "#D2C9DB", "#D8E2CE", "#9B8DD4", "#B89C69", "#BAA9AE", "#92A1B2", 
    "#A3B5AA", "#B0C1E2", "#A899B3", "#D7E4E9", "#ACAD83"),
  c("#D4C66B", "#6670CF", "#B8CBDD", "#9A797B", "#D6E1C8", "#B4A4C2", 
    "#BEA79F", "#B89F6C", "#927894", "#92A1B2", "#A6B274", "#A2ACE2", 
    "#A3B9B6", "#9584C7", "#A5856C", "#D7E4E3", "#AABEA1", "#CDBFCA", 
    "#E9E0AC", "#A8969E", "#CCC7EA", "#D7CEC4", "#CCBD9C", "#758BC8"
  ),
  c("#6670CF", "#D7C56C", "#CFDAE8", "#9A797B", "#A8C3A3", "#9D92B1", 
    "#C9B6A5", "#A6B274", "#A5856C", "#AFBCE3", "#CFC4C9", "#D4D9A6", 
    "#92A1B2", "#BFA571", "#B4A0AB", "#E1DAC9", "#ACB8B3", "#A0C0CE", 
    "#758BC8", "#DBE6E0", "#8F7EBE", "#927894", "#BFB4CE", "#A4A3DF", 
    "#B0968F"),
  c("#D7C56C", "#6170CC", "#C8D4E9", "#977568", "#A997BB", "#C0D3BC", 
    "#BDABA0", "#92A1B2", "#F0E8C5", "#A1C0C1", "#A6B274", "#A0AEE1", 
    "#BFA571", "#927894", "#D7E4E3", "#B4A0AB", "#D5CDD4", "#7C8DC6", 
    "#8F7EBE", "#D3C3A2", "#988286", "#C5BCD9", "#A58E75", "#A7AF96", 
    "#DDD6CB", "#CFD49C"),
  c("#DBCA75", "#6170CC", "#B6CCDD", "#977568", "#A7C1A6", "#C2AEAC", 
    "#ABA2DD", "#907AA0", "#AE946E", "#92A1B2", "#DED2C3", "#997D88", 
    "#D1CCEB", "#7B92D0", "#D7E4E3", "#AFBC6D", "#A5B5B5", "#E8E0B3", 
    "#D5CDD4", "#BCADC7", "#D6E1C8", "#A7B5DD", "#8E80CE", "#C9B59B", 
    "#BEC39C", "#AA8F88", "#9D919F"),
  c("#6670CF", "#D7C56C", "#A6C4D3", "#977568", "#A895AC", "#D7CEC4", 
    "#97AAD9", "#A7AF96", "#BBA897", "#907CB3", "#F0E8C5", "#D4E4E7", 
    "#D1CCD6", "#997D88", "#A2BDB2", "#D0DDC7", "#BFA571", "#758BC8", 
    "#92A1B2", "#A59CD6", "#A6B274", "#BFB4CE", "#BCC6E8", "#AF999A", 
    "#D3C3A2", "#CFD49C", "#C4B5B8", "#A48A70"),
  c("#DFD083", "#6170CC", "#977568", "#D1CBE1", "#A9C5B9", "#BBA897", 
    "#907AA0", "#92A4C5", "#D6CCC6", "#997D88", "#ADBE74", "#EBE7C1", 
    "#A7AF96", "#B89C69", "#D8E4D5", "#B5B2E3", "#C0B2BB", "#9EB2B7", 
    "#CBD3D7", "#758BC8", "#BFD0A9", "#A895AC", "#B2C5DF", "#AF999A", 
    "#A292C7", "#8977C3", "#A48A70", "#D6C8A9", "#BEB479"),
  c("#D4C66B", "#6170CC", "#D2CEE0", "#977568", "#A2B39B", "#927894", 
    "#ABA2DD", "#BAA6A5", "#92A1B2", "#F0E8C5", "#B2987D", "#DAE4E5", 
    "#A5B5B5", "#819CCF", "#C5D3A9", "#ADB8DC", "#B6CCDD", "#AEBC76", 
    "#D2C5C7", "#C2D2C3", "#9D92B1", "#BAAABA", "#BEAC75", "#BCB09D", 
    "#907CB3", "#D6C8A0", "#DBD1C3", "#9D8B96", "#8488D1", "#9E8A83"
  ))[[n-1]]
}

#' Qualitative color palette - tritanopy
#' 
#' Preset qualitative color palettes from qualpalr
#' 
#' @param n numeric, the number of colors. Must be between 2 and 30 inclusive.
#' @returns a color palette with n colors, optimized for tritanopy  
#' @details
#'     These palettes were generated using the \code{qualpalr} package by 
#'     Johan Larsson, and are optimized for colors that are maximally distant
#'     in a perceptual color space. They are listed here to avoid a dependency
#'     but it is important to credit the original author of the code that 
#'     generated these palettes.
#'     
#'     These palettes simulate tritanopy (color vision deficiency) at a severity
#'     of 0.5.

.cpal_qual_tritan <- function(n){
  
  if(n < 2) stop("n must be at least 2")
  if(n > 30) stop("n cannot be greater than 30")
  
  list(
    c("#6FC98B", "#C76DAC"),
    c("#A6CB90", "#CB6B6A", "#6881B8"),
    c("#6FC98B", "#CB6B6A", "#716DB1", "#DCCEDE"),
    c("#71C9CB", "#CB6B6A", "#716DB1", "#CEC686", "#E7CAE1"),
    c("#6CBFC3", "#C76DAC", "#C8D39F", "#D5C7E0", "#6F77B6", "#D0947D"
    ),
    c("#75CFBE", "#CB6B6A", "#716DB1", "#C3C885", "#B0B9D8", "#E9C3BF", 
      "#CA8DC3"),
    c("#6FC98B", "#C973A2", "#7098B7", "#CADFE1", "#D2CB8D", "#7F6CB3", 
      "#E9C5D5", "#C88F78"),
    c("#ABC185", "#A471B2", "#76B0C7", "#CB6B6A", "#C6E9DE", "#D2A887", 
      "#E0A6BF", "#CBC8DE", "#7483B5"),
    c("#95C888", "#C76DAC", "#7098B7", "#B6A08B", "#7F6CB3", "#E3B7CC", 
      "#EDE8CF", "#CB6B6A", "#6DCAC2", "#C4CFDE"),
    c("#6DCAC2", "#D07079", "#716DB1", "#EDE8CF", "#EBC4D2", "#C9A282", 
      "#95C888", "#B878B5", "#6CA0C1", "#C5D3DE", "#A9A4CD"),
    c("#79CCD0", "#CB6B6A", "#716DB1", "#C8C388", "#C0B1D6", "#7098B7", 
      "#E2CCC7", "#77BA87", "#B370B4", "#CFE9CF", "#CB9B7A", "#DD9CB2"
    ),
    c("#6BC89F", "#C57194", "#6881B8", "#CFAD81", "#BFC5E0", "#EDE8CF", 
      "#73B2C0", "#9773B5", "#C2E7E4", "#ABC185", "#D0A5D0", "#E7BEBE", 
      "#C48477"),
    c("#C76DAC", "#6FC98B", "#7098B7", "#C9A282", "#E1CAE0", "#CB6B6A", 
      "#E9DDCD", "#B5BE86", "#6C72AE", "#6CBFC3", "#C0E6D5", "#AC94BD", 
      "#DFA2AE", "#BBCCDF"),
    c("#79CFAB", "#C76DAC", "#EAD0C7", "#91A4CB", "#CB6B6A", "#C4D7E4", 
      "#6C72AE", "#6AB5C2", "#DDA3B5", "#D9C6E1", "#CE9877", "#D8E5CF", 
      "#D1BF8B", "#A3C683", "#A68CAE"),
    c("#6FC98B", "#B370B4", "#C3D4E5", "#C59F91", "#6D8CB7", "#E6D9CE", 
      "#CB6B6A", "#A99BC8", "#E1CAE0", "#C2E3D1", "#D497B1", "#D2BB84", 
      "#6DCAC2", "#716DB1", "#B9D097", "#81A3AE"),
    c("#76C685", "#C973A2", "#7293BF", "#D4B58A", "#D4C7DB", "#6AB5C2", 
      "#EDE8CF", "#CB6B6A", "#716DB1", "#BBE3CC", "#A99BC8", "#ABB590", 
      "#BED6E3", "#E9CBC6", "#BD958F", "#DEA7C4", "#A471B2"),
    c("#B370B4", "#7DCF99", "#CFAD81", "#6881B8", "#C4D7E4", "#CB6B6A", 
      "#81A3AE", "#E6DCC6", "#A8ABD2", "#DDCAD3", "#D0A5D0", "#7BD0CC", 
      "#D6A6A3", "#C3C885", "#CAEACF", "#D37EA1", "#7F6CB3", "#98B090"
    ),
    c("#6DC595", "#C57194", "#8793C4", "#DECCA6", "#E7C9E1", "#BBE9E7", 
      "#A471B2", "#CF7E72", "#ABC185", "#6CBFC3", "#B4C5D8", "#D6E5CA", 
      "#C49CA2", "#D49BCB", "#B6A08B", "#799FB3", "#E9CBC6", "#716DB1", 
      "#B4A7C2"),
    c("#7DCF99", "#B370B4", "#D0947D", "#BABCD3", "#6881B8", "#C3C885", 
      "#E3D5CE", "#7AA7BE", "#D2E7CA", "#C1DDE2", "#6DCAC2", "#A395C1", 
      "#E7CAE1", "#7F6CB3", "#DFA9B2", "#CB799A", "#D4B58A", "#98B090", 
      "#CB6B6A", "#D49BCB"),
    c("#73BF8F", "#C76DAC", "#6881B8", "#E8CEC0", "#CADFE1", "#CB6B6A", 
      "#D0A5D0", "#71C9CB", "#9773B5", "#DCCEDE", "#D9A092", "#D1BF8B", 
      "#6CA0C1", "#B88C9A", "#A09CC7", "#B1C0D9", "#E4E6CD", "#A9BF8B", 
      "#B5E1C1", "#E8BBC6", "#B6A08B"),
    c("#6DCAC2", "#CB6B6A", "#7F6CB3", "#E5C7DA", "#D2BB84", "#99A7CB", 
      "#76C685", "#BDE7C8", "#CE9877", "#E3D5CE", "#6881B8", "#AC7DB1", 
      "#ABB590", "#DB9DB8", "#76B0C7", "#C8CEE1", "#BFA8D2", "#B8E0E4", 
      "#DDB4B0", "#C57194", "#D9DDB3", "#8DB49E"),
    c("#C76DAC", "#76C685", "#749FC0", "#E6D9CE", "#C88F78", "#7F6CB3", 
      "#CBE6E5", "#DEABD3", "#E7B7B8", "#DCCEDE", "#B88C9A", "#6FBECA", 
      "#D07079", "#C2AD94", "#ADC3D7", "#9EA0C4", "#AF8BBE", "#6881B8", 
      "#ABC185", "#CBE1C3", "#DCD0A0", "#97DAC3", "#8DB49E"),
    c("#C76DAC", "#76C685", "#7293BF", "#E4C9C0", "#8FCDD5", "#CB6B6A", 
      "#C3C885", "#E1B1CF", "#9773B5", "#B1C0D9", "#C3AB85", "#EDE8CF", 
      "#6C72AE", "#BCDDB5", "#81A3AE", "#DCCEDE", "#CADED9", "#A09CC7", 
      "#C88F78", "#98B090", "#B6828E", "#E1A8AC", "#BE8EC2", "#89D4B8"
    ),
    c("#B370B4", "#87D299", "#6CA0C1", "#CFAD81", "#CB6B6A", "#E9C5D5", 
      "#C8DDE2", "#E9DDCD", "#74C4C8", "#7F6CB3", "#9EA0C4", "#6881B8", 
      "#C57194", "#ABC185", "#C19FC8", "#E0B2A9", "#DD9CB2", "#CAC1D8", 
      "#AFE3CE", "#8DB49E", "#C88F78", "#D0DFBD", "#A4BDD2", "#DCD0A0", 
      "#C3B7A4"),
    c("#C76DAC", "#76C685", "#7293BF", "#CB9B7A", "#DCCEDE", "#82C8CE", 
      "#E4E6CD", "#716DB1", "#BB98C6", "#B88C9A", "#C7B886", "#B4CADE", 
      "#B5E1C1", "#E4D0C7", "#CF7E72", "#DEABA8", "#98B090", "#C2E8E8", 
      "#81A3AE", "#B9D097", "#E1B1CF", "#B8B2D5", "#7BBAA3", "#918FBE", 
      "#A471B2", "#C96F8A"),
    c("#C76DAC", "#76C685", "#7293BF", "#EAD0C7", "#C2E0E2", "#BA827D", 
      "#C19ECD", "#D1BF8B", "#D0CAD9", "#6C72AE", "#B9D097", "#81A3AE", 
      "#E7BFD7", "#EDE8CF", "#E1A8AC", "#82C8CE", "#CAEACF", "#98B090", 
      "#ADC3D7", "#CE9877", "#89D4B8", "#A9A9C5", "#BC8DA7", "#A471B2", 
      "#968BC1", "#C96F8A", "#C5AD9B"),
    c("#6FC98B", "#B370B4", "#CE9877", "#6CA0C1", "#DCCBDF", "#C2E8E8", 
      "#D4D7A7", "#9694C6", "#CB6B6A", "#6DCAC2", "#7F6CB3", "#C57194", 
      "#E9D4C1", "#E8BBC6", "#D694C2", "#C4CFDE", "#8DB49E", "#D69299", 
      "#DAE1D1", "#D2BB84", "#C8AEA2", "#C0A8CD", "#ABB590", "#96C4D3", 
      "#6881B8", "#9CAAC0", "#B4D7A9", "#A1DDC3"),
    c("#6FCDAE", "#CB6B6A", "#6881B8", "#D9CAC4", "#B370B4", "#A0B6CA", 
      "#A3C683", "#C6AFD3", "#E1A8AC", "#B6A08B", "#7F6CB3", "#C57194", 
      "#D2CB8D", "#6AB5C2", "#B5E5D9", "#D694C2", "#D0947D", "#989ECA", 
      "#EDE8CF", "#D0CAD9", "#E9C5D5", "#9F89B4", "#7098B7", "#C4DABC", 
      "#C5DBE4", "#8FB997", "#96BCB5", "#B88C9A", "#DDBD9D"),
    c("#76C685", "#B370B4", "#7098B7", "#C88F78", "#E4C2DF", "#E3E0BB", 
      "#B2E4E2", "#6C72AE", "#A796B2", "#CCD3DF", "#B88C9A", "#6AB5C2", 
      "#CB6B6A", "#98B090", "#E3D5CE", "#D49BCB", "#B2B4D8", "#B6A08B", 
      "#D1DECF", "#9CBCD2", "#D2B6B9", "#8A91BD", "#B4D7A9", "#DCC2A4", 
      "#C3C885", "#937CB5", "#6FCDAE", "#96BCB5", "#C57194", "#DEA2A0"
    )
  )[[n-1]]
}


#' Sequential color palette - parula
#' 
#' Preset continuous color palettes from pals
#' 
#' @details
#'     These palettes were generated using the \code{pals} package by 
#'     Kevin Wright.

.cpal_seq_parula <- function() {
  c("#352A87", "#3439A8", "#214DC8", "#0E5FDB", "#056EDE", "#0F79D9", 
    "#1283D4", "#0D8FD1", "#089BCE", "#06A5C7", "#0BACBC", "#1CB1AE", 
    "#33B7A0", "#4EBB91", "#6DBE81", "#8ABE75", "#A3BD6A", "#BBBC60", 
    "#D1BA58", "#E6B94E", "#F9BD3F", "#FBC831", "#F8D626", "#F5E71A", 
    "#F9FB0E")
}

#' Sequential color palette - turbo
#' 
#' Preset continuous color palettes from pals
#' 
#' @details
#'     These palettes were generated using the \code{pals} package by 
#'     Kevin Wright.

.cpal_seq_turbo <- function() {
  c("#30123B", "#3B3081", "#424DB8", "#456ADF", "#4485F6", "#39A0FA", 
    "#29BAEA", "#1AD1D1", "#1AE4B6", "#31F197", "#56F975", "#7EFD55", 
    "#A1FB3E", "#BDF235", "#D7E335", "#EBCF39", "#FABA39", "#FD9D2D", 
    "#FA7F20", "#F06114", "#E3460B", "#D03005", "#B81E01", "#9C1001", 
    "#7A0403")
}


#' Sequential color palette - sunset
#' 
#' Preset continuous color palettes from colorspace
#' 
#' @details
#'     These palettes were generated using the \code{colorspace} package by 
#'     Ross Ihaka, Paul Murrell, Kurt Hornik, Hason Fisher, Reto Stauffer, 
#'     Claus Wilke, Claire McWhite, and Achim Zeileis.

.cpal_seq_sunset <- function(){
  c("#704D9E", "#7F4EA2", "#8D4FA5", "#9A51A7", "#A653A8", "#B156A9", 
    "#BC5AA9", "#C65EA8", "#CF63A6", "#D868A3", "#E06EA0", "#E7759C", 
    "#ED7C97", "#F18592", "#F48E8E", "#F6978A", "#F7A086", "#F9A983", 
    "#F9B282", "#F9BB81", "#F9C483", "#F8CD86", "#F6D68A", "#F5DE91", 
    "#F3E79A")
}

#' Sequential color palette - YlGnBu truncated
#' 
#' Preset continuous color palettes from colorspace
#' 
#' @details
#'     These palettes were generated using the \code{colorspace} package by 
#'     Ross Ihaka, Paul Murrell, Kurt Hornik, Hason Fisher, Reto Stauffer, 
#'     Claus Wilke, Claire McWhite, and Achim Zeileis.

.cpal_seq_ylgnbu <- function(){
  c("#C3ECC7", "#BAE8C4", "#B0E5C1", "#A5E1BE", "#9ADCBB", "#8DD8B9", 
    "#80D4B7", "#71CFB4", "#60CAB3", "#4BC5B1", "#30C0B0", "#00BAAF", 
    "#00B5AE", "#00AFAE", "#00A9AE", "#00A3AE", "#009CAF", "#0095AF", 
    "#008FB0", "#0087B1", "#0080B3", "#0076AF", "#006CA9", "#0062A2", 
    "#00579A", "#004D91", "#004288", "#00387E", "#152E74", "#21236A", 
    "#26185F")
}

#' Sequential color palette - Heat 2
#' 
#' Preset continuous color palettes from colorspace
#' 
#' @details
#'     These palettes were generated using the \code{colorspace} package by 
#'     Ross Ihaka, Paul Murrell, Kurt Hornik, Hason Fisher, Reto Stauffer, 
#'     Claus Wilke, Claire McWhite, and Achim Zeileis.

.cpal_seq_heat2 <- function(){
  c("#E2E6BD", "#E3E37C", "#E4DD6B", "#E5D75D", "#E6D050", "#E7CA45", 
    "#E8C33C", "#E9BD34", "#E9B62D", "#EAAF29", "#EAA828", "#E9A129", 
    "#E99A2C", "#E89330", "#E78C36", "#E6853B", "#E57E41", "#E37747", 
    "#E1704C", "#DF6852", "#DD6157", "#DB595C", "#D85161", "#D54865", 
    "#D33F6A")
}
