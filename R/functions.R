#' @importFrom forcats fct_collapse

# check sex contains only 1/M/B/TRUE or 2/F/G
test_sex <- function(sex) {
  fsex <- toupper(substr(sex, 1, 1))
  fsex <- factor(fsex, levels = c(1:2, 'M', 'F', 'B', 'G', 'T'))
  levels(fsex) <- c(rep(1:2, 3), 1)
  fsex <- fct_collapse(fsex, boys = '1', girls = '2')
  droplevels(fsex)
}

# return sitar code for reference
ref_sitar <- function(x) {
  x <- unique(x)
  stopifnot('cutoff not unique' = length(x) == 1L)
  x <- sub('^(.*) .*$', '\\1', x) # drop any cutoff
  f <- factor(x, levels = c('CDC', 'IOTF', 'WHO'))
  levels(f) <- c('cdc2000', 'iotf', 'who0607')
  as.character(f)
}

