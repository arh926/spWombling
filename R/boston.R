#' Corrected Boston Housing Data
#'
#' The boston.c data frame has 506 rows and 20 columns. It contains the Harrison and Rubinfeld (1978) data corrected for a few minor errors and augmented with the latitude and longitude of the observations. Gilley and Pace also point out that MEDV is censored, in that median values at or over USD 50,000 are set to USD 50,000. The original data set without the corrections is also included in package mlbench as BostonHousing. In addition, a matrix of tract point coordinates projected to UTM zone 19 is included as boston.utm, and a sphere of influence neighbours list as boston.soi.
#'
#' @format Refer to spData for full documentation
#' @usage data(boston.c)
#' @source spData
"boston.c"