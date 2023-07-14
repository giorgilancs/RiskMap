##' @title Malaria Transmission in the Western Kenyan Highlands
##' @description The dataset contains information on 82014 individuals enrolled
##' in concurrent school and community cross-sectional surveys, conducted in 46
##' school clusters in the western Kenyan highlands. Malaria was assessed by
##' rapid diagnostic test (RDT).
##'
##' The variables are as follows:
##' \itemize{
##'  \item Cluster: unique ID for each of the 46 school clusters.
##'  \item Long: longitude coordinate of the household location.
##'  \item Lat: latitude coordinate of the household location.
##'  \item RDT: binary variable indicating the outcome of the RDT:
##'  1, if positive, and 0, if negative.
##'  \item Gender: factor variable indicating the gender of the sampled individual.
##'  \item Age: age in years of the sampled individual.
##'  \item NetUse: binary variable indicating whether the sampled individual
##'  slept under a bed net the previous night: 1, if yes, 0, if no.
##'  \item MosqCntl: binary variable indicating whether the household has used some kind
##'   of mosquito control, such as sprays and coils: 1, if yes, 0, if no.
##'  \item IRS: binary variables in indicating whether there has been indoor
##'  residual spraying (IRS) in the house in the last 12 months: 1, if yes, 0,
##'  if no.
##'  \item Travel:  binary variable indicating whether the sampled individual
##'  has travelled outside the village in the last three months: 1, if yes, 0,
##'  if no.
##'  \item SES: ordinal variable indicating the socio-economic status (SES) of
##'   the household. The variables is an integer score from 1(=poor) to 5(=rich).
##'  \item District: factor variable indicating the village of the sampled
##'  individual, "Kisii Central" or "Rachuonyo".
##'  \item Survey: factor variables indicating the survey in which the
##'  participant was enrolled, "community" or "school".
##' }
##' @docType data
##' @keywords datasets
##' @name malkenya
##' @usage data(malkenya)
##' @format A data frame with 82014 rows and 13 variables
##' @source Stevenson, J.C., Stresman, G.H., Gitonga, C.W., Gillig, J.,
##' Owaga, C., et al. (2013). Reliability of School Surveys in Estimating Geographic
##' Variation in Malaria Transmission in the Western Kenyan Highlands.
##' PLOS ONE 8(10): e77641. doi: 10.1371/journal.pone.0077641
NULL
