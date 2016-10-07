\name{primers}
\alias{primers}
\docType{data}
\title{
primers
}
\description{
primers designed using Primer3 for sample data
}
\usage{data("primers")}
\format{
  A data frame with 5 observations on the following 28 variables.
  \describe{
    \item{\code{i}}{a numeric vector}
    \item{\code{PRIMER_LEFT_SEQUENCE}}{a character vector}
    \item{\code{PRIMER_RIGHT_SEQUENCE}}{a character vector}
    \item{\code{PRIMER_LEFT_TM}}{a numeric vector}
    \item{\code{PRIMER_RIGHT_TM}}{a numeric vector}
    \item{\code{PRIMER_LEFT_pos}}{a numeric vector}
    \item{\code{PRIMER_LEFT_len}}{a numeric vector}
    \item{\code{PRIMER_RIGHT_pos}}{a numeric vector}
    \item{\code{PRIMER_RIGHT_len}}{a numeric vector}
    \item{\code{PRIMER_PAIR_PENALTY}}{a numeric vector}
    \item{\code{PRIMER_LEFT_PENALTY}}{a numeric vector}
    \item{\code{PRIMER_RIGHT_PENALTY}}{a numeric vector}
    \item{\code{PRIMER_LEFT_GC_PERCENT}}{a numeric vector}
    \item{\code{PRIMER_RIGHT_GC_PERCENT}}{a numeric vector}
    \item{\code{PRIMER_LEFT_SELF_ANY_TH}}{a numeric vector}
    \item{\code{PRIMER_RIGHT_SELF_ANY_TH}}{a numeric vector}
    \item{\code{PRIMER_LEFT_SELF_END_TH}}{a numeric vector}
    \item{\code{PRIMER_RIGHT_SELF_END_TH}}{a numeric vector}
    \item{\code{PRIMER_LEFT_HAIRPIN_TH}}{a numeric vector}
    \item{\code{PRIMER_RIGHT_HAIRPIN_TH}}{a numeric vector}
    \item{\code{PRIMER_LEFT_END_STABILITY}}{a numeric vector}
    \item{\code{PRIMER_RIGHT_END_STABILITY}}{a numeric vector}
    \item{\code{PRIMER_LEFT_TEMPLATE_MISPRIMING}}{a numeric vector}
    \item{\code{PRIMER_RIGHT_TEMPLATE_MISPRIMING}}{a numeric vector}
    \item{\code{PRIMER_PAIR_COMPL_ANY_TH}}{a numeric vector}
    \item{\code{PRIMER_PAIR_COMPL_END_TH}}{a numeric vector}
    \item{\code{PRIMER_PAIR_PRODUCT_SIZE}}{a numeric vector}
    \item{\code{PRIMER_PAIR_TEMPLATE_MISPRIMING}}{a numeric vector}
  }
}
\value{
Dataframe of primer design results
}
\examples{
data(primers)
}
\keyword{datasets}
