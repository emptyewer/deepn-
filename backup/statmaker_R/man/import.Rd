\name{import}
\alias{import}
\title{Imports data from CSV files}
\description{
Accepts csv filenames and reads them into a structured list for use in
downstream analyses.
}
\usage{
import(vn, vs, bn, bs)
}
\arguments{
  \item{vn}{Vector of filenames for empty (v)ector under (n)onselective
  conditions}
  \item{vs}{Vector of filenames for empty (v)ector under (s)elective
  conditions}
  \item{bn}{Vector of filenames for (b)ait under (n)onselective conditions}
  \item{bs}{Vector of filenames for (b)ait under (s)selective conditions}
}
\value{
This function is essentially a wrapper for the \code{\link{readDEEPN}} function,
returning a list of objects:
  \item{Vector}{Empty vector counts}
  \item{Bait}{Bait counts}
  \item{multiBait}{Logical indicator for whether multiple baits are involved}
  \item{vtr}{Total mapped reads (vector)}
  \item{btr}{Total mapped reads (bait)}
}
\author{Patrick Breheny}
\seealso{
\code{\link{readDEEPN}},
\code{\link{applyFilter}},
\code{\link{overdisp}}
}
\examples{
# Should set up demo data
\dontrun{
vn <- c("file1.csv", "file2.csv")
vs <- c("file3.csv", "file4.csv")
bn <- "file5.csv"
bs <- "file6.csv"
Data <- import(vn, vs, bn, bs)}
}
