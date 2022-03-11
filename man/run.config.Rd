\name{run.config}
\alias{run.config}
\title{run.config}
\usage{
run.config(path.to.config.tsv='config.tsv', path.to.mzMLs='mzMLneg/', chrom.pattern="std|standard", chrom.match.mode="omit", shorten.col.names=TRUE)
}
\arguments{
    \item{path.to.config.tsv}{path to a tab-delimited table of config data, see below}
    \item{path.to.mzMLs}{a directory containing the mzML files to be analysed}
    \item{chrom.pattern}{pattern to match for including or omitting samples in peak detection phase}
    \item{chrom.match.mode}{whether to omit or include samples whose name matches chrom.pattern}
    \item{shorten.col.names}{whether to attempt semi-intelligent shortening of sample names}
}
\description{
a shortcut for calling run.config.tics and run.config.peaks
}
\examples{
p <- run.config("tca-cycle.tsv", "mzMLneg/")
write.csv(p, "results.csv")
}
