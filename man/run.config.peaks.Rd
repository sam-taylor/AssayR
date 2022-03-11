\name{run.config.peaks}
\alias{run.config.peaks}
\title{run.config.peaks}
\usage{
run.config.peaks(path.to.config.tsv='config.tsv', path.to.tics, Interactive=TRUE, chrom.pattern="std|standard", chrom.match.mode="omit", shorten.col.names=TRUE)
}
\arguments{
    \item{path.to.config.tsv}{path to a tab-delimited table of config data, see below}
    \item{path.to.tics}{path to the tics files generated during run.config.tics}
    \item{Interactive}{optionally run interactively}
    \item{chrom.pattern}{pattern to match for including or omitting samples in peak detection phase}
    \item{chrom.match.mode}{whether to omit or include samples whose name matches chrom.pattern}
    \item{shorten.col.names}{whether to attempt semi-intelligent shortening of sample names}
}
\description{
run.config.peaks take the config definitions and the path to
extracted ion chromatograms, does interactive peak detection
and quantifies isotopes in each sample.
The results are returned in a data.frame.
}
\examples{
p <- run.config("tca-cycle.tsv", "20170526-161533/")
write.csv(p, "results.csv")
}
