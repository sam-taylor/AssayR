\name{run.config.tics}
\alias{run.config.tics}
\title{run.config.tics}
\usage{
run.config.tics(path.to.config.tsv='config.tsv', path.to.mzMLs='mzMLneg/')
}
\arguments{
    \item{path.to.config.tsv}{path to a tab-delimited table of config data, see below}
    \item{path.to.mzMLs}{a directory containing the mzML files to be analysed}
}
\description{
run.config.tics calls msaccess on mzML files to exctract ion
chromatograms according to the definitions in mzML.
Returns the new path that the mzMLs were written to (which can
then be passed to run.config.peaks)
}
\examples{
path.to.tics <- run.config.tics("tca-cycle.tsv", "mzMLneg/")
}
