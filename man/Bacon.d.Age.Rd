% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc.R
\name{Bacon.d.Age}
\alias{Bacon.d.Age}
\title{Output all depths for a single age.}
\usage{
Bacon.d.Age(
  age,
  set = get("info"),
  BCAD = set$BCAD,
  its = set$output,
  na.rm = FALSE
)
}
\arguments{
\item{age}{The age estimate for which depths are to be returned. Has to be a single age.}

\item{set}{Detailed information of the current run, stored within this session's memory as variable info.}

\item{BCAD}{The calendar scale of graphs and age output-files is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}.}

\item{its}{The set of MCMC iterations to be used. Defaults to the entire MCMC output, \code{its=set$output}.}

\item{na.rm}{Whether or not to remove NA values (ages within slumps)}
}
\value{
Outputs all MCMC-derived depths for a given age.
}
\description{
Output all MCMC-derived depth estimates for a single given age.
}
\details{
Obtaining an age-depth model is often only a step towards a goal, e.g., plotting a core's
fossil series ('proxies') against calendar time. Bacon.d.Age can be used to list all MCMC-derived depths belonging to a given (single) age, for example to calculate mean depths belonging to a modelled depth. 
This function was kindly written and provided by Timon Netzel (Bonn University). See also Bacon.Age.d, which calculates the ages for a single depth.
}
\examples{
\dontrun{
  Bacon(run=FALSE, coredir=tempfile())
  agedepth(age.res=50, d.res=50, d.by=10)
  ages.d20 = Bacon.Age.d(20)
  mean(ages.d20)
}
}
\author{
Maarten Blaauw, J. Andres Christen
}
