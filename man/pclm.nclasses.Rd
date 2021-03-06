% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pclm.R
\name{pclm.nclasses}
\alias{pclm.nclasses}
\title{Calculate the number of PCLM internal (raw) classes}
\usage{
pclm.nclasses(x, control = list())
}
\arguments{
\item{x}{Vector with start of the interval for age/time classes.}

\item{control}{List with additional parameters. See \code{\link{pclm.control}}.}
}
\description{
Calculate the number of PCLM internal (raw) classes.
}
\examples{
\dontrun{
# Use a simple data set
x <- c(0, 1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
dx <- c(38, 37, 17, 104, 181, 209, 452, 1190, 2436, 3164, 1852, 307, 13)

# Define the open interval by zero.class.frac
control.1 = list(x.div = 5, zero.class.frac = 0.2, zero.class.end = NULL)
pclm.nclasses(x, control = control.1) #calculate number of raw classes
AU10p.1A <- pclm.default(x, dx, control = control.1)
length(AU10p.1A$raw$x) # the number of raw classes after fit
plot(AU10p.1A)

# Define the open interval by zero.class.end
control.2 = list(x.div = 5, zero.class.end = 109)
pclm.nclasses(x, control = control.2) #calculate the number of raw classes
AU10p.1B <- pclm.default(x, dx, control = control.2)
length(AU10p.1B$raw$x) # the number of raw classes after fit
plot(AU10p.1B)
}
}
\author{
Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
}
\seealso{
\code{\link{pclm.control}}, \code{\link{pclm.interval.multiple}},
}

