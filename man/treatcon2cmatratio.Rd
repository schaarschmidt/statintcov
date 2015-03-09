\name{treatcon2cmatratio}
\alias{treatcon2cmatratio}

\title{Check or construct a matrices defining multuiple ratios of interest.}

\description{
Helper function to check a given object defining (between-treatment) comparisons as ratios or contruct a suitable list of two matrices of linear combinations
for numerator and denominator given one of a character string available in \code{contrMatRatio} in package mratios. For internal use.
}

\usage{
treatcon2cmatratio(treatcon, base, ntab)
}

\arguments{
  \item{treatcon}{either a list of 2 matrices with numeric entries, or a single character string; if a list of 2 matrices: number of columns must be = number of treatment levels in the factor variable rows specifies contrasts between treatments; if a character string: character string will be passed to the type-argument in \code{contrMatRatio} (package mratios) see \code{?contrMatRatio} for details
}
  \item{base}{a single integer, indicating which treatment level shoudl serve as a control, when passed to \code{contrMat} (package multcomp) used if treatcon is a characted string, ignored, if treatcon is a matrix
}
  \item{ntab}{a (named) vector with sample sizes in the treatment groups, to be passed to argument \code{n} in \code{contrMat}; ignored if treatcon is character (makes a difference for pooling contrasts only, i.e. "GrandMean", "AVE", "Williams" etc.)
}
}

\value{
a list containing two matrices to be used in \code{gsci.ratio} (package mratios)
}


\keyword{misc}
