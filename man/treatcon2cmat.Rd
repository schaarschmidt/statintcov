\name{treatcon2cmat}
\alias{treatcon2cmat}

\title{Check or construct between-treatment contrast matrix as defined in contrMat, package multcomp. For internal use.}

\description{Helper function to check a given (between-treatment) contrast matrix or contruct a suitable matrix given one of the character string available in contrMat(package multcomp). For internal use}

\usage{
treatcon2cmat(treatcon, base, ntab)
}

\arguments{
  \item{treatcon}{either a matrix with numeric entries, where each row specifies one contrasts between treatments, or a single character string; if a matrix: number of columns must equal the number of treatment levels in the factor variable
  if a character string: chracter string will be passed to the type-argument in \code{contrMat} (package multcomp); see ?contrMat for details
}
  \item{base}{a single integer, indicating which treatment level shoudl serve as a control, when passed to \code{contrMat} (package multcomp) used if treatcon is a characted string, ignored, if treatcon is a matrix
}
  \item{ntab}{a (named) vector with sample sizes in the treatment groups, to be passed to argument \code{n} in \code{contrMat}; ignored if treatcon is character (makes a difference for pooling contrasts only, i.e. "GrandMean", "AVE", "Williams" etc.)
}
}

\value{
a matrix, for internal use
}

\keyword{misc}
