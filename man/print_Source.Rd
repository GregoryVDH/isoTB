\name{print_Source}
\alias{print_Source}
\title{Print the isotopic composition of a source type object}
\usage{
print_Source(source)
}
\arguments{
\item{source}{a source type object created with \code{\link{create_Source}} function}
}
\value{
NULL
}
\description{
This function prints a dataframe containg the isotopic composition of the input source, displaying the fractional abundances, isotope ratios and delta notations
}
\examples{
S = create_Source("S", element = "Mg", size = 2.5)
S = TEIF(source = S, mass = 26, epsilon = 5, fraction = "p", Vr = 1, K = 0.001)
print_Source(S)
}
