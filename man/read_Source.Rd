\name{read_Source}
\alias{read_Source}
\title{Extract the isotopic composition from a source type object}
\usage{
read_Source(source, mass, notation)
}
\arguments{
\item{source}{a source type object created with \code{\link{create_Source}} function}
\item{mass}{the atomic mass for the desired output stable isotope}
\item{notation}{the desired output isotopic composition notation. Possible values are:
	\itemize{
		\item "F" (fractional abundance)
		\item "R" (Isotope Ratio)
		\item "delta" (delta notation)
	}
}
}
\value{
a numeric value
}
\description{
This function extracts the isotopic composition of a selected isotope and selected isotopic composition notation
}
\examples{
S = create_Source("S", element = "Mg", size = 2.5)
S = TEIF(source = S, mass = 26, epsilon = 5, fraction = "p", Vr = 1, K = 0.001)
read_Source(S, mass = 26, notation = "delta")
read_Source(S, mass = 26, notation = "R")
read_Source(S, mass = 26, notation = "F")
read_Source(S, mass = 25, notation = "delta")
}
