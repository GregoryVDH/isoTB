\name{Frac_from_Delta}
\alias{Frac_from_Delta}
\title{Isotope fractionation calculated from a target isotopic composition difference from a source (i.e. \eqn{\Delta})}
\usage{
Frac_from_Delta(source, mass, Delta, frac.type)
}
\arguments{
\item{source}{a source object created with the \code{\link{create_Source}} function containing the isotopic composition of the source being fractionated}
\item{mass}{the atomic mass of the nominator stable isotope in the isotope ratio (commonly the heavy isotope)}
\item{Delta}{The isotopic difference between the isotopic composition before and after isotope fractionation (unit ‰)}
\item{frac.type}{a string containing the type of isotopic fractionation process. Possible values are: 
		\itemize{
			\item "TEIF" (thermodynamic equilibrium isotopic fractionation)
			\item "FOKIF" (first order kinetic isotope fractionation)
			}
		}
}
\value{
a source object created with the \code{\link{create_Source}} function containing the isotopic composition of the fractionated source
}
\description{
This function calculates the isotopic composition of a reservoir after the occurrence of an isotope fractionation process from an input isotopic difference (i.e. \eqn{\Delta}) between the reservoir and its original source (before fractionation). This function is similar to the \code{\link{TEIF}} and \code{\link{FOKIF}} functions but calculates the istope fractionation from the isotopic difference instead of the isotope fractionation/discrimination factor.
Isotopic composition is calculated for all stable isotopes of the selected element.
}
\examples{
S = isoTB::create_Source("S", element = "Mg", size = 1)
Frac_from_Delta(source = S, mass = 26, Delta = + 5, frac.type = "TEIF")
}
