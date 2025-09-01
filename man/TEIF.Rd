\name{TEIF}
\alias{TEIF}
\title{Simple thermodynamic equilibrium isotope fractionation}
\usage{
TEIF(source, mass, epsilon, fraction, Vr, K)
}
\arguments{
\item{source}{a source object created with the \code{\link{create_Source}} function containing the isotopic composition of the source being fractionated}
\item{mass}{The atomic mass of the stable isotope to which the isotope discrimination factor (below) refers to}
\item{epsilon}{the isotope discrimination factor of the isotope fractionation processes, expressed in \eqn{\unicode{x2030}}}
\item{fraction}{a character ("r" or "p") determining whether to output the isotopic composition of the reactant or product of the isotope fractionation reaction at thermodynamic equilibrium}
\item{Vr}{the ratio of the volume dimension of the reaction products over the volume dimension of the reactants}
\item{K}{the thermodynamic equilibrium constant of the reference stable isotope (most commonly the lighter isotope, denominator in the isotop eratio)}
}
\value{
a source object (such as created with the \code{\link{create_Source}} function) containg the isotopic composition of fractionated input source.  
}
\description{
This function calculates the isotopic composition of a source after the occurrence of a thermodynamic equilibrium isotope fractionation process assuming that no reaction products are initially present in the system. Isotope fractionation of all stable isotopes of the element is calculated assuming mass-dependent isotope fractionation.
}
\examples{
# Create a Magnesium source object:
S = create_Source(name = "Source_Mg", element = "Mg", size = 10)
# Thermodynamic equilirbium isotope fractionation of source "S"
S2 = TEIF(source = S, mass = 26, epsilon = 1, fraction = "p", Vr = 1, K = 0.001)
S3 = TEIF(source = S, mass = 26, epsilon = 1, fraction = "r", Vr = 1, K = 0.001)
}