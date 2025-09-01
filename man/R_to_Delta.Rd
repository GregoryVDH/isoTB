\name{R_to_Delta}
\alias{R_to_Delta}
\title{Convert isotope ratio to delta notation}
\usage{
R_to_Delta(R, element, mass, Rref)
}
\arguments{
\item{R}{The isotope ratio \eqn{{}^{m_2/m_1}X} to be converted to delta notation}

\item{element}{a string containing the chemical element name (e.g. C, N, O, Ca, Mg, etc.) which is used to extract the stable isotope data from the \code{\link{isoTab}} package dataset.}
\item{mass}{the atomic mass of the nominator stable isotope in the isotope ratio (commonly the heavy isotope)}
\item{Rref}{The isotope ratio of the international standard}
}
\value{
Isotope ratio relative deviation from the international standard, expressed in \eqn{\delta} notation (unit ‰).
}
\description{
Converts the isotope ratio \eqn{{}^{m_2/m_1}X} to the delta notation (i.e. relative deviation from the isotope ratio of the international standard) noted \eqn{\delta^{m_2/m_1}X} and expressed in ‰. The value of the international standard isotope ratio can be an input to the function (Rref) and thus set manually (Example 1) or it can be derived from the \code{\link{isoTab}} package dataset (Example 2). 
}
\examples{
# Example 1: Convert Magnesium 26/24 isotope ratio to delta notation 
# by setting manually the reference isotope ratio:
R_to_Delta(R = 0.1394, Rref = 0.13932)
# Example 2: Convert Magnesium 26/24 isotope ratio to delta notation 
# using the reference isotope ratio in the isoTab dataset
R_to_Delta(R = 0.1394, element = "Mg", mass = 26)
}
