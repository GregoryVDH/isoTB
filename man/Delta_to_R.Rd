\name{Delta_to_R}
\alias{Delta_to_R}
\title{Convert delta notation to an isotope ratio}
\usage{
Delta_to_R(delta, element, mass, Rref)
}
\arguments{
\item{delta}{the isotopic composition expressed in \eqn{\delta} notation (unit ‰) to convert. Accepted input formats are single values and arrays of values}
\item{element}{a string containing the chemical element name (e.g. C, N, O, Ca, Mg, etc.) which is used to extract the stable isotope data from the \code{\link{isoTab}} package dataset.}
\item{mass}{the atomic mass of the nominator stable isotope in the isotope ratio (commonly the heavy isotope)}
\item{Rref}{Standard reference material isotope ratio (optional argument). If not supplied, \code{Delta_to_R()} calculates isotope ratios from the standard reference material isotope ratio in \code{\link{isoTab}}. If \code{Rref} is supplied, \code{element} and \code{mass} arguments are not required}
}
\value{
Isotope ratio of the selected stable isotope
}
\description{
Converts the delta notation (i.e. relative deviation from the isotope ratio of the international standard) noted \eqn{\delta^{m_2/m_1}X} and expressed in ‰ into its equivalent isotope ratio
}
\examples{
# Magnesium
Delta_to_R(delta = 100, element = "Mg", mass = 26) 
Delta_to_R(delta = c(0, 100, 500), element = "Mg", mass = 26) 
Delta_to_R(delta = 100, Rref = 0.140)
Delta_to_R(delta = c(0, 100, 500), Rref = 0.140)
# Calcium
Delta_to_R(delta = 450, element = "Ca", mass = 44) 
}
