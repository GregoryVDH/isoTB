\name{Delta_to_F}
\alias{Delta_to_F}
\title{Convert delta notation to fractional abundance for enriched/spiked isotopic composition}
\usage{
Delta_to_F(delta, element, mass)
}
\arguments{
\item{delta}{the isotopic composition expressed in \eqn{\delta} notation (unit ‰) to convert. Accepted values are single value or an array of values or a dataframe (see description for the required dataframe structure)}
\item{element}{a string containing the chemical element name (e.g. C, N, O, Ca, Mg, etc.) which is used to extract the stable isotope data from the \code{\link{isoTab}} package dataset.}
\item{mass}{the atomic mass of the nominator stable isotope in the isotope ratio (commonly the heavy isotope). If a dataframe is input as the \code{delta} argument, the `mass` is no longer required}
}
\value{
Fractional abundance. Outputs vary depending on the format of the input \code{delta} argument:
\itemize{
	\item{a single value is output if \code{delta} is a single value}
	\item{an array of values is output if \code{delta} is an array}
	\item{a dataframe of values is output if \code{delta} is a dataframe. The fractional abundance of each stable isotope is output in individual colums named "F_XX", where XX corresponds to the atomic mass of the stable isotope}
}
}
\description{
Converts the delta notation (i.e. relative deviation from the isotope ratio of the international standard) noted \eqn{\delta^{m_2/m_1}X} and expressed in ‰ into a fractional abundance (between 0 and 1). The function can be used in two different ways:
		\enumerate{
			\item{to calculate the fractional abundance of a selected enriched/spiked isotope. Delta-notation isotope enrichments for this isotope are input as a single value or an array of values through the \code{delta} argument and fractional abundances are calculated assuming that the isotope ratios of all aother stable isotopes are equal to their isotope ratios in the standard reference material.}
			\item{to calculate the fractional abundance of all stable isotopes of a given element (best adapted for natural variations of isotope ratios). The \code{delta} argument requires in this case a dataframe containing in different columns the delta-notation isotope enrichments for each stable isotope ratio (see example 2) and ordered by increasing isotope mass (e.g. 25/24, 26/24 or 42/40, 43/40, 44/40, 46/40, 48/40). Rows in the dataframe correspond to different samples/measurements/reservoirs. If the format of the input dataframe is incorrect, the function throws back and "invalid delta argument" error.}
		}
}
\examples{
# Example 1:
	# Magnesium
	Delta_to_F(delta = 100, element = "Mg", mass = 26) 
	Delta_to_F(delta = c(0, 100, 500), element = "Mg", mass = 26) 
	Delta_to_F(delta = c(0, 50, 1000), element = "Mg", mass = 25) 
	# Calcium
	Delta_to_F(delta = 450, element = "Ca", mass = 44) 

# Example 2:
	# Magnesium:
	df <- data.frame(delta25_24 = c(0, 10, 100), delta26_24 = c(0, 50, 1000))
	Delta_to_F(delta = df, element = "Mg")
}
