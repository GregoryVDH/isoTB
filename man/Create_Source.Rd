\name{create_Source}
\alias{create_Source}
\title{Create a source type list object}
\usage{
create_Source(name, element, size)
}
\arguments{
\item{name}{a string containing the name of the source object}
\item{element}{a string containing the chemical element name (e.g. C, N, O, Ca, Mg, etc.) which is used to extract the stable isotope data from the \code{\link{isoTab}} package dataset.}
\item{size}{a numeric value containing the size of the quantity of the selected element in the source/reservoir (unitless). Ensure the homogeniety of units used for size arguments when creating multiple source objects}
}
\value{
A Source type object
}
\description{
Creates a source type list object
}
\examples{
# Create a Magnesium source object:
create_Source(name = "Source_Mg", element = "Mg", size = 10)
# Create a Calcium source object:
create_Source(name = "Source_Ca", element = "Ca", size = 0.725)
# Create a carbon source object:
create_Source(name = "Source_C", element = "C", size = 1.46*1E3)
}
