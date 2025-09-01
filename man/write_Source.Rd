\name{write_Source}
\alias{write_Source}
\title{Modify an existing source type list object}
\usage{
write_Source(source, value, mass, notation)
}
\arguments{
\item{source}{the source object created with the \code{\link{create_Source}} function to modify} 
\item{value}{the new value to overwrite the previous value in the source object. The notation and mass arguments define which values in the source object are modified}
\item{mass}{The atomic mass of the stable isotope to modify in the source object}
\item{notation}{A string to define which  isotopic composition notation to modify in the source object. Can only be one of the three following values:
		\itemize{
			\item{"F" - Fractional abundance}
			\item{"R" - Isotope ratio}
			\item{"delta" - Delta notation}
		}
}
}
\value{
A Source type object
}
\description{
Manually modify the isotopic composition of a source object. Use this function to create an articially enriched/spiked tracer source.
}
\examples{
## Create a Magnesium source object:
S = create_Source(name = "Source_Mg", element = "Mg", size = 10)
print(S$F)
print(S$R)
print(S$delta)
S1 = write_Source(source = S, value = 50, mass = 26, notation = "delta")
print(S1$F)
print(S1$R)
print(S1$delta)
S2 = write_Source(source = S, value = 0.150, mass = 25, notation = "R")
print(S2$F)
print(S2$R)
print(S2$delta)
S3 = write_Source(source = S, value = 0.25, mass = 26, notation = "F")
print(S3$F)
print(S3$R)
print(S3$delta)
}
