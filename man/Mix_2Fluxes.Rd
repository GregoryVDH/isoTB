\name{Mix_2Fluxes}
\alias{Mix_2Fluxes}
\title{2-source isotope mixing model of two fluxes}
\usage{
Mix_2Fluxes(Source1, Source2, Flux1, Flux2)
}
\arguments{
\item{Source1}{Source 1 in the isotope mixing model containg the isotopic composition of Flux1: a source object created with the \code{\link{create_Source}} function
}
\item{Source2}{Source 2 in the isotope mixing model containg the isotopic composition of Flux2: a source object created with the \code{\link{create_Source}} function
}
\item{Flux1}{Size of Flux1 (a unitless numerical value). Ensure the homogeniety of units used for Flux1 and Flux 2 size arguments}
\item{Flux2}{Size of Flux2 (a unitless numerical value). Ensure the homogeniety of units used for Flux1 and Flux 2 size arguments}
}

\value{
a source object as that created with the \code{\link{create_Source}} function containing the isotopic composition of the isotope mixing of Flux1 and Flux2. The size property of the output object is the sum of Flux1 and Flux2.
}

\description{
Stable isotope mixing model of two sources (1 and 2) with mixing proportions defined by two values of fluxes (1 and 2). The isotope mixing model is based on the fractional abundances mixing equations.
}
\examples{
S1 <- create_Source("S1", element = "Mg", size = 1)
S1 <- write_Source(source = S1, value = -5, mass = 26, notation = "delta")
S2 <- create_Source("S2", element = "Mg", size = 1)
S2 <- write_Source(source = S2, value = 100, mass = 26, notation = "delta")
Mix_2Fluxes(Source1 = S1, Source2 = S2, Flux1 = 10, Flux2 = 25)
}
\seealso{
	\code{\link{Mix_Source}}
}

