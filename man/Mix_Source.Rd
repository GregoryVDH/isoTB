\name{Mix_Source}
\alias{Mix_Source}
\title{2-source isotope mixing model}
\usage{
Mix_Source(Source1, Source2, prop.source1, by.size)
}
\arguments{
\item{Source1}{Source 1 in the isotope mixing model containg the isotopic composition of Flux1: a source object created with the \code{\link{create_Source}} function
}
\item{Source2}{Source 2 in the isotope mixing model containg the isotopic composition of Flux2: a source object created with the \code{\link{create_Source}} function
}
\item{prop.source1}{a numeric value between 0 and 1 giving the proportion of source 1 in the isotope mixing model. If \code{by.size = FALSE}, source 1 and 2 are mixed according to the \code{prop.Source1} argument
}
\item{by.size}{a logical value to determine whether to calculate the isotope mixing model from the `prop.Source1` argument (`by.size = FALSE`) or from the size properties of Source1 and Source2 objects (`by.size = TRUE`). By default, `by.size = FALSE`
}
}

\value{
a source object as that created with the \code{\link{create_Source}} function containing the isotopic composition of the isotope mixing of Flux1 and Flux2. The size property of the output object is the sum of Source1 and Source2 size properties if by.size is set to TRUE. When \code{by.size = FALSE}, The size property of the output object is set by default to the size property of Source1.
}

\description{
Stable isotope mixing model of two sources (1 and 2) with mixing proportions defined either by a proportion argument or by the size property of each of the two source objects. The isotope mixing model is based on the fractional abundances mixing equations.
}
\seealso{
	\code{\link{Mix_2Fluxes}}
}