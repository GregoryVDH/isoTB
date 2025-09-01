\name{Tracer_prop}
\alias{Tracer_prop}
\title{Proportion of a reservoir originating from an enriched tracer source}
\usage{
Tracer_prop(delta.samp, delta.nat, delta.tracer, element, mass, F.tracer) 
}
\arguments{
\item{delta.samp}{the isotope enrichment of the sampled reservoir expressed in \eqn{delta} notation, unit \eqn{\unicode{x2030}}}
\item{delta.nat}{the isotope enrichment of the natural abundance source in the isotope mixing model expressed in \eqn{delta} notation, unit \eqn{\unicode{x2030}}}
\item{delta.tracer}{the isotope enrichment of the enriched tracer source in the isotope mixing model expressed in \eqn{delta} notation, unit \eqn{\unicode{x2030}}. If the `delta.tracer` argument is not supplied then the isotope enrichment of the enriched tracer source needs to be supplied as a fractional abundance with argument `F.tracer`}
\item{element}{a string containing the chemical element name (e.g. C, N, O, Ca, Mg, etc.) which is used to extract the stable isotope data from the \code{\link{isoTab}} package dataset.}
\item{mass}{the atomic mass of the enriched stable isotope}
\item{F.tracer}{fractional abundance of the enriched isotope in the tracer source. Argument used if supplied. If absent in the list of arguments  the `delta.tracer` argument is used}
}
\value{
a numerical value between 0 and 1.
}
\description{
Solves a two-source isotope mixing model (natural abundance source and enriched tracer source) to determine the proportion of a given reservoir that originates from th e enriched tracer source.
}
\examples{
Tracer_prop(delta.samp = 100, delta.nat = 0, delta.tracer = 500, element = "Mg", mass = 26)
Tracer_prop(delta.samp = 100, delta.nat = 0, element = "Mg", mass = 26, F.tracer = 0.3)
}
