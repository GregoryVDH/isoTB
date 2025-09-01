\name{min_Enrich}
\alias{min_Enrich}
\title{Calculate the minimum enrichment of a tracer source to be able to neglect isotopic fractionation}
\usage{
min_Enrich(er.tar, prop.T, Source_N, mass, alpha, frac.type, reac.par)
}
\arguments{
\item{er.tar}{the target error}
\item{prop.T}{proportion of the tracer source in the isotopic mixing model}
\item{Source_N}{The natural abundance source in the isotope mixing model: a source object created with the \code{\link{create_Source}} function}
\item{mass}{The atomic mass of the enriched stable isotope}
\item{alpha}{a vector of numeric values (dim = 3) containing the isotopic discrimination factors for the isotopic fractionation:
		\enumerate{
			\item the natural abundance source
			\item the enriched tracer source
			\item the mix of the natural abundance and enriched tracer sources 
		} 
	}
\item{frac.type}{an array of strings (dim = 3) containing the types of isotopic fractionation (thermodynamic equilibrium or first order kinetic) for each of the three isotopic fractionation processes: 
		\enumerate{
			\item the natural abundance source
			\item the enriched tracer source
			\item the mix of the natural abundance and enriched tracer sources 
		}
		Possible values are: 
		\itemize{
			\item "TEIF" (thermodynamic equilibrium isotopic fractionation)
			\item "FOKIF" (first order kinetic isotope fractionation)
			}
		}
\item{reac.par}{an array of numeric values (dim = 3) containing the reaction parameters for each of the three isotopic fractionation processes: 
		\enumerate{
			\item the natural abundance source
			\item the enriched tracer source
			\item the mix of the natural abundance and enriched tracer sources 
		}
		For thermodynamic equilibrium isotopic fractionation, input the thermodynamic equilibrium constant of the reference isotope. For first order kinetic isotopic fractionation, input the fraction f of reactants remaining in the system.
		}
}
\value{
The minimum isotope enrichment, expressed in \eqn{\delta} notation (i.e. isotope ratio relative deviation from the international standard), in order to ensure that neglecting isotopic fractionation causes isotope mixing model errors below the target error set as an argument to the function 
}

\description{
In a system with two sources N and T that mix together with a proportion of prop.T of source T and (1 - prop.T) of source Nin which three isotopic fractionation processes are taken into account:
	\enumerate{
		\item fractionation of the natural abundance source
		\item fractionation of the enriched tracer source
		\item fractionation of the mix of the natural abundance and enriched tracer sources 
	}
This function is solver algorythm to find the minimum enrichment of the tracer source T for which the error induced by neglecting isotopic fractionation processes (calculated with \code{\link{error_calc}}) is inferior to a user define error target.
}
\seealso{
	\code{\link{error_calc}}, \code{\link{min_Enrich}}
}

