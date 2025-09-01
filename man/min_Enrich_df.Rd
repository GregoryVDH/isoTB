\name{min_Enrich_df}
\alias{min_Enrich_df}
\title{Calculate the minimum enrichment of a tracer source to be able to neglect isotopic fractionation over a range of isotope mixing proportions between the natural abundance and enriched tracer sources}
\usage{
min_Enrich_df(er.tar, Source_N, mass, alpha, frac.type, reac.par, range, res)
}
\arguments{
\item{er.tar}{the target error}
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
	\item{range}{an array of numeric values (dim = 2) comprised between 0 and 1 setting the tracer source mixing proportion minimum and maximum values of the range over which to calculate the minimum enrichment of the tracer source}
	\item{res}{a numeric value defining the resolution/step to discreticize the range of defined above}
}

\value{
A data frame containing for different mixing proportions between the natural abundance source and the tracer source the minimum isotope enrichment, expressed in \eqn{\delta} notation (i.e. isotope ratio relative deviation from the international standard), in order to ensure that neglecting isotopic fractionation causes isotope mixing model errors below the target error set as an argument to the function. The data frame contains two columns named:
		\item{Tracer.source}{the proportion of the enriched tracer source in the isotope mixing model (ranging between 0 and 1)}
		\item{min.Enrich}{the minimum enrichment calculated for the given proportion of isotope mixing (Tracer.Source) to ensure an error inferior to er.tar when neglecting isotopic fractionation}
}

\description{
In a system with two sources N and T that mix together with a proportion of prop.T of source T and (1 - prop.T) of source Nin which three isotopic fractionation processes are taken into account:
	\enumerate{
		\item fractionation of the natural abundance source
		\item fractionation of the enriched tracer source
		\item fractionation of the mix of the natural abundance and enriched tracer sources 
	}
This function calculates the minimum enrichment of the tracer source T for which the error induced by neglecting isotopic fractionation processes (calculated with \code{link{error_calc}} is inferior to a user define error target over a range of isotope mixing proportions between the sources N and T.
}
\seealso{
	\code{\link{error_calc}}, \code{\link{min_Enrich}}
}

