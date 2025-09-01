\name{error_calc}
\alias{error_calc}
\title{Calculate the error when estimating the contribution of an enriched tracer source to a sink due to neglecting isotopic fractionation processes}
\usage{
error_calc(prop.T, Source_N, Source_T, mass, alpha, frac.type, reac.par) 
}
\arguments{
\item{prop.T}{proportion of the tracer source in the isotopic mixing model}
\item{Source_N}{The natural abundance source in the isotope mixing model: a source object created with the \code{\link{create_Source}} function}
\item{Source_T}{The enriched tracer source in the isotope mixing model: a source object created with the \code{\link{create_Source}} function}
\item{mass}{The atomic mass of the enriched stable isotope}
\item{alpha}{a vector of numeric values (dim = 3) containing the isotopic discrimination factors for the isotopic fractionation processes affecting:
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
(numeric value) The isotope mixing model error due to neglecting isotopic fractionation processes. Positive values indicate an overestimation of the contribution of the tracer source in the isotope mixing model, negative values indicate and underestimation. 
}

\description{
For two sources N and T that mix together with a proportion of prop.T of source T and (1 - prop.T) of source N, the isotopic composition of the sink resulting from the mixing of N and T is calculated taking into account three isotopic fractionation processes affecting:
		\enumerate{
			\item the natural abundance source
			\item the enriched tracer source
			\item the mix of the natural abundance and enriched tracer sources 
		}
		The contribution of the tracer source in the mixing model is then calculated neglecting all fractionation processes and the induced error (i.e. the difference between the calculated contribution and the prop.T input) is calculated.
}

