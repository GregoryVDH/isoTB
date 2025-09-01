\name{R_to_F}
\alias{R_to_F}
\title{Convert isotope ratio to fractional abundance}
\usage{
R_to_F(R, Rsum)
}
\arguments{
\item{R}{The isotope ratio \eqn{{}^{m_2/m_1}X} to be converted to fractional abundance}
\item{Rsum}{The sum of all stable isotope ratios of the considered element}
}
\value{
Fractional abundance of the stable isotope
}
\description{
Converts isotope ratios into fractional abundances. For a chemical element with n different stable isotopes, the fractional abundance of isotope i is calculated as follows:
\deqn{F = \frac{{}^{i}X/{}^{1}X}{\sum_{j=1}^{n} ({}^{j}X/{}^{1}X)}} 
where the isotope numbered 1 is the reference isotope in isotope ratios (most commonly the lighter and most abundant isotope).
}
\examples{
# Convert Magnesium 26/24 isotope ratio to fractional abundance
# by setting manually the reference isotope ratio:
R_24_24_Mg = 1
R_25_24_Mg = 0.126742
R_26_24_Mg = 0.1394
R_to_F(R = R_26_24_Mg, Rsum = sum(c(R_24_24_Mg, R_25_24_Mg, R_26_24_Mg)))
}
