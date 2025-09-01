# Simple TEIF with no products present in the system
TEIF <- function(source, mass, epsilon, fraction, Vr, K) {
  if (sum(fraction == c("r", "p")) != 1) {
		base::print("ERROR: invalid fraction argument")
		return()
	}
	mref <- isoTB::isoTab[[source$element]]$denomIso
	beta <- (1 / mref - 1 / source$masses) / (1 / mref - 1 / mass)
	alpha <- (epsilon * 1e-3 + 1)**beta
	if (fraction == "p") {
		source$delta <- alpha * (1 + Vr * K) / (1 + Vr * K * alpha) * (source$delta + 1000) - 1000
	} else {
		source$delta <- (1 + Vr * K) / (1 + Vr * K * alpha) * (source$delta + 1000) - 1000
	}
	source$R <- (source$delta * 1E-3 + 1) * isoTB::isoTab[[source$element]]$Rref
	source$F <- source$R / sum(source$R)
	return(source)
}

# Simple FOKIF with no products present in the system
FOKIF <- function(source, mass, epsilon, fraction, f) {
	if (sum(fraction == c("r", "p")) != 1) {
		base::print("ERROR: invalid fraction argument")
		return()
	}
	mref <- isoTB::isoTab[[source$element]]$denomIso
	beta <- log(mref / source$masses) / log(mref / mass)
	alpha <- (epsilon * 1e-3 + 1)**beta
	if (fraction == "p") {
		source$delta <- (1 - f**alpha) / (1 - f) * (source$delta + 1000) - 1000
	} else {
		source$delta <- f**(alpha - 1) * (source$delta + 1000) - 1000
	}
	source$R <- (source$delta * 1E-3 + 1) * isoTB::isoTab[[source$element]]$Rref
	source$F <- source$R / sum(source$R)
	return(source)
}

Frac_from_Delta <- function(source, mass, Delta, frac.type) {
	alpha <- Delta / (source$delta + 1000) + 1
	if (frac.type == "TEIF") {
		source <- isoTB::TEIF(source, mass, (alpha - 1) * 1000, fraction = "p", Vr = 1, K = 1e-4)
	} else {
		source <- isoTB::FOKIF(source, mass, (alpha - 1) * 1000, fraction = "p", f = 0.999)
	}
	return(source)
}
