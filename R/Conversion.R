# Conversion.R
#
# Description:
#   contains the different functions to convert between the different isotopic notations:
#     fractional abundance
#     isotope ratio
#     delta notation
#

R_to_F <- function(R, Rsum) {
  return(R / (Rsum))
}

R_to_Delta <- function(R, element, mass, Rref) {
  if (hasArg(Rref)) {
    return((R / Rref - 1) * 1E3)
  } else {
    Rref <- isoTB::isoTab[[element]]$Rref[isoTB::isoTab[[element]]$masses == mass]
    return((R / Rref - 1) * 1E3)
  }
}

Delta_to_R <- function(delta, element, mass, Rref) {
	if (! hasArg(Rref)) {
  	Rref <- isoTB::isoTab[[element]]$Rref[isoTB::isoTab[[element]]$masses == mass]
	}
  return((delta * 1e-3 + 1) * Rref)
}

Delta_to_F <- function(delta, element, mass) {
	if (is.null(nrow(delta))) {
		R <- isoTB::Delta_to_R(delta, element, mass)
		Rsum <- sum(isoTB::isoTab[[element]]$Rref[isoTB::isoTab[[element]]$masses != mass])
		return(isoTB::R_to_F(R = R, Rsum = Rsum + R))
	} else {
		if (ncol(delta) != (isoTB::isoTab[[element]]$n - 1)) {
			base::print("ERROR: invalid delta argument")
			return()
		}
		R <- as.data.frame(matrix(nrow = nrow(delta), ncol = ncol(delta)))
		colnames(R) <- paste("R", isoTB::isoTab[[element]]$masses[2:isoTB::isoTab[[element]]$n], "_", isoTB::isoTab[[element]]$denomIso, "_", element, sep = "")
		for (i in seq(1, ncol(R), 1)) {
			R[, i] <- isoTB::Delta_to_R(delta[, 1], Rref = isoTB::isoTab[[element]]$Rref[i + 1])
		}
		Frac <- data.frame(v = 1, R)
		colnames(Frac) <- paste("F_", isoTB::isoTab[[element]]$masses, sep = "")
		Rsum <- rowSums(Frac)
		for (i in seq(1, ncol(Frac), 1)) {
			Frac[, i] <- Frac[, i] / (Rsum)
		}
		return(Frac)
	}
}