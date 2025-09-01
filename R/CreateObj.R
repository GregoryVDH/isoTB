create_Source <- function(name, element, size) {
	if (sum(element == names(isoTB::isoTab)) > 0) {
		# Element name provided in argument list exists in isoTB::isoTab
		obj <- list(name = name,
								element = element,
								size = size,
								masses = isoTB::isoTab[[element]]$masses,
								F = array(dim = isoTB::isoTab[[element]]$n, data = 0),
								R = isoTB::isoTab[[element]]$Rref,
								delta = array(dim = isoTB::isoTab[[element]]$n, data = 0)
								)
		obj$delta <- (obj$R / isoTB::isoTab[[element]]$Rref - 1) * 1000
		obj$F <- obj$R / sum(obj$R)
		return(obj)
	} else {
		print("ERROR: invalid element argument")
		return()
	}
}

write_Source <- function(source, value, mass, notation) {
  if (sum(notation == c("F", "R", "delta")) != 1) {
    print("ERROR: invalid notation argument")
    return()
  }
  if (sum(mass == source$masses) != 1) {
    print("ERROR: invalid mass argument")
    return()
  }
  if (notation == "F" && value > 1) {
    print("ERROR: invalid value argument")
    return()
  }
  if (notation == "F" && value < 0) {
    print("ERROR: invalid value argument")
    return()
  }
  if (notation == "R" && value < 0) {
    print("ERROR: invalid value argument")
    return()
  }
  if (notation == "delta" && value < -1000) {
    print("ERROR: invalid value argument")
    return()
  }

  if (notation == "F") {
    source$F[source$masses == mass] <- value
    source$R <- source$F / source$F[isoTB::isoTab[[source$element]]$denomIso == source$masses]
    source$delta <- (source$R / isoTB::isoTab[[source$element]]$Rref - 1) * 1000
  } else if (notation == "R") {
    source$R[source$masses == mass] <- value
    source$F <- source$R / sum(source$R)
    source$delta <- (source$R / isoTB::isoTab[[source$element]]$Rref - 1) * 1000
  } else {
    source$delta[source$masses == mass] <- value
    source$R <- (source$delta * 1E-3 + 1) * isoTB::isoTab[[source$element]]$Rref
    source$F <- source$R / sum(source$R)
  }
  return(source)
}

print_Source <- function(source) {
	df <- data.frame(notation = c("F", "R", "delta"), m1 = NA)
	n <- isoTB::isoTab[[source$element]]$n
	for (i in seq(1, n, 1)) {
		df[, (i + 1)] <- c(source$F[i], source$R[i], source$delta[i])
	}
	colnames(df)[2:(n + 1)] <- paste("mass_", source$masses, sep = "")
	base::cat(paste("Source Name: ", source$name, "\n", sep =""))
	base::cat(paste("Element: ", source$element, "\n",sep =""))
	base::cat(paste("Size: ", source$size, "\n", sep =""))
	base::cat(paste("Isotopic composition: ", "\n", sep =""))
	base::print(df)
	return()
}

read_Source <- function(source, mass, notation) {
	if (sum(notation == c("F", "R", "delta")) == 0) {
		base::print("ERROR invalid notation argument")
		return()
	} else if (sum(mass == source$masses) == 0) {
		base::print("ERROR invalid mass argument")
		return()
	} else {
		return(source[[notation]][source$masses == mass])
	}
}