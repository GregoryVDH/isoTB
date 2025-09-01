Mix_Source <- function(Source1, Source2, prop.source1, by.size = FALSE) {
  if (by.size) {
    mix <- Source1
    mix$size <- Source1$size + Source2$size
    mix$F <- (Source1$F * Source1$size + Source2$F * Source2$size) / mix$size
    mix$R <- mix$F / mix$F[isoTB::isoTab[[mix$element]]$denomIso == mix$masses]
    mix$delta <- (mix$R / isoTB::isoTab[[mix$element]]$Rref - 1) * 1000
  } else {
    mix <- Source1
    mix$F <- Source1$F * prop.source1 + Source2$F * (1 - prop.source1)
    mix$R <- mix$F / mix$F[isoTB::isoTab[[mix$element]]$denomIso == mix$masses]
    mix$delta <- (mix$R / isoTB::isoTab[[mix$element]]$Rref - 1) * 1000
  }
  return(mix)
}


Mix_2Fluxes <- function(Source1, Source2, Flux1, Flux2) {
  mix <- Source1
  mix$size <- Flux1 + Flux2
  mix$F <- (Source1$F * Flux1 + Source2$F * Flux2) / mix$size
  mix$R <- mix$F / mix$F[isoTB::isoTab[[mix$element]]$denomIso == mix$masses]
  mix$delta <- (mix$R / isoTB::isoTab[[mix$element]]$Rref - 1) * 1000
  return(mix)
}


error_calc <- function(prop.T, Source_N, Source_T, mass, alpha, frac.type, reac.par) {
  # Source A fractionation
  if (frac.type[1] == "TEIF") {
    Source_N_f <- isoTB::TEIF(source = Source_N, mass = mass, epsilon = (alpha[1] - 1) * 1000, fraction = "p", Vr = 1, K = reac.par[1])
  } else {
    Source_N_f <- isoTB::FOKIF(source = Source_N, mass = mass, epsilon = (alpha[1] - 1) * 1000, fraction = "p", f = reac.par[1])
  }

  # Source B fractionation
  if (frac.type[2] == "TEIF") {
    Source_T_f <- isoTB::TEIF(source = Source_T, mass = mass, epsilon = (alpha[2] - 1) * 1000, fraction = "p", Vr = 1, K = reac.par[2])
  } else {
    Source_T_f <- isoTB::FOKIF(source = Source_T, mass = mass, epsilon = (alpha[2] - 1) * 1000, fraction = "p", f = reac.par[2])
  }

  # Mixing of Source A and B
  mix <- isoTB::Mix_Source(Source_T_f, Source_N_f, prop.source1 = prop.T, by.size = FALSE)

  # fractionation of the mix
  if (frac.type[3] == "TEIF") {
    mix_f <- isoTB::TEIF(source = mix, mass = mass, epsilon = (alpha[3] - 1) * 1000, fraction = "p", Vr = 1, K = reac.par[3])
  } else {
    mix_f <- isoTB::FOKIF(source = mix, mass = mass, epsilon = (alpha[3] - 1) * 1000, fraction = "p", f = reac.par[3])
  }
  retVal <- ((mix_f$F - Source_N$F) / (Source_T$F - Source_N$F) - prop.T)[mix_f$masses == mass]
  return(retVal)

}


min_Enrich <- function(er.tar, prop.T, Source_N, mass, alpha, frac.type, reac.par) {
  mat <- matrix(nrow = 2, ncol = 2)
  delta <- 10**seq(0, 7, 0.5)
  er <- array(dim = length(delta), data = NA)
  for (i in seq(1, length(delta), 1)){
    Source_T <- isoTB::create_Source("Source_T", element = "Mg", size = 5)
    Source_T <- isoTB::write_Source(Source_T, value = delta[i], mass = 26, notation = "delta")
    er[i] <- isoTB::error_calc(prop.T, Source_N, Source_T, mass, alpha, frac.type, reac.par) - er.tar
  }
  if (sum(er > 0) == length(er) || sum(er < 0) == length(er)) {
    return(NA)
  } else {
    mat[1, ] <- c(delta[er < 0][1], er[er < 0][1])
    mat[2, ] <- c(tail(delta[er > 0], n = 1), tail(er[er > 0], n = 1))

    while (abs(diff(mat[, 1])) > 0.1) {
      delta <- mean(mat[, 1])
      Source_T <- isoTB::create_Source("Source_T", element = "Mg", size = 5)
      Source_T <- isoTB::write_Source(Source_T, value = delta, mass = 26, notation = "delta")
      er <- isoTB::error_calc(prop.T, Source_N, Source_T, mass, alpha, frac.type, reac.par) - er.tar
      mat[sign(mat[, 2]) == sign(er), ] <- c(delta, er)
    }

    return(mean(mat[, 1]))
  }
}


min_Enrich_df <- function(er.tar, Source_N, mass, alpha, frac.type, reac.par, range = c(0, 1), res = 0.01) {
  df <- data.frame(Tracer.Source = seq(range[1], range[2], res), min.Enrich = NA)
  for (i in seq(1, nrow(df), 1)){
    df[i, 2] <- isoTB::min_Enrich(er.tar = 0.003, prop.T = df[i, 1], Source_N = Source_N, mass = 26, alpha = c(1, 1.01, 1.), frac.type = c("TEIF", "TEIF", "TEIF"), reac.par = c(0.001, 0.001, 0.001))
  }
  return(df)
}


Tracer_prop <- function(delta.samp, delta.nat, delta.tracer, element, mass, F.tracer) {
  if (length(delta.samp) == 1) {
    # Convert delta to F
    samp <- isoTB::create_Source("Sample", element, size = 1)
    samp <- isoTB::write_Source(source = samp, value = delta.samp, mass = mass, notation = "delta")

    nat <- isoTB::create_Source("nat", element, size = 1)
    if (delta.nat != 0) {
      nat <- isoTB::Frac_from_Delta(nat, mass, Delta = delta.nat, type = "TEIF")
    }
    tracer <- isoTB::create_Source("Tracer", element, size = 1)
    if (hasArg(F.tracer)) {
      tracer <- isoTB::write_Source(source = tracer, value = F.tracer, mass = mass, notation = "F")
    } else {
      tracer <- isoTB::write_Source(source = tracer, value = delta.tracer, mass = mass, notation = "delta")
    }

    return(((samp$F - nat$F) / (tracer$F - nat$F))[isoTB::isoTab[[element]]$masses == mass])
  } else {
    out <- array(dim = length(delta.samp), data = NA)
    # Convert delta to F
    nat <- isoTB::create_Source("nat", element, size = 1)
    if (delta.nat != 0) {
      nat <- isoTB::Frac_from_Delta(nat, mass, Delta = delta.nat, type = "TEIF")
    }
    if (hasArg(F.tracer)) {
      tracer <- isoTB::write_Source(source = tracer, value = F.tracer, mass = mass, notation = "F")
    } else {
      tracer <- isoTB::write_Source(source = tracer, value = delta.tracer, mass = mass, notation = "delta")
    }
    for (i in seq(1, length(delta.samp), 1)) {
      samp <- isoTB::create_Source("Sample", element, size = 1)
      samp <- isoTB::write_Source(source = samp, value = delta.samp[i], mass = mass, notation = "delta")
      out[i] <- ((samp$F - nat$F) / (tracer$F - nat$F))[isoTB::isoTab[[element]]$masses == mass]
    }
    return(out)
  }

}