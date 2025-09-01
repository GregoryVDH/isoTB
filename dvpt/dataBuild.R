isoTab <- list()
isoTab[[1]] <- list(name = "Mg", n = 3, masses = c(24, 25, 26), Rref = c(1, 0.126742, 0.13932), denomIso = 24, IS = "NIST SRM 980")
names(isoTab) <- "Mg"
isoTab[[2]] <- list(name = "Ca", n = 6, masses = c(40, 42, 43, 44, 46, 48), Rref = c(1, 0.006674, 0.001393, 0.021222869, 4.13e-5, 0.001929), denomIso = 40, IS = "NIST 915a")
names(isoTab)[2] <- "Ca"

isoTab[[3]] <- list(name = "C", n = 2, masses = c(12, 13), Rref = c(1, 1.1237 * 1e-2), denomIso = 12, IS = "VPDB")
names(isoTab)[3] <- "C"

isoTab[[4]] <- list(name = "N", n = 2, masses = c(14, 15), Rref = c(1, 3.677 * 1e-3), denomIso = 14, IS = "atmosphere N2")
names(isoTab)[4] <- "N"

isoTab[[5]] <- list(name = "H", n = 2, masses = c(1, 2), Rref = c(1, 1.5575 * 1e-4), denomIso = 1, IS = "VSMOW")
names(isoTab)[5] <- "H"

isoTab[[6]] <- list(name = "O", n = 3, masses = c(16, 17, 18), Rref = c(1, 379.9  * 1e-6, 2.0052 * 1e-3), denomIso = 16, IS = "VSMOW")
names(isoTab)[6] <- "O"

isoTab[[7]] <- list(name = "K", n = 3, masses = c(39, 40, 41), Rref = c(1, 0.0001276478, 0.07216778), denomIso = 39, IS = "SRM 985")
names(isoTab)[7] <- "K"

isoTab[[8]] <- list(name = "B", n = 2, masses = c(10, 11), Rref = c(1, 4.04362), denomIso = 10, IS = "NBS 951")
names(isoTab)[8] <- "B"

isoTab[[9]] <- list(name = "Li", n = 2, masses = c(6, 7), Rref = c(1, 0.08215), denomIso = 6, IS = "LSVEC")
names(isoTab)[9] <- "Li"

usethis::use_data(isoTab, overwrite = TRUE)