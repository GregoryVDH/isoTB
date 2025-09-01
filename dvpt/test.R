Source_A <- Create_Source("Source_A", element = "Mg", size = 10)
Source_B <- Create_Source("Source_B", element = "Mg", size = 5)
Source_B <- write_Source(Source_B, value = 100, mass = 26, notation = "delta")

Source_B <- TEIF(source = Source_B, mass = 26, epsilon = 5, fraction = "p", Vr = 1, K = 0.001)

Mix_Source(Source_B, Source_A, element = "Mg", 0.5, TRUE)



Source_N <- Create_Source("Source_A", element = "Mg", size = 10)
Source_T <- Create_Source("Source_B", element = "Mg", size = 5)
Source_T <- write_Source(source = Source_B, value = 100, mass = 26, notation = "delta")

error_calc(prop.T = 0.5, Source_N = Source_N, Source_T = Source_T, mass = 26, alpha = c(1, 1.01, 1), frac.type = c("TEIF", "TEIF", "TEIF"), reac.par = c(0.001, 0.001, 0.001))

prop.T <- 0.5
mass <- 26
alpha <- c(1, 1.01, 1)
frac.type <- c("TEIF", "TEIF", "TEIF")
reac.par <- c(0.001, 0.001, 0.001)


Source_N <- Create_Source("Source_A", element = "Mg", size = 10)
min_Enrich(er.tar = 0.003, prop.T = 0.2, Source_N = Source_N, mass = 26, alpha = c(1, 1.01, 1.01), frac.type = c("TEIF", "TEIF", "TEIF"), reac.par = c(0.001, 0.001, 0.001))

df <- data.frame(x = seq(0.01, 0.99, 0.01), delta = NA)
for (i in 1:99){
  df[i, 2] <- min_Enrich(er.tar = 0.003, prop.T = df[i, 1], Source_N = Source_N, mass = 26, alpha = c(1, 1.01, 1.), frac.type = c("TEIF", "TEIF", "TEIF"), reac.par = c(0.001, 0.001, 0.001))
}
plot(df$x, log10(df$delta), type = "l")

df <- min_Enrich_df(er.tar = 0.01, Source_N = Source_N, mass = 26, alpha = c(1, 1.01, 1.0), frac.type = c("TEIF", "TEIF", "TEIF"), reac.par = c(0.001, 0.001, 0.001), res = 0.001)

points(df$Tracer.Source, log10(df$min.Enrich))

source <- Create_Source("S", element = "Mg", size = 10)


if (TRUE) {
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
}
cat("Hello")

S = Create_Source("S", element = "Mg", size = 2.5)
S = TEIF(source = S, mass = 26, epsilon = 5, fraction = "p", Vr = 1, K = 0.001)
mass = 26
notation = "delta"

source=S
if (sum(notation == c("F", "R", "delta")) == 0) {
  base::print("ERROR invalid notation argument")

} else if (sum(mass == source$masses) == 0) {
  base::print("ERROR invalid mass argument")

} else {
  source[[notation]][source$masses == mass]
}
