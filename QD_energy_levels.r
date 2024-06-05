# Подсчёт числа уровней в кубической КТ
me <- 0.067
h2m0 <- 0.0762 # eV * nm^2
a <- 10 #nm
Chi <- 5 # eV
E0 <- h2m0 * pi^2 / 2 / me / a^2
E0
En <- function(nx, ny, nz){E0 * (nx^2 + ny^2 + nz^2)}
nmax <- floor(sqrt(Chi/E0 - 2))
# Надо сделать посчёт электронов для каждого значения энергии (с учётом вырождения)
Elevels <- NULL
Ne <- 0
for (jx in 1:nmax) {
for (jy in 1:nmax) {
for (jz in 1:nmax) {
	Ej <- En(jx, jy, jz)
	if (Ej < Chi) {
		Elevels <- c(Elevels, Ej)
		Ne <- Ne + 1
	}
}
}
}

Elevels <- sort(Elevels)

degenerate_levels <- 1
Elevels_new <- Elevels[1]
jlevels <- 1

for (je in 2:Ne) {
	if (Elevels[je] == Elevels[je-1]) {
		degenerate_levels[jlevels] <- degenerate_levels[jlevels] + 1
	} else {
		jlevels <- jlevels + 1
		degenerate_levels <- c(degenerate_levels, 1)
		Elevels_new <- c(Elevels_new, Elevels[je])
	}
}

plot(0, 0, col = "white", xlim = c(0, Chi), ylim = c(0, 1), xlab = "Энергия, эВ", ylab = " ")

for (je in 1:jlevels) {
	abline(v = Elevels_new[je], lwd = degenerate_levels[je], col="blue")
}

#plot(Elevels_new, degenerate_levels, pch = 16, cex = 2, col = "purple", xlab = "Энергия, эВ", ylab = "Степень вырождения")





