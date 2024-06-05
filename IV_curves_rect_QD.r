# Моделирование ВАХ прямоугольных КТ

# Параметры, которые надо задавать

h2m0 <- 0.0762 # eV*nm^2
me <- 0.067 
w <- 1 # nm
H <- 4 # eV
Emin <- 0 # eV
Emax <- 5 # eV
Npoints <- 1000
dE <- (Emax-Emin)/Npoints
Epoints <- Emin + (1:Npoints)*dE
nmax <- 10 # сколько уровней энергии мы учитываем по каждой координате

# Температура
kT <- 0.026 # eV

fermi <- function(E){1/(1+exp(E/kT))}

# Размеры КТ (будут меняться)
a <- 10 #nm
b <- 10 #nm
cc <- 10 #nm

d <- a+2*w

# Уширение по энергии
p <- me/h2m0*H*w
K <- h2m0/2/me
k <- function(E){sqrt(E/K)}
broadening <- function(E){K*k(E)^3/((k(E)^2+2*p^2)*d)}

gamma_points <- broadening(Epoints)

par(mar = c(5, 5, 5, 5))
plot(Epoints, gamma_points, type = "l", col = "blue", xlab = "Энергия, эВ", ylab = "Уширение, эВ", lwd = 3, cex.lab = 2, cex.axis = 1.5)

write.table(data.frame(Epoints, gamma_points), "Gamma_vs_E.dat", row.names = FALSE, col.names = FALSE)

# Вычисление плотности состояний
En <- function(nx, ny, nz){pi^2*K * (nx^2/a^2 + ny^2/b^2 + nz^2/cc^2)}

Dpoints <- rep(0, Npoints)
for(nE in 1:Npoints){
	result <- 0
	for(nx in 1:nmax){
	for(ny in 1:nmax){
	for(nz in 1:nmax){
		result <- result + gamma_points[nE]^2/(gamma_points[nE]^2 + (Epoints[nE]-En(nx,ny,nz))^2)
	}
	}
	}
	Dpoints[nE] <- result
}

par(mar = c(5, 5, 5, 5))
plot(Epoints, Dpoints, type = "l", col = "purple", xlab = "Энергия, эВ", ylab = "Плотность состояний, эВ", lwd = 3, cex.lab = 2, cex.axis = 1.5)

write.table(data.frame(Epoints, Dpoints), "DensityOfStates_vs_E.dat", row.names = FALSE, col.names = FALSE)

# Окончательное выражение для ВАХ
current <- function(voltage){
	V <- voltage
	sum(Dpoints*(fermi(Epoints-V)-fermi(Epoints)))*dE
}

NV <- 200
Vmin <- -1
Vmax <- 4
dV <- (Vmax-Vmin)/NV
Vn <- Vmin + (1:NV)*dV
In <- Vn

for(nV in 1:NV){
	In[nV] <- current(Vn[nV])
}

par(mar = c(5, 5, 5, 5))
plot(Vn, In, type = "l", col = "red", xlab = "Напряжение, В", ylab = "Ток, условные ед.", lwd = 3, cex.lab = 2, cex.axis = 1.5)

write.table(data.frame(Vn, In), "Current_vs_Voltage.dat", row.names = FALSE, col.names = FALSE)
