# PLOTTING PARAMETERS, DON'T CHANGE
col_main <- "#282828"
lwd_main <- 5
axis_label_size <- 1.8
legend_font_size <- 1.5
tick_font_size <- 1.3
title_font_size <- 1.5
lwd_main <- 1
axis_label_size <- 14/12
legend_font_size <- 14/12
tick_font_size <- 1
title_font_size <- 14/12
line_width <- 2
psz <- 1

# =====================================
h2m0 <- 0.0762
me <- 0.067
h2me <- h2m0/me
K <- 3*pi^2*h2me/2
d <- 8
H <- 2

# Create an array of QD with normally distributed sizes
# Mean or average size
a0 <- 4.25
# Dispersion
b <- 0.75
# NUMBER OF PARTICLES!!!
N_particles <- 500

distrib <- function(a){exp(-(a-a0)^2/b^2)/sqrt(pi)/b}

energy_cube <- function(a){K/a^2}
gamma_value <- function(a,E){
	k <- sqrt(2*E/h2me)
	w <- (d-a)/2
	p <- H*w/h2me
	h2me/2*k^3/a/(k^2+2*p^2)
}

# Amplitude of the distribution
f0 <- distrib(a0)

# Smallest and largest size
a_min <- a0-3*b
a_max <- 6+3*b
E0_min <- energy_cube(a_max) - 5*gamma_value(a_max,energy_cube(a_max))
E0_max <- energy_cube(a_min) + 5*gamma_value(a_min,energy_cube(a_min))

# x points for plotting
Na <- 1000
an <- a_min + (a_max-a_min)*(1:Na)/Na

# Algorithm start

particles <- NULL

Np <- 0
while(Np < N_particles){
	# Two random numbers between 0 and 1
	r1 <- runif(1)
	r2 <- runif(1)
	# Pick a random point between a_min and a_max
	a <- a_min + (a_max-a_min)*r1
	f <- distrib(a)
	if(f > r2*f0){
		particles <- c(particles, a)
		Np <- Np+1
	}
}

# ENERGY DISTRIBUTION

particles_energy <- energy_cube(particles)
particles_gamma <- gamma_value(particles,particles_energy)

# =================
# PLOTTING
# =================

filename <- paste("size_distribution_Np", Np, "_a0", a0, "_b", b, sep="")
x_name <- "а, нм"
y_name <- "Распределение, усл.ед."

cairo_pdf(paste(filename, ".pdf", sep = ""), width = 6, height = 4)
# par(mar = c(bottom, left, top, right))
par(mar = c(3.5, 3.5, 0.6, 0.6), mgp = c(2, 0.6, 0), lwd = lwd_main,
    bg = "#ffffff", col = col_main, col.lab = col_main, col.axis = col_main, col.main = col_main)


hist(particles, breaks = "FD", freq = FALSE, 
	xlab = x_name, 
    ylab = y_name, 
	main = NA,
    cex.axis = tick_font_size,
    cex.main = title_font_size,
    cex.lab = axis_label_size,
	col = "burlywood2"
)

axis(1, col.ticks = col_main, lwd = lwd_main, lwd.ticks = lwd_main, col = col_main, labels = FALSE)
axis(2, col.ticks = col_main, lwd = lwd_main, lwd.ticks = lwd_main, col = col_main, labels = FALSE)

#grid(nx = 5, ny = 5, col = "lightgray", lty = 2, lwd = 0.5*line_width, equilogs = FALSE)

lines(an, distrib(an), lwd = line_width, col = "red")

dev.off()

# ================

filename <- paste("energy_distribution_Np", Np, "_a0", a0, "_b", b, sep="")
x_name <- "E0, эВ"
y_name <- "Распределение, усл.ед."

cairo_pdf(paste(filename, ".pdf", sep = ""), width = 6, height = 4)
# par(mar = c(bottom, left, top, right))
par(mar = c(3.5, 3.5, 0.6, 0.6), mgp = c(2, 0.6, 0), lwd = lwd_main,
    bg = "#ffffff", col = col_main, col.lab = col_main, col.axis = col_main, col.main = col_main)


hist(particles_energy, breaks = "FD", freq = FALSE, 
	xlab = x_name, 
    ylab = y_name, 
	main = NA,
    cex.axis = tick_font_size,
    cex.main = title_font_size,
    cex.lab = axis_label_size,
	col = "lightblue"
)

axis(1, col.ticks = col_main, lwd = lwd_main, lwd.ticks = lwd_main, col = col_main, labels = FALSE)
axis(2, col.ticks = col_main, lwd = lwd_main, lwd.ticks = lwd_main, col = col_main, labels = FALSE)

#grid(nx = 5, ny = 5, col = "lightgray", lty = 2, lwd = 0.5*line_width, equilogs = FALSE)

dev.off()

# ================

filename <- paste("gamma_distribution_Np", Np, "_a0", a0, "_b", b, sep="")
x_name <- expression(paste(gamma, ", эВ"))
y_name <- "Распределение, усл.ед."

cairo_pdf(paste(filename, ".pdf", sep = ""), width = 6, height = 4)
# par(mar = c(bottom, left, top, right))
par(mar = c(3.5, 3.5, 0.6, 0.6), mgp = c(2, 0.6, 0), lwd = lwd_main,
    bg = "#ffffff", col = col_main, col.lab = col_main, col.axis = col_main, col.main = col_main)


hist(particles_gamma, breaks = "FD", freq = FALSE, 
	xlab = x_name, 
    ylab = y_name, 
	main = NA,
    cex.axis = tick_font_size,
    cex.main = title_font_size,
    cex.lab = axis_label_size,
	col = "seagreen"
)

axis(1, col.ticks = col_main, lwd = lwd_main, lwd.ticks = lwd_main, col = col_main, labels = FALSE)
axis(2, col.ticks = col_main, lwd = lwd_main, lwd.ticks = lwd_main, col = col_main, labels = FALSE)

#grid(nx = 5, ny = 5, col = "lightgray", lty = 2, lwd = 0.5*line_width, equilogs = FALSE)

dev.off()

# ====================
# Вычисление ВАХ
# ====================

Vmin <- 0
Vmax <- 4
Npoints <- 50
Vn <- Vmin + (Vmax-Vmin)*(1:Npoints)/Npoints
In1 <- Vn
In2 <- Vn
In3 <- Vn
kT <- 0.026

fermi <- function(E){1/(1+exp(E/kT))}

density_of_states <- function(a,E0,E){
	gam <- gamma_value(a,E)
	gam^2/(gam^2 + (E-E0)^2)
}

function_to_integrate <- function(a,E0,Vg,V,E){
	density_of_states(a,E0,E-Vg)*(fermi(E-V)-fermi(E))
}

current <- function(Vg,V){
	curr <- 0
	for(np in 1:Np){
		a <- particles[np]
		E0 <- particles_energy[np]
		to_integrate <- function(E){function_to_integrate(a,E0,Vg,V,E)}
		curr <- curr + integrate(to_integrate,E0_min+Vg,E0_max+Vg)$value
	}
	curr/Np
}

for(n in 1:Npoints){
	In1[n] <- current(0,Vn[n])
	In2[n] <- current(0.5,Vn[n])
	In3[n] <- current(1,Vn[n])
}

# =================
# PLOTTING
# =================

filename <- paste("IV_curve_Np", Np, "_a0", a0, "_b", b, sep="")
x_name <- "Напряжение, В"
y_name <- "Сила тока, усл.ед."
title <- paste(N_particles, " КТ, сред. диаметр ", a0, " нм, дисп. ", b, "нм.")

number_of_lines <- 3
nms <- c("Vg = 0","Vg = 0.5","Vg = 1")
cols <- c("blue", "red", "purple")
ltys <- c(1, 1, 1)
pchs <- c(NA, NA, NA)

x <- Vn
x_min <- Vmin
x_max <- Vmax
y1 <- In1
y2 <- In2
y3 <- In3

y_min <- min(y1,y2,y3)
y_max <- max(y1,y2,y3)

cairo_pdf(paste(filename, ".pdf", sep = ""), width = 6, height = 4)
# par(mar = c(bottom, left, top, right))
par(mar = c(3.5, 3.5, 3.6, 0.6), mgp = c(2, 0.6, 0), lwd = lwd_main,
    bg = "#ffffff", col = col_main, col.lab = col_main, col.axis = col_main, col.main = col_main)


plot(0, 0, pch = NA, # log = "y",
    ylim = c(y_min, y_max), xlim = c(x_min, x_max),
    xlab = x_name, 
    ylab = y_name, 
	main = title,
    cex.axis = tick_font_size,
    cex.main = title_font_size,
    cex.lab = axis_label_size
)

axis(1, col.ticks = col_main, lwd = lwd_main, lwd.ticks = lwd_main, col = col_main, labels = FALSE)
axis(2, col.ticks = col_main, lwd = lwd_main, lwd.ticks = lwd_main, col = col_main, labels = FALSE)

grid(nx = 5, ny = 5, col = "lightgray", lty = 2, lwd = 0.5*line_width, equilogs = FALSE)

lines(x, y1, lwd = line_width, col = cols[1], lty = ltys[1])
lines(x, y2, lwd = line_width, col = cols[2], lty = ltys[2])
lines(x, y3, lwd = line_width, col = cols[3], lty = ltys[3])

legend("bottomright", inset = 0.02, 
    legend = nms, 
    lwd = rep(line_width, number_of_lines), 
    lty = ltys,
    pch = pchs, 
    col = cols, 
    cex = legend_font_size,
    box.col = col_main
)

dev.off()
