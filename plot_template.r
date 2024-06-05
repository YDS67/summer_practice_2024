
# PDF изображения

# Импорт данных из файла

dat <- read.table("dataset_fermi.dat")
filename <- paste("fermi_distributions_", sep="")

x_name <- "Энергия, эВ"
y_name <- "Распределение Ферми"

x <- dat$V1
z1 <- dat$V2
z2 <- dat$V3
z3 <- dat$V4

number_of_lines <- 3
nms <- c("kT = 0","kT = 0.1","kT = 0.2")
cols <- c("black", "#ff3333", "darkgreen")
ltys <- c(2, 1, 1)
pchs <- c(NA, NA, NA)

x_max <- max(x)
x_min <- min(x)
x_length <- x_max-x_min
z_max <- max(z1,z2,z3)
z_min <- min(z1,z2,z3)

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

cairo_pdf(paste(filename, ".pdf", sep = ""), width = 6, height = 4)

# par(mar = c(bottom, left, top, right))
par(mar = c(3.5, 3.5, 0.6, 0.6), mgp = c(2, 0.6, 0), lwd = lwd_main,
    bg = "#ffffff", col = col_main, col.lab = col_main, col.axis = col_main, col.main = col_main)

line_width <- 2
psz <- 1

plot(0, 0, pch = NA, # log = "y",
    xlim = c(x_min, x_max), ylim = c(z_min, z_max), 
    xlab = x_name, 
    ylab = y_name, 
    cex.axis = tick_font_size,
    cex.main = title_font_size,
    cex.lab = axis_label_size)

axis(1, col.ticks = col_main, lwd = lwd_main, lwd.ticks = lwd_main, col = col_main, labels = FALSE)
axis(2, col.ticks = col_main, lwd = lwd_main, lwd.ticks = lwd_main, col = col_main, labels = FALSE)

grid(nx = 5, ny = 5, col = "lightgray", lty = 2, lwd = 0.5*line_width, equilogs = FALSE)

lines(x, z1, lwd = line_width, col = cols[1], lty = ltys[1])
lines(x, z2, lwd = line_width, col = cols[2], lty = ltys[2])
lines(x, z3, lwd = line_width, col = cols[3], lty = ltys[3])

legend("topright", inset = 0.02, 
    legend = nms, 
    lwd = rep(line_width, number_of_lines), 
    lty = ltys,
    pch = pchs, 
    col = cols, 
    cex = legend_font_size,
    box.col = col_main
)

dev.off()