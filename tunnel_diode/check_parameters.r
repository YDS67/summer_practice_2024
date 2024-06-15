#==============
# CONSTANTS
#==============
k_boltzmann <- 8.62e-5
e_h <- 2.43e-4
h2_m0 <- 0.0762
e2_4pieps0 <- 1.44
#==============
# INITIAL PARAMETERS
#==============
band_gap <- 1.12
eff_mass_e <- 0.19
eff_mass_h <- 0.49
dielectric <- 11.7
temperature <- 300
donor_conc <- 0.1
accept_conc <- 0.2
#==============
# SECONDARY PARAMETERS
#==============
n_c <- 2*(eff_mass_e*k_boltzmann*temperature/2/pi/h2_m0)^(3/2)
n_v <- 2*(eff_mass_h*k_boltzmann*temperature/2/pi/h2_m0)^(3/2)
fermi_n <- k_boltzmann*temperature*log(donor_conc/n_c)
fermi_p <- k_boltzmann*temperature*log(n_v/accept_conc)
delta_phi <- fermi_n - fermi_p + band_gap
delta_phi


