from math import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as spi
upper_limit = 1.0
lower_limit = 0.0
k_boltzmann=8.62e-5
e_h=2.43e-4
h2_m0=0.0762
e2_4pieps0=1.44
band_gap=1.12
eff_mass_e=0.19
eff_mass_h=0.49
dielectric=11.7
temperature=300
donor_conc=0.3
accept_conc=0.3
voltage=np.linspace(-0.5, 1, 200)
n_c=2*abs(((eff_mass_e*k_boltzmann*temperature)/(2*pi*h2_m0))**1.5)
n_v=2*abs(((eff_mass_h*k_boltzmann*temperature)/(2*pi*h2_m0))**1.5)
fermi_n=-k_boltzmann*temperature*log(n_c/donor_conc)
fermi_p=k_boltzmann*temperature*log(n_v/accept_conc)
delta_phi=fermi_n-fermi_p+band_gap
transmission_parameter=pi*e2_4pieps0*h2_m0*k_boltzmann*temperature*donor_conc*accept_conc/(donor_conc+accept_conc)/eff_mass_e/dielectric
richardson_constant=(e_h*eff_mass_e*(k_boltzmann**2)*(temperature**2))/(2*(pi**2)*h2_m0)
w_0=exp((fermi_n/(k_boltzmann*temperature)))
b_upper_limit=(delta_phi - band_gap - voltage)/(k_boltzmann*temperature)
s0=exp((fermi_n - delta_phi)/(k_boltzmann*temperature))
tunnel_current = voltage
diode_current = voltage
for j in range(len(voltage)):
    diode_current[j]=richardson_constant*s0*(exp(voltage[j]/(k_boltzmann*temperature))-1)
#    def inter(r):
#        return (r*log(1+w_0*exp(-b_upper_limit[j]*r)))/(b_upper_limit[j]*r + ((delta_phi-voltage[j])**3)/transmission_parameter)
#    if b_upper_limit[j]>0:
#        tunnel_current[j]=richardson_constant*(b_upper_limit[j]**2)*spi.quad(func=inter,a=lower_limit,b=upper_limit,full_output=0)
#    else:
#        tunnel_current[j]=0.0

fig, ax = plt.subplots()
ax.plot(voltage, diode_current)

ax.set(xlabel='Напряжение, В', ylabel='Плотность тока, А/см$^2$',
       title='Туннельный ток')
ax.grid()

fig.savefig("tunnel_current.png")
plt.show()
