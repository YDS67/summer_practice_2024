use crate::consts::*;
use crate::integrate::eval_integral;

pub struct State {
    pub initial_parameters: InitialParameters,
    pub secondary_parameters: SecondaryParameters,
    pub variable_data: VariableData,
    pub calculation_resuts: CalculationResults,
}

impl State {
    pub fn fill_default() -> State {
        State {
            initial_parameters: InitialParameters {
                band_gap: 1.12,
                eff_mass_e: 0.19,
                eff_mass_h: 0.49,
                dielectric_eps: 11.7,
                temperature: 300.0,
                donor_conc: 0.2,
                accept_conc: 0.3,
                voltage_min: -0.4,
                voltage_max: 1.1,
                voltage_points: 100,
                plot_max: 0.0,
            },
            secondary_parameters: SecondaryParameters {
                n_c: 0.0,
                n_v: 0.0,
                fermi_n: 0.0,
                fermi_p: 0.0,
                delta_phi: 0.0,
                transmission_parameter: 0.0,
                richardson_constant: 0.0,
                w_0: 0.0,
                s_0: 0.0,
            },
            variable_data: VariableData {
                voltage: 0.0,
                lower_limit_a: 0.0,
                upper_limit_b: 0.0,
            },
            calculation_resuts: CalculationResults {
                voltage: vec![0.0; 100],
                total_current: vec![0.0; 100],
                diode_current: vec![0.0; 100],
                tunnel_current: vec![0.0; 100],
            },
        }
    }
    pub fn fill_initial(&mut self, parameters: &Vec<f64>) {
        self.initial_parameters.band_gap = parameters[0];
        self.initial_parameters.eff_mass_e = parameters[1];
        self.initial_parameters.eff_mass_h = parameters[2];
        self.initial_parameters.dielectric_eps = parameters[3];
        self.initial_parameters.temperature = parameters[4];
        self.initial_parameters.donor_conc = parameters[5];
        self.initial_parameters.accept_conc = parameters[6];
        self.initial_parameters.voltage_min = parameters[7];
        self.initial_parameters.voltage_max = parameters[8];
        self.initial_parameters.voltage_points = parameters[9] as usize;
        self.initial_parameters.plot_max = parameters[10];
        self.calculation_resuts = CalculationResults {
            voltage: vec![0.0; self.initial_parameters.voltage_points],
            total_current: vec![0.0; self.initial_parameters.voltage_points],
            diode_current: vec![0.0; self.initial_parameters.voltage_points],
            tunnel_current: vec![0.0; self.initial_parameters.voltage_points],
        };
    }
    pub fn fill_secondary(&mut self) {
        self.secondary_parameters.n_c = 2.0*(
            self.initial_parameters.eff_mass_e*K_BOLTZMANN*self.initial_parameters.temperature
            / 2.0 / PI / H2_M0
        ).powf(1.5);
        self.secondary_parameters.n_v = 2.0*(
            self.initial_parameters.eff_mass_h*K_BOLTZMANN*self.initial_parameters.temperature
            / 2.0 / PI / H2_M0
        ).powf(1.5);
        self.secondary_parameters.fermi_n = K_BOLTZMANN*self.initial_parameters.temperature
            * (self.initial_parameters.donor_conc / self.secondary_parameters.n_c).ln();
        self.secondary_parameters.fermi_p = - K_BOLTZMANN*self.initial_parameters.temperature
            * (self.initial_parameters.accept_conc / self.secondary_parameters.n_v).ln();
        self.secondary_parameters.delta_phi = self.secondary_parameters.fermi_n
            - self.secondary_parameters.fermi_p + self.initial_parameters.band_gap;
        self.secondary_parameters.transmission_parameter = PI*E2_4PIEPS0*H2_M0*K_BOLTZMANN
            / self.initial_parameters.dielectric_eps
            * self.initial_parameters.temperature * self.initial_parameters.donor_conc
            * self.initial_parameters.accept_conc / self.initial_parameters.eff_mass_e
            / (self.initial_parameters.donor_conc + self.initial_parameters.accept_conc);
        // Richardson's constant for kA/cm^2
        self.secondary_parameters.richardson_constant = 1e11*E_H*self.initial_parameters.eff_mass_e
            * (K_BOLTZMANN * self.initial_parameters.temperature).powi(2)
            / 2.0 / PI.powi(2) / H2_M0;
        self.secondary_parameters.w_0 = (
            self.secondary_parameters.fermi_n
            /K_BOLTZMANN/self.initial_parameters.temperature
        ).exp();
        self.secondary_parameters.s_0 = (
            (self.secondary_parameters.fermi_n - self.secondary_parameters.delta_phi)
            /K_BOLTZMANN/self.initial_parameters.temperature
        ).exp();
    }
    pub fn fill_voltage(&mut self) {
        self.calculation_resuts.voltage = (0..self.initial_parameters.voltage_points)
            .map(|n| 
                self.initial_parameters.voltage_min + (n as f64)
                *(self.initial_parameters.voltage_max - self.initial_parameters.voltage_min)
                /(self.initial_parameters.voltage_points as f64)
            ).collect()
    }
    pub fn fill_diode_current(&mut self) {
        self.calculation_resuts.diode_current = self.calculation_resuts.voltage.clone()
            .into_iter().map(|v| 
                self.secondary_parameters.richardson_constant * self.secondary_parameters.s_0
                * ((v / K_BOLTZMANN / self.initial_parameters.temperature).exp() - 1.0)
            ).collect()
    }
    // Tunnel current calculations
    fn find_integration_limits(&mut self) {
        self.variable_data.upper_limit_b = (self.secondary_parameters.delta_phi - self.initial_parameters.band_gap - self.variable_data.voltage)
            / K_BOLTZMANN / self.initial_parameters.temperature;
        if self.variable_data.upper_limit_b < 0.0 {
            self.variable_data.upper_limit_b = 0.0
        }
        //self.variable_data.lower_limit_a = 0.0
        self.variable_data.lower_limit_a = (self.secondary_parameters.delta_phi - self.initial_parameters.band_gap - self.variable_data.voltage + self.secondary_parameters.fermi_p)
        / K_BOLTZMANN / self.initial_parameters.temperature - 5.0;
        if self.variable_data.lower_limit_a < 0.0 {
            self.variable_data.lower_limit_a = 0.0
        }
    }
    pub fn fill_tunnel_current(&mut self) {
        for jv in 0..self.initial_parameters.voltage_points {
            self.variable_data.voltage = self.calculation_resuts.voltage[jv];
            self.find_integration_limits();
            self.calculation_resuts.tunnel_current[jv] = self.secondary_parameters.richardson_constant
                * eval_integral(&self, 
                    &function_to_integrate,
                    self.variable_data.lower_limit_a,
                    self.variable_data.upper_limit_b
                )
        }
    }
    pub fn fill_total_current(&mut self) {
        for jv in 0..self.initial_parameters.voltage_points {
            self.calculation_resuts.total_current[jv] = 
            self.calculation_resuts.diode_current[jv] + self.calculation_resuts.tunnel_current[jv]
        }
    }
}

pub struct InitialParameters {
    pub band_gap: f64,
    pub eff_mass_e: f64,
    pub eff_mass_h: f64,
    pub dielectric_eps: f64,
    pub temperature: f64,
    pub donor_conc: f64,
    pub accept_conc: f64,
    pub voltage_min: f64,
    pub voltage_max: f64,
    pub voltage_points: usize,
    pub plot_max: f64,
}

pub struct SecondaryParameters {
    pub n_c: f64,
    pub n_v: f64,
    pub fermi_n: f64,
    pub fermi_p: f64,
    pub delta_phi: f64,
    pub transmission_parameter: f64,
    pub richardson_constant: f64,
    pub w_0: f64,
    pub s_0: f64,
}

pub struct VariableData {
    pub voltage: f64,
    pub lower_limit_a: f64,
    pub upper_limit_b: f64,
}

pub struct CalculationResults {
    pub voltage: Vec<f64>,
    pub total_current: Vec<f64>,
    pub diode_current: Vec<f64>,
    pub tunnel_current: Vec<f64>,
}

fn function_to_integrate(state: &State, x: f64) -> f64 {
    x * (1.0 + state.secondary_parameters.w_0 * (-x).exp()).ln()
        / (
            x + (state.secondary_parameters.delta_phi - state.variable_data.voltage).powi(3)
            / state.secondary_parameters.transmission_parameter
        )
}