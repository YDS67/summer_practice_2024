// Prevents additional console window on Windows in release, DO NOT REMOVE!!
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

mod plot;
use std::path::PathBuf;

use num::complex::Complex;
use num::range;
use serde::{Deserialize, Serialize};

use tauri::Manager;

#[derive(Serialize, Deserialize, Debug)]
struct MyStruct {
    s: String,
}

// Learn more about Tauri commands at https://tauri.app/v1/guides/features/command
#[tauri::command]
fn spectrum(
    eg_inp: &str,
    emin_inp: &str,
    emax_inp: &str,
    hg_inp: &str,
    r_inp: &str,
    dr_inp: &str,
) -> Result<MyStruct, String> {
    let hg = hg_inp.parse::<f64>().map_err(|err| err.to_string())?;
    let r = r_inp.parse::<f64>().map_err(|err| err.to_string())?;
    let dr = dr_inp.parse::<f64>().map_err(|err| err.to_string())?;
    let eg = eg_inp.parse::<f64>().map_err(|err| err.to_string())?;
    let emin = emin_inp.parse::<f64>().map_err(|err| err.to_string())?;
    let emax = emax_inp.parse::<f64>().map_err(|err| err.to_string())?;
    let solver: AbsorptionOf2dExcitons =
        AbsorptionOf2dExcitons::new(eg, emin, emax, hg, r, dr);
    let x = solver.hw.clone();
    let y1 = solver.get_exact_spectrum_solution();
    let y2 = solver.get_time_dependent_part();
    let flnm = format!("Spectrum_2D");
    let title = format!("2D exciton spectrum");
    //set parameters
    let mut plot_par = plot::PlotPar::new(
        1600, 1080,
        "E, eV", 
        "Absorption, a.u.", 
        &title, 
        &flnm,
        vec![
            format!("Exact solution"), 
            format!("Numerical solution"),
        ],
    );
    plot_par.colors[0] = plot::COLORS[1];
    plot_par.colors[1] = plot::COLORS[2];
    plot::line_plot(&vec![x; 2], &vec![y1, y2], &plot_par);
    let current_dir = std::env::current_dir().unwrap();
    let s = format!(r"<p>Saved files {}.pdf and {}.png</p> <p>Current folder path is: {}</p>", &plot_par.flnm, &plot_par.flnm, current_dir.display());
    let path: PathBuf = [current_dir.display().to_string(), format!("{}.png", plot_par.flnm)].iter().collect();
    open::that(path).unwrap();
    Ok(MyStruct { s })
}

fn main() {
    tauri::Builder::default()
    .invoke_handler(tauri::generate_handler![spectrum])
    .run(tauri::generate_context!())
    .expect("error while running tauri application");
}

#[derive(Debug)]
struct AbsorptionOf2dExcitons {
    zm: Complex<f64>,
    h2m0: f64,
    meh: f64,
    energy_gap: f64,
    hw: Vec<f64>,
    hg: f64,
    rex: f64,
    radius_outer: f64,
    rydberg_energy: f64,
    number_of_exciton_states: i32,
    pre: f64,
    sf: f64,
    dr0: f64,
}

impl AbsorptionOf2dExcitons {
    fn get_enf(&self, j: f64, coef: f64) -> f64 {
        self.energy_gap - self.rydberg_energy / (j - coef).powi(2)
    }

    fn get_sd(&self, hwg: &Vec<Complex<f64>>) -> Vec<Complex<f64>> {
        let mut sd: Vec<Complex<f64>> = Vec::new();
        let complex_zero: Complex<f64> = Complex::new(0.0, 0.0);

        for _ in 0..hwg.len() {
            sd.push(complex_zero);
        }

        for j in 1..(self.number_of_exciton_states + 1) {
            for idx in range(0, sd.len()) {
                sd[idx] += 1.0
                    / (j as f64 - 0.5).powi(3)
                    / self.get_enf(j as f64, 0.5)
                    / ((self.get_enf(j as f64, 0.5)).powi(2) - (hwg[idx]).powi(2));
            }
        }

        for idx in 0..sd.len() {
            sd[idx] = 2.0 / self.rex.powi(2) * hwg[idx] * sd[idx];
        }

        return sd;
    }

    fn psi(&self, z: Complex<f64>) -> Complex<f64> {
        let mut k: f64 = 0.0;

        if z.re < self.zm.re {
            k = (self.zm.re - z.re).floor();
        }

        let z_k_1_value: Complex<f64> = z + k + 1.0;
        let log_value: Complex<f64> = z_k_1_value.ln();

        let val_1: Complex<f64> = 1.0 / 2.0 / z_k_1_value;
        let val_2: Complex<f64> = 1.0 / 6.0 / z_k_1_value;
        let val_3: Complex<f64> = 1.0 / 10.0 / z_k_1_value.powi(2);
        let val_4: Complex<f64> = 10.0 / 21.0 / z_k_1_value.powi(2);

        let mut idx: f64 = 0.0;
        let mut sum_value: Complex<f64> = Complex::new(0.0, 0.0);

        while idx <= k {
            sum_value += 1.0 / (z + idx);

            idx += 1.0;
        }

        let psi_value: Complex<f64> =
            log_value - val_1 * (1.0 + val_2 * (1.0 - val_3 * val_4)) - sum_value;

        return psi_value;
    }

    fn fcont(&self, hw: f64) -> Complex<f64> {
        let complex_value: Complex<f64> = Complex::new(0.0, self.hg);
        let u1: Complex<f64> =
            (self.rydberg_energy / (self.energy_gap + hw + complex_value)).sqrt();
        let u2: Complex<f64> =
            (self.rydberg_energy / (self.energy_gap - hw - complex_value)).sqrt();
        let u3: Complex<f64> = Complex::new((self.rydberg_energy / self.energy_gap).sqrt(), 0.0);

        let fcont_value: Complex<f64> =
            u1.ln() + u2.ln() - 2.0 * u3.ln() + self.psi(u1 + 0.5) + self.psi(u2 + 0.5)
                - 2.0 * self.psi(u3 + 0.5);

        return fcont_value;
    }

    fn get_alpha_ex1_2_d(&self, hwg: &Vec<Complex<f64>>, sd: &Vec<Complex<f64>>) -> Vec<f64> {
        let mut complex_value: Complex<f64>;
        let mut alpha_ex1: Vec<f64> = Vec::new();

        for idx in range(0, sd.len()) {
            complex_value = sd[idx] + self.fcont(self.hw[idx]) * self.meh / self.h2m0 / hwg[idx];

            alpha_ex1.push(self.pre * complex_value.im);
        }

        return alpha_ex1;
    }

    fn get_chit(
        &self,
        psi_r_0: Vec<Complex<f64>>,
        eh_energies: f64,
        delta_t: f64,
        dr: f64,
    ) -> Vec<f64> {
        let mut chit: Vec<f64> = Vec::new();
        let mut complex_value: Complex<f64> = Complex::new(-self.hg, 0.0);
        let mut val: Complex<f64>;
        let mut sum: Complex<f64> = Complex::new(0.0, 0.0);

        for jw in 0..self.hw.len() {
            sum.re = 0.0;
            sum.im = 0.0;

            for idx in 0..psi_r_0.len() {
                complex_value.im = self.hw[jw] - self.energy_gap - eh_energies;
                val = (complex_value * delta_t * idx as f64).exp();

                sum += val * psi_r_0[idx];
            }

            chit.push(self.pre * 4.0 / dr.powi(2) / self.hw[jw] * delta_t * sum.re);
        }

        return chit;
    }

    fn get_exact_spectrum_solution(&self) -> Vec<f64> {
        let hwg: Vec<Complex<f64>> =
            AbsorptionOf2dExcitons::add_complex_num_to_vect(&self.hw, Complex::new(0.0, self.hg));
        let sd: Vec<Complex<f64>> = self.get_sd(&hwg);
        self.get_alpha_ex1_2_d(&hwg, &sd)
    }

    fn get_time_dependent_part(&self) -> Vec<f64> {
        let radius_outer_nm: f64 = self.radius_outer * self.rex;
        let dr: f64 = self.dr0 * self.rex;
        let number_of_grid_points: usize = (radius_outer_nm / dr).floor() as usize;

        let kinetic_energy: f64 = self.h2m0 / 2.0 / self.meh / dr.powi(2);
        let coulomb_energy: f64 = self.h2m0 / self.meh / self.rex / dr;
        let coulomb_potential: Vec<f64> = (0..number_of_grid_points).map(|j| if j==0 {4.0} else {1.0 / j as f64}).collect();
        let delta_t: f64 = 0.1;
        let number_of_time_steps = (self.sf / delta_t).floor() as usize;

        // setting tmp complex val as coefficient in formulas
        let complex_zero: Complex<f64> = Complex::new(0.0, 0.0);
        let real_unit: Complex<f64> = Complex::new(1.0, 0.0);
        let imaginary_unit: Complex<f64> = Complex::new(0.0, 1.0);

        let chit: Vec<f64>;

        // Initial state
        let mut psi_0: Vec<Complex<f64>> = vec![complex_zero; number_of_grid_points];
        let mut psi_1: Vec<Complex<f64>> = vec![complex_zero; number_of_grid_points];
        let mut psi_2: Vec<Complex<f64>> = vec![complex_zero; number_of_grid_points];
        // let mut scratch: Vec<Complex<f64>> = vec![complex_zero; number_of_grid_points];
        // let mut a: Vec<Complex<f64>> = vec![complex_zero; number_of_grid_points];
        // let mut b: Vec<Complex<f64>> = vec![real_unit; number_of_grid_points];
        // let mut c: Vec<Complex<f64>> = vec![complex_zero; number_of_grid_points];
        // let mut d: Vec<Complex<f64>> = vec![real_unit; number_of_grid_points];

        psi_0[0] = real_unit;

        // for j in 1..number_of_grid_points {
        //     a[j].im = -kinetic_energy*(1.0-0.5/j as f64)*delta_t/2.0;
        //     c[j].im = -kinetic_energy*(1.0+0.5/j as f64)*delta_t/2.0;
        //     b[j].im = (2.0*kinetic_energy - coulomb_energy*coulomb_potential[j])*delta_t/2.0;
        // }
        // c[0].im = -kinetic_energy*4.0*delta_t/2.0;
        // b[0].im = (4.0*kinetic_energy - coulomb_energy*coulomb_potential[0])*delta_t/2.0;

        // Сохраняем значение волновой функции в нуле
        let mut psi_r_0: Vec<Complex<f64>> = Vec::new();
        psi_r_0.push(psi_0[0]);

        psi_1[0] = psi_0[0]
            - imaginary_unit
                * delta_t
                * ((4.0 * kinetic_energy - coulomb_energy * coulomb_potential[0]) * psi_0[0]
                    - 4.0 * kinetic_energy * psi_0[1]);
        for j in 1..(number_of_grid_points - 1) {
            psi_1[j] = psi_0[j]
                - imaginary_unit
                    * delta_t
                    * ((2.0 * kinetic_energy - coulomb_energy * coulomb_potential[j]) * psi_0[j]
                        - (1.0 - 0.5 / j as f64) * kinetic_energy * psi_0[j - 1]
                        - (1.0 + 0.5 / j as f64) * kinetic_energy * psi_0[j + 1]);
        }
        let j = number_of_grid_points - 1;
        psi_1[j] = psi_0[j]
            - imaginary_unit
                * delta_t
                * ((2.0 * kinetic_energy - coulomb_energy * coulomb_potential[j]) * psi_0[j]
                    - (1.0 - 0.5 / j as f64) * kinetic_energy * psi_0[j - 1]);

        for _ in 1..number_of_time_steps {
            psi_2[0] = psi_0[0]
                - 2.0
                    * imaginary_unit
                    * delta_t
                    * ((4.0 * kinetic_energy - coulomb_energy * coulomb_potential[0]) * psi_1[0]
                        - 4.0 * kinetic_energy * psi_1[1]);
            for j in 1..(number_of_grid_points - 1) {
                psi_2[j] = psi_0[j]
                    - 2.0
                        * imaginary_unit
                        * delta_t
                        * ((2.0 * kinetic_energy - coulomb_energy * coulomb_potential[j])
                            * psi_1[j]
                            - (1.0 - 0.5 / j as f64) * kinetic_energy * psi_1[j - 1]
                            - (1.0 + 0.5 / j as f64) * kinetic_energy * psi_1[j + 1]);
            }
            let j = number_of_grid_points - 1;
            psi_2[j] = psi_0[j]
                - 2.0
                    * imaginary_unit
                    * delta_t
                    * ((2.0 * kinetic_energy - coulomb_energy * coulomb_potential[j]) * psi_1[j]
                        - (1.0 - 0.5 / j as f64) * kinetic_energy * psi_1[j - 1]);

            psi_0 = psi_1.clone();
            psi_1 = psi_2.clone();

            // d[0] = real_unit-imaginary_unit*((4.0*kinetic_energy - coulomb_energy*coulomb_potential[0])*psi_0[0]
            //     -4.0*kinetic_energy*psi_0[1])*delta_t/2.0;
            // for j in 1..(number_of_grid_points-1) {
            //     d[j] = real_unit-imaginary_unit*((2.0*kinetic_energy - coulomb_energy*coulomb_potential[j])*psi_0[j]
            //     -(1.0-0.5/j as f64)*kinetic_energy*psi_0[j-1]
            //     -(1.0+0.5/j as f64)*kinetic_energy*psi_0[j+1])*delta_t/2.0;
            // }
            // let j = number_of_grid_points-1;
            // d[j] = real_unit-imaginary_unit*((2.0*kinetic_energy - coulomb_energy*coulomb_potential[j])*psi_0[j]
            //     -(1.0-0.5/j as f64)*kinetic_energy*psi_0[j-1])*delta_t/2.0;

            // thomas_algorithm(&mut psi_0, &mut scratch, &a, &b, &c, &mut d);

            psi_r_0.push(psi_0[0]);
        }

        // get chit
        chit = self.get_chit(psi_r_0.clone(), 0.0, delta_t, dr);

        chit
    }

    fn add_complex_num_to_vect(vec: &Vec<f64>, complex_value: Complex<f64>) -> Vec<Complex<f64>> {
        let mut result_vect: Vec<Complex<f64>> = Vec::new();

        for item in vec {
            result_vect.push(item + complex_value);
        }

        return result_vect;
    }

    fn new(
        eg: f64,
        emin: f64,
        emax: f64,
        hg: f64,
        r: f64,
        dr: f64,
    ) -> AbsorptionOf2dExcitons {
        let zm: Complex<f64> = Complex::new(16.0, 0.0);
        let h2m0: f64 = 0.07619964231;
        let e24pi0: f64 = 1.439964547;
        let al0: f64 = 1.0 / 137.0;
        let me: f64 = 0.044;
        let mh: f64 = 0.308;
        let meh: f64 = me * mh / (me + mh);
        let energy_gap: f64 = eg;
        let optical_matrix_parameter: f64 = 23.84;
        let eps: f64 = 13.9;
        let nb: f64 = eps.sqrt();
        let hg: f64 = hg;
        let rex0: f64 = h2m0 / e24pi0;
        let rex: f64 = rex0 * eps / meh;
        let radius_outer: f64 = r;
        let rydberg_energy: f64 = h2m0 / 2.0 / meh / rex.powi(2);
        let number_of_exciton_states: i32 = 100;
        let pre: f64 = 4.0 / 3.0 * al0 / nb * optical_matrix_parameter * h2m0 / rex * 1000.0;
        let sf: f64 = 7.0 / hg;
        let dr0: f64 = dr;
        let mut jw = 0.0;

        AbsorptionOf2dExcitons {
            zm,
            h2m0,
            meh,
            energy_gap,
            hw: [0.0; 2000]
                .into_iter()
                .map(|_| {
                    jw += 1.0;
                    emin + jw / 2000. * (emax - emin)
                })
                .collect(),
            hg,
            rex,
            radius_outer,
            rydberg_energy,
            number_of_exciton_states,
            pre,
            sf,
            dr0,
        }
    }
}

fn thomas_algorithm(
    x: &mut Vec<Complex<f64>>,
    s: &mut Vec<Complex<f64>>,
    a: &Vec<Complex<f64>>,
    b: &Vec<Complex<f64>>,
    c: &Vec<Complex<f64>>,
    d: &mut Vec<Complex<f64>>,
) {
    s[0] = c[0] / b[0];
    d[0] = d[0] / b[0];
    let n = x.len();

    /* loop from 1 to X - 1 inclusive */
    for ix in 1..n {
        s[ix] = c[ix] / (b[ix] - a[ix] * s[ix - 1]);
        d[ix] = (d[ix] - a[ix] * d[ix - 1]) / (b[ix] - a[ix] * s[ix - 1]);
    }

    x[n - 1] = d[n - 1];

    /* loop from X - 2 to 0 inclusive */
    for ix in (n - 2)..=0 {
        x[ix] = d[ix] - s[ix] * x[ix + 1]
    }
}

use std::f64::consts::PI;
const TOLERANCE: f64 = 1e-16;
const H2M0: f64 = 0.0762;
const E24PI0: f64 = 1.44;
const EPS: f64 = 13.9;

fn add_vectors(vec1: &Vec<f64>, vec2: &Vec<f64>) -> Vec<f64> {
    let mut result_vect: Vec<f64> = Vec::new();

    for je in 0..vec1.len() {
        result_vect.push(vec1[je] + vec2[je]);
    }

    return result_vect;
}

const MINUS_CI: [f64; 40] = [
    0.02256066174634607,
    0.006116639131919573,
    0.002769358610857675,
    0.001568553021292536,
    0.001007172104752763,
    0.0007006894745077421,
    0.0005153576529232389,
    0.0003948533143785074,
    0.0003121365636133246,
    0.0002529199146935864,
    0.0002090794372250259,
    0.0001757198260382837,
    0.0001497490190791835,
    0.0001291362492747358,
    0.0001125032190398279,
    0.00009888784174648064,
    0.00008760207012103075,
    0.00007814331019930193,
    0.00007013752554872027,
    0.00006330170927334961,
    0.00005741854044895157,
    0.00005231890537746918,
    0.00004786961336164939,
    0.00004396461492350048,
    0.00004051862608885234,
    0.00003746243385373541,
    0.00003473939467765655,
    0.00003230279169944803,
    0.00003011381817630760,
    0.00002814002313657675,
    0.00002635410202251499,
    0.00002473294751330401,
    0.00002325689847169071,
    0.00002190914112532553,
    0.00002067522821346019,
    0.00001954269026991978,
    0.00001850071940551399,
    0.00001753991053877051,
    0.00001665204844990988,
    0.00001582993161395283,
];
const SI: [f64; 40] = [
    1.418151576132628,
    1.492161225584460,
    1.518033961467181,
    1.531131284990666,
    1.539029079577564,
    1.544307522339321,
    1.548083269679786,
    1.550917632803829,
    1.553123463303973,
    1.554888871044745,
    1.556333738245307,
    1.557538071093344,
    1.558557302244981,
    1.559431050403327,
    1.560188383041399,
    1.560851108998710,
    1.561435910729740,
    1.561955766742638,
    1.562420925373798,
    1.562839586736321,
    1.563218390009691,
    1.563562767945477,
    1.563877208964840,
    1.564165453777924,
    1.564430644845112,
    1.564675441367728,
    1.564902108739276,
    1.565112588838288,
    1.565308555784678,
    1.565491460549714,
    1.565662566935343,
    1.565822980809969,
    1.565973674030547,
    1.566115504144548,
    1.566249230715517,
    1.566375528928525,
    1.566495000989941,
    1.566608185727732,
    1.566715566715164,
    1.566817579176233,
];

const KRONROD_NODES_21: [f64; 21] = [
    -0.9956571630258080807355273,
    -0.9739065285171717200779640,
    -0.9301574913557082260012072,
    -0.8650633666889845107320967,
    -0.7808177265864168970637176,
    -0.6794095682990244062343274,
    -0.5627571346686046833390001,
    -0.4333953941292471907992659,
    -0.2943928627014601981311266,
    -0.1488743389816312108848260,
    0.0000000000000000000000000,
    0.1488743389816312108848260,
    0.2943928627014601981311266,
    0.4333953941292471907992659,
    0.5627571346686046833390001,
    0.6794095682990244062343274,
    0.7808177265864168970637176,
    0.8650633666889845107320967,
    0.9301574913557082260012072,
    0.9739065285171717200779640,
    0.9956571630258080807355273,
];

const KRONROD_WEIGHTS_21: [f64; 21] = [
    0.0116946388673718742780644,
    0.0325581623079647274788190,
    0.0547558965743519960313813,
    0.0750396748109199527670431,
    0.0931254545836976055350655,
    0.1093871588022976418992106,
    0.1234919762620658510779581,
    0.1347092173114733259280540,
    0.1427759385770600807970943,
    0.1477391049013384913748415,
    0.1494455540029169056649365,
    0.1477391049013384913748415,
    0.1427759385770600807970943,
    0.1347092173114733259280540,
    0.1234919762620658510779581,
    0.1093871588022976418992106,
    0.0931254545836976055350655,
    0.0750396748109199527670431,
    0.0547558965743519960313813,
    0.0325581623079647274788190,
    0.0116946388673718742780644,
];

const EULER_GAMMA: f64 = 0.5772156649015329;

fn electron_hole_wavefunction(z: &Vec<f64>, qw_width: f64, lambda: usize) -> Vec<f64> {
    let lambdaf = lambda as f64;
    let wavefunction = z
        .clone()
        .into_iter()
        .map(|z| (PI * lambdaf * z / qw_width).sin().powi(2) / qw_width)
        .collect();
    wavefunction
}

fn eval_potential_fd(s: f64, z: &Vec<f64>, dz: f64, wavefunction: &Vec<f64>) -> f64 {
    let nz = z.len();
    let mut potential = 0.0;
    for ne in 0..nz {
        for nh in 0..nz {
            potential += dz.powi(2) * wavefunction[ne] * wavefunction[nh]
                / ((z[ne] - z[nh]).powi(2) + s.powi(2)).sqrt()
        }
    }
    potential
}

fn eval_potential_fd_zero(ds: f64, z: &Vec<f64>, dz: f64, wavefunction: &Vec<f64>) -> f64 {
    let nz = z.len();
    let mut potential = 0.0;
    for ne in 0..nz {
        for nh in 0..nz {
            if ne == nh {
                potential += dz.powi(2) * wavefunction[ne] * wavefunction[nh] * 4.0 / ds
            } else {
                potential +=
                    dz.powi(2) * wavefunction[ne] * wavefunction[nh] / (z[ne] - z[nh]).abs()
            }
        }
    }
    potential
}

fn grid_potential_fd(
    s: &Vec<f64>,
    ds: f64,
    z: &Vec<f64>,
    dz: f64,
    qw_width: f64,
    lambda: usize,
) -> Vec<f64> {
    let wavefunction = electron_hole_wavefunction(z, qw_width, lambda);
    let potential = s
        .clone()
        .into_iter()
        .map(|s| {
            if s > TOLERANCE {
                -eval_potential_fd(s, z, dz, &wavefunction)
            } else {
                -eval_potential_fd_zero(ds, z, dz, &wavefunction)
            }
        })
        .collect();
    potential
}

// fn to_integrate(x: f64, lambda: usize, l: usize, s: f64) -> f64 {
//     let sign = if l % 2 == 0 {1.0} else {-1.0};
//     let lf = l as f64;
//     let lambdaf = lambda as f64;
//     (3.0*sign*(PI*x).sin()-PI*(2.0*lambdaf-lf-x)*(1.0-sign*(PI*x).cos()))
//         /((x+lf).powi(2)+(2.0*lambdaf*s).powi(2)).sqrt()
// }

// fn eval_integral(lambda: usize, s: f64) -> f64 {
//     let lambdaf = lambda as f64;
//     let mut int = 0.0;
//     for l in 0..(2*lambda) {
//         for n in 0..KRONROD_NODES_21.len() {
//             int += to_integrate((1.0+KRONROD_NODES_21[n])/2.0, lambda, l, s)*KRONROD_WEIGHTS_21[n]/2.0
//         }
//     }
//     int/(8.0*PI*lambdaf)
// }

// fn eval_integral_zero(lambda: usize) -> f64 {
//     let lambdaf = lambda as f64;
//     3.0*SI[lambda-1]/(8.0*PI*lambdaf) + 1.0/4.0*(1.0-EULER_GAMMA-MINUS_CI[lambda-1]-(2.0*PI*lambdaf).ln())
// }

// fn eval_potential(s: f64, integral: f64) -> f64 {
//     integral + 3.0/4.0*((1.0/s).asinh() - 1.0/(s+(1.0+s*s).sqrt()))
// }

// // This is important! If we use a different value at 0, the potential bends visibly
// fn eval_potential_zero(lambda: usize, ds: f64) -> f64 {
//     eval_integral_zero(lambda) + 3.0/4.0*((2.0/ds).asinh() + 1.0/2.0)
// }

// fn grid_potential(s: &Vec<f64>, ds: f64, lambda: usize) -> Vec<f64> {
//     let potential = s.clone().into_iter().map(|s| {
//         if s > TOLERANCE {
//             let integral = eval_integral(lambda, s);
//             -eval_potential(s, integral)
//         } else {
//             -eval_potential_zero(lambda, ds)
//         }
//     }).collect();
//     potential
// }
