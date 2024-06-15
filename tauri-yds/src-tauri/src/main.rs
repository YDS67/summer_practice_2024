// Prevents additional console window on Windows in release, DO NOT REMOVE!!
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

mod consts;
mod state;
mod integrate;
mod plot;
mod file;

use std::path::PathBuf;
use serde::{Deserialize, Serialize};
use state::State;

//use tauri::Manager;

#[derive(Serialize, Deserialize, Debug)]
struct InitialMessage {
    num: usize,
    values: [f64; 10],
    labels: [String; 10],
}

#[derive(Serialize, Deserialize, Debug)]
struct FinalMessage {
    s: String,
}

#[tauri::command]
fn setup() -> Result<InitialMessage, String> {
    let input_values = [
        1.12,
        0.19,
        0.49,
        11.7,
        300.0,
        0.1,
        0.3,
        -0.5,
        1.1,
        500.0,
    ];
    let input_labels = [
        r"\( E_g, ~eV \)".to_string(),
        r"\( m_e, ~m_0 \)".to_string(),
        r"\( m_h, ~m_0 \)".to_string(),
        r"\( \varepsilon \)".to_string(),
        r"\( T, ~K \)".to_string(),
        r"\( N_d, ~nm^{-3} \)".to_string(),
        r"\( N_a, ~nm^{-3} \)".to_string(),
        r"\( V_{min}, ~V \)".to_string(),
        r"\( V_{max}, ~V \)".to_string(),
        r"\( N_{points} \)".to_string(),
    ];
    Ok(InitialMessage {
        num: 10,
        values: input_values,
        labels: input_labels,
    })
}

#[tauri::command]
fn calculate(
    data_from_user: Vec<&str>,
) -> Result<FinalMessage, String> {

    let mut parameters: Vec<f64> = vec![0.0; data_from_user.len()];

    for jp in 0..data_from_user.len() {
        parameters[jp] = data_from_user[jp].parse::<f64>().map_err(|err| err.to_string())?;
    }

    let mut state = State::fill_default();
    state.fill_initial(&parameters);
    state.fill_secondary();
    state.fill_voltage();
    state.fill_diode_current();
    state.fill_tunnel_current();
    state.fill_total_current();

    let x = state.calculation_resuts.voltage.clone();
    let y1 = state.calculation_resuts.diode_current.clone();
    let y2 = state.calculation_resuts.tunnel_current.clone();
    let y3 = state.calculation_resuts.total_current.clone();

    let flnm = format!("Tunnel_Diode_IV");
    let title = format!("Tunnel Diode Model");
    //set parameters
    let mut plot_par = plot::PlotPar::new(
        1600, 1080,
        "Voltage, V", 
        "Current density, mA/cm<sup>2</sup>", 
        &title, 
        &flnm,
        vec![
            format!("Diode current"), 
            format!("Tunneling current"),
            format!("Total current"),
        ],
    );
    plot_par.legend_al = plot::LegendAl::TopLeft;
    plot_par.colors[0] = plot::COLORS[1];
    plot_par.colors[1] = plot::COLORS[2];
    plot_par.colors[2] = plot::COLORS[0];
    plot_par.dashes[2] = plot::DASHTYPES[1].clone();
    file::save_columns_to_file(&vec![x.clone(), y1.clone(), y2.clone(), y3.clone()], &format!("{}.dat", &plot_par.flnm));
    plot::line_plot(&vec![x; 3], &vec![y1, y2, y3], &plot_par);
    let mut s;
    let current_dir_result = std::env::current_dir();
    let current_dir;
    match current_dir_result {
        Ok(path) => {
            current_dir = path.display().to_string();
            s = format!(r"<p>Calculation complete</p> <p>Saved files {}.pdf, {}.png and {}.dat</p> <p>Current directory is: {}</p>", 
                &plot_par.flnm, &plot_par.flnm, &plot_par.flnm, current_dir);
        },
        Err(_) => {
            current_dir = ".".to_string();
            s = "<p>Error: Can't find the current directory</p>".to_string()
        },
    }
    let path: PathBuf = [current_dir, format!("{}.png", &plot_par.flnm)].iter().collect();
    let open_result = open::that(path);
    match open_result {
        Ok(_) => {},
        Err(_) => {s = "<p>Error: Can't open the .png file</p>".to_string()},
    }
    Ok(FinalMessage { s })
}

fn main() {
    tauri::Builder::default()
    .invoke_handler(tauri::generate_handler![setup, calculate])
    .run(tauri::generate_context!())
    .expect("error while running tauri application");
}

