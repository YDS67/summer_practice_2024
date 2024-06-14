const { invoke } = window.__TAURI__.tauri;

let Eg;
let Emin;
let Emax;
let Hg;
let R;
let dR;

let result;

let EgEl;
let EminEl;
let EmaxEl;
let HgEl;
let REl;
let dREl;

let resultEl;

async function spectrum() {
  // Learn more about Tauri commands at https://tauri.app/v1/guides/features/command
  Eg = EgEl.value;
  Emin = EminEl.value;
  Emax = EmaxEl.value;
  Hg = HgEl.value;
  R = REl.value;
  dR = dREl.value;
  result = await invoke("spectrum", {
    egInp: Eg,
    eminInp: Emin,
    emaxInp: Emax,
    hgInp: Hg,
    rInp: R,
    drInp: dR,
  });

  resultEl.innerHTML = result.s;
}

window.addEventListener("DOMContentLoaded", () => {
  EgEl = document.querySelector("#input-eg");
  EmaxEl = document.querySelector("#input-emax");
  EminEl = document.querySelector("#input-emin");
  HgEl = document.querySelector("#input-hg");
  REl = document.querySelector("#input-r");
  dREl = document.querySelector("#input-dr");

  resultEl = document.querySelector("#result");

  document.querySelector("#input-form").addEventListener("submit", (e) => {
    e.preventDefault();
    spectrum();
  });
});
