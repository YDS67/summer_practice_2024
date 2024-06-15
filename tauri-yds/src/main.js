const { invoke } = window.__TAURI__.tauri;

let initialParameters = new Array();
let inputElements = new Array();
let inputLabels = new Array();

let setupSent;
let result;
let resultEl;
let paramsNum = 10;
let dummyString;

async function setup() {
  setupSent = await invoke("setup");
  for (let i = 0; i < paramsNum; i++) {
    inputElements[i].value = setupSent.values[i];
    inputLabels[i].innerHTML = setupSent.labels[i];
  }

  resultEl.innerHTML = "<p>Default values loaded</p>";
}

async function spectrum() {
  for (let i = 0; i < paramsNum; i++) {
    initialParameters[i] = inputElements[i].value;
  }
  result = await invoke("calculate", {
    dataFromUser: initialParameters,
  });

  resultEl.innerHTML = result.s;
}

window.addEventListener("DOMContentLoaded", () => {

  for (let i = 1; i < paramsNum+1; i++) {
    dummyString = "#input-" + i;
    inputElements[i-1] = document.querySelector(dummyString);
    dummyString = "#label-" + i;
    inputLabels[i-1] = document.querySelector(dummyString);
  }
  resultEl = document.querySelector("#result");

  setup();

  document.querySelector("#input-form").addEventListener("submit", (e) => {
    e.preventDefault();
    spectrum();
  });
});
