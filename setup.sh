#!/bin/bash

## create enviroment with x86 arch
create_x86_env(){
CONDA_SUBDIR=osx-64 conda env create -n mitnanex --file environment_mac_x86.yml
conda activate mitnanex
conda config --env --set subdir osx-64
}

## Load rust module with maturin
build_rust_mod(){
maturin develop -m src/utils_rs/Cargo.toml
}

create_x86_env && build_rust_mod
