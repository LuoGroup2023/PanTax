#!/bin/bash
set -e

arg1=${1:-}
arg2=${2:-}
basedir=$(pwd)

chmod +x $basedir/scripts/pantax-utils $basedir/scripts/pantax-rg

# conda env create -f environment.yaml -y

cd $basedir/pantax
# rm -f Cargo.lock

arg1_lc="${arg1,,}"
arg2_lc="${arg2,,}"

if [[ -z "$arg1_lc" ]]; then
    echo "all solvers"
    RUSTFLAGS="-C link-args=-Wl,-rpath,${CONDA_PREFIX}/lib" \
        cargo build -r
    cp target/release/pantax "$basedir/scripts/pantax"

elif [[ "$arg1_lc" == "free" ]]; then
    echo "free solvers"
    RUSTFLAGS="-C link-args=-Wl,-rpath,${CONDA_PREFIX}/lib" \
        cargo build -r --features free --no-default-features --bin pantax-free
    cp target/release/pantax-free "$basedir/scripts/pantax-free"

elif [[ "$arg1_lc" == "gb" || "$arg1_lc" == "gurobi" ]]; then
    echo "solver gurobi"
    RUSTFLAGS="-C link-args=-Wl,-rpath,${CONDA_PREFIX}/lib" \
        cargo build -r --features gb --no-default-features --bin pantax-gb
    cp target/release/pantax-gb "$basedir/scripts/pantax-gb"

elif [[ ( "$arg1_lc" == "cp" || "$arg1_lc" == "cplex" ) && -n "$arg2" ]]; then
    echo "solver cplex"
    RUSTFLAGS="-C link-args=-Wl,-rpath,${CONDA_PREFIX}/lib" \
        CPLEX_PATH="$arg2" \
        cargo build -r --features cp --no-default-features --bin pantax-cp
    cp target/release/pantax-cp "$basedir/scripts/pantax-cp"

else
    echo "Unknown solver option: $arg1"
    exit 1
fi