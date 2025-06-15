#!/bin/bash
set -e

arg1="${1:-}"

basedir="$(pwd)"

chmod +x "$basedir/scripts/pantax" "$basedir/scripts/pantax_utils" "$basedir/scripts/data_preprocessing"

if [ -z "$arg1" ]; then
    # conda env create -f environment.yaml -y
    cd $basedir/tools/fastix
    cargo install fastix --root ./
else
    cd "$basedir/pantaxr"
    rm -f Cargo.lock
    cargo clean
    RUSTFLAGS="-C link-args=-Wl,-rpath,${CONDA_PREFIX}/lib" cargo build -r
    cp target/release/pantaxr "$basedir/scripts/pantaxr"
fi
