basedir=`pwd`

conda env create -f environment.yaml
cd $basedir/tools/fastix
cargo install fastix --root ./