basedir=`pwd`

# conda env create -f environment.yaml -y
chmod +x $basedir/scripts/pantax
chmod +x $basedir/scripts/data_preprocessing
cd $basedir/tools/fastix
cargo install fastix --root ./