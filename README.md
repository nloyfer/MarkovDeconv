# MarkovDeconv

# Installation
```bash
git clone https://github.com/nloyfer/MarkovDeconv.git
cd ./MarkovDeconv/counter/
make
cd ..
```

Choose a marker file, a group file, and a corresponding reference data directory (`pat.gz*` files), e.g
```bash
groups.csv
markers.bed
./reference_data/
```

### train:
```bash
python train.py markers.bed -g groups.csv -f -v -o ./my_workind_dir
```

### deconvolve:
```bash
python deconvolve.py hep.markers -v --target Liver-Hep --pats /path/to/test/files/*pat.gz -o ./my_workind_dir
```

