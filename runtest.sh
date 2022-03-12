python3 dump.py 1e+4 1e-10 ir_nlambda4_ndigit10.dat
python3 mk_preset.py --nlambda 1 2 3 4 --ndigit 10 > sparse_ir_preset.f90
make test
