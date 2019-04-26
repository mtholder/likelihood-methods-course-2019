#!/bin/bash

python continuous-mcmc.py 1000000 10000000 100 100 sample-data.csv > huge-window-out2.tsv
python continuous-mcmc.py 1000000 10000000 -100 100 sample-data.csv > huge-window-out2.tsv
python continuous-mcmc.py 1000000 10 100 100 sample-data.csv > out2.tsv
python continuous-mcmc.py 1000000 10 -100 100 sample-data.csv > out2.tsv
python continuous-mcmc.py 1000000 0.0001 100 100 sample-data.csv > tiny-window-out2.tsv
python continuous-mcmc.py 1000000 0.0001 -100 100 sample-data.csv > tiny-window-out2.tsv 
