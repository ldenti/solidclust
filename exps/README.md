``` sh
snakemake -c4 -p --use-conda [-n]

# assuming $WD to be the output directory (as in config/config.yml)
python3 scripts/plot.py $WD
python3 scripts/get_accuracy_table.py $WD
python3 scripts/get_efficiency_table.py $WD
```
