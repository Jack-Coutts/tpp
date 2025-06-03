echo "Running thermal proteome profiling script...."

poetry run python3 thermal_proteome_profiling/main.py -d data/3965_SJSA1_spectronaut_Report.csv -o outputs/

# poetry run python3 thermal_proteome_profiling/main.py -d data/3965_U2OS_spectronaut_Report.csv -o outputs/
