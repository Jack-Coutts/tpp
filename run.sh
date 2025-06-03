echo "Running thermal proteome profiling script...."

# filtering all mode
poetry run python3 thermal_proteome_profiling/main.py -d data/3965_U2OS_spectronaut_Report.csv -o outputs/U2OS_filtered/

# single gene mode
# poetry run python3 thermal_proteome_profiling/main.py -d data/3965_U2OS_spectronaut_Report.csv -o outputs/ -g CPSF3

# gene list mode
# poetry run python3 thermal_proteome_profiling/main.py -d data/3965_U2OS_spectronaut_Report.csv -o outputs/U2OS_gene_list/ -gl data/gene_list.txt
