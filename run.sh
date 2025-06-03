echo "Running thermal proteome profiling script...."

# filtering all mode
# poetry run python3 thermal_proteome_profiling/main.py -d data/3965_SJSA1_spectronaut_Report.csv -o outputs/

# single gene mode
# poetry run python3 thermal_proteome_profiling/main.py -d data/3965_SJSA1_spectronaut_Report.csv -o outputs/ -g CPSF3

# gene list mode
poetry run python3 thermal_proteome_profiling/main.py -d data/3965_SJSA1_spectronaut_Report.csv -o outputs/SJSA1_gene_list/ -gl data/gene_list.txt
