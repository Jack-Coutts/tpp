for file in data/*; do

    # skip the gene list file
    [ "$file" = "data/gene_list.txt" ] && continue

    echo "Running thermal proteome profiling script on $file ...."

    # single gene mode
    # poetry run python3 thermal_proteome_profiling/main.py -d $file -o outputs/ -g NOL6

    # gene list mode
    # poetry run python3 thermal_proteome_profiling/main.py -d $file -o outputs/gene_list -gl data/gene_list.txt

    # filtering all mode
    poetry run python3 thermal_proteome_profiling/main.py -d $file -o outputs/filtered/

done





