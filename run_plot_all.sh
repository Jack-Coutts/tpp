for file in data/*; do

    # skip the gene list file
    [ "$file" = "data/gene_list.txt" ] && continue
    [ "$file" = "data/gene_list_from_filtering.txt" ] && continue

    echo "Running thermal proteome profiling script on $file ...."

    # plot all mode
    poetry run python3 thermal_proteome_profiling/main.py -d $file -o outputs/filtered/

done





