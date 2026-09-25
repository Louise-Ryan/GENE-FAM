{
    echo -e "Species\tRefSeq\tAUGUSTUS\tTotal"
    for dir in *outfiles; do
        augustus=$(grep -c "^>Hit" "$dir"/*cds_prot.fa)
        total=$(grep -c "^>" "$dir"/*cds_prot.fa)
        refseq=$((total - augustus))
        species="${dir%_outfiles}"
        printf "%s\t%s\t%s\t%s\n" "$species" "$refseq" "$augustus" "$total"
    done
} | tee sequence_counts.tsv
