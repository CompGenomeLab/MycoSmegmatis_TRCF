

rule nucleotide_table_xr_sim:
    input:
        "results/{samples}/simulation/{samples}.fa",
    output:
        dinuc=temp("results/{samples}/simulation/{samples}_dinucleotideTable.txt"),
        nuc=temp("results/{samples}/simulation/{samples}_nucleotideTable.txt"),
        filt=temp("results/{samples}/simulation/{samples}_filt.fa"),
    log:
        "logs/rule/analysis/{samples}/simulation/{samples}_nucleotide_table.log",
    benchmark:
        "logs/rule/analysis/{samples}/simulation/{samples}_nucleotide_table.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    shell:  
        """
        len="$(awk 'BEGIN{{RS=">";ORS=""}}{{print length($2)"\\n"}}' {input} | sort | uniq -c | sort -k1,1n | tail -1 | awk '{{print $2}}')"

        (echo "`date -R`: Filtering by most occurred read length..." &&
        awk -v num="$len" \
        'BEGIN{{RS=">";ORS=""}}length($2)==num{{print ">"$0}}' \
        {input} > {output.filt} &&
        echo "`date -R`: Success!" ||
        echo "`date -R`: Filtering is failed...") \
        > {log} 2>&1
        
        (echo "`date -R`: Calculating dinucleotide abundance table..." &&
        workflow/scripts/fa2kmerAbundanceTable.py \
        -i {output.filt} \
        -k 2 \
        -o {output.dinuc} &&
        echo "`date -R`: Success! Dinucleotide abundance table is calculated." ||
        echo "`date -R`: Dinucleotide abundace table cannot be calculated...") \
        >> {log} 2>&1

        (echo "`date -R`: Calculating nucleotide abundance table..." &&
        workflow/scripts/fa2kmerAbundanceTable.py \
        -i {output.filt} \
        -k 1 \
        -o {output.nuc}  &&
        echo "`date -R`: Success! Nucleotide abundance table is calculated." ||
        echo "`date -R`: Nucleotide abundace table cannot be calculated...") \
        >> {log} 2>&1
        """