
rule length_mode:
    input:
        bed= rules.sort_filter.input,
        ld=  rules.sort_filter.output,
    output:
        temp("results/{sample}/{sample}_lengthMode.bed"),
    log:
        "logs/rule/analysis/{sample}/{sample}_length_mode.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_length_mode.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    shell:  
        """
        length="$(awk -v m=0 '{{if(m<$2){{m=$2;l=$1}}}}END{{print l}}' \
        {input.ld})" 

        (echo "`date -R`: Filtering the reads by the lengths..." &&
        awk -v num="$length" '{{ if ($3-$2 == num) {{ print }} }}' {input.bed} \
        > {output} &&
        echo "`date -R`: Success! Reads are filtered." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) > {log} 2>&1
        """