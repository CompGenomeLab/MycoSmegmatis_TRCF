
rule reposition:
    input:
        plus="results/{sample}/{sample}_sorted_plus.bed",
        minus="results/{sample}/{sample}_sorted_minus.bed",
        index="resources/ref_genomes/{build}/genome.fa.fai",
    output:
        plus=temp("results/{sample}/{sample}_sorted_plus_10.bed"),
        minus=temp("results/{sample}/{sample}_sorted_minus_10.bed"), 
    log:
        "logs/rule/analysis/{sample}/{sample}_reposition.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_reposition.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    conda:
        "../envs/bedtools.yaml"
    shell:  
        """
        (echo "`date -R`: Centering {input.plus} at damage site..." &&
        bedtools flank \
        -i {input.plus} \
        -g {input.index} \
        -l 6 \
        -r 0 |& 
        bedtools slop \
        -g {input.index} \
        -l 0 \
        -r 4 |& 
        awk '{{ if ($3-$2 == 10) {{ print }} }}' \
        > {output.plus} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.plus}; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Centering {input.minus} at damage site..." &&
        bedtools flank \
        -i {input.minus} \
        -g {input.index} \
        -l 0 \
        -r 6 |& 
        bedtools slop \
        -g {input.index} \
        -l 4 \
        -r 0 |& 
        awk '{{ if ($3-$2 == 10) {{ print }} }}' \
        > {output.minus} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.minus}; exit 1; }}  ) >> {log} 2>&1
        """