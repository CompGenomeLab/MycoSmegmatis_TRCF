
rule removeDuplicatesAtFastq:
    input:
        rules.fastq_sort.output,
    output:
        "results/{sample}/{sample}_removedDup.fastq"
    log:
        "logs/rule/analysis/{sample}/{sample}_removedDup.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}removedDup.benchmark.txt",
    resources:
        memory="32GB",
        cpu=1
    conda:
        "../envs/removeDuplicatesAtFastq.yaml"
    shell:  
        re.sub(' +', ' ',
        """
        (echo "`date -R`: Removing duplicates of the sorted fastq file..." &&
        python3 workflow/scripts/removeDuplicatesAtFastq.py \
        -i {input} \
        > {output} &&
        echo "`date -R`: Success! Removing duplicates is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1
        """
        )
