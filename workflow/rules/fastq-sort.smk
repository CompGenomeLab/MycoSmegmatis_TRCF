
rule fastq_sort:
    input:
        "resources/samples/{sample}.fastq",
    output:
        "results/{sample}/{sample}_sorted.fastq"
    params:
        option='-s',
    log:
        "logs/rule/analysis/{sample}/{sample}_fastqSort.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_fastqSort.benchmark.txt",
    resources:
        memory="32GB",
        cpu=1
    conda:
        "../envs/fastq-sort.yaml"
    shell:  
        re.sub(' +', ' ', 
        """
        (echo "`date -R`: Sorting fastq file..." &&
        fastq-sort \
        --temporary-directory tmp \
        --seq \
        {input} \
        > {output} &&
        echo "`date -R`: Success! Fastq-sorting is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1
        """
        )
