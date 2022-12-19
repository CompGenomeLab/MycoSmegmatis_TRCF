rule sra:
    params:
        name="{sample}",
        sra_id=lambda w: getSRA_ID(config, w.sample),
    output:
        "resources/samples/{sample}.fastq", 
    log:
        "logs/rule/analysis/{sample}/{sample}_sra.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_sra.benchmark.txt",
    resources:
        memory="16GB",
        cpu=6
    conda:
        "../envs/sra.yaml"
    shell:
        re.sub(' +', ' ',
        """
            (echo "`date -R`: Downloading SRA file {params.sra_id}..." &&
            prefetch {params.sra_id} \
            -O resources/samples/ &&
            vdb-validate resources/samples/{params.sra_id} &&
            fastq-dump \
            resources/samples/{params.sra_id} \
            --outdir resources/samples/ &&
            mv {params.sra_id}.fastq {params.name}.fastq &&
            echo "`date -R`: Download is successful!" || 
            {{ echo "`date -R`: Process failed..."; exit 1; }} ) \
            >> {log} 2>&1 
        """
        )