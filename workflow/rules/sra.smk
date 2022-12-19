rule sra:
    output:
        "resources/samples/{sample}.fastq", 
    params:
        sra_id=getSRA_ID(config, w.sample),
        name="{sample}",
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
        """
            (echo "`date -R`: Downloading SRA file {params.sra_id}..." &&
            prefetch {params.sra_id} \
            -O resources/samples/ &&
            vdb-validate resources/samples/$srr &&
            fastq-dump \
            resources/samples/$srr \
            --outdir resources/samples/ &&
            mv {paramd.sra_id}.fastq {params.name}.fastq \
            echo "`date -R`: Download is successful!" || 
            {{ echo "`date -R`: Process failed..."; exit 1; }} ) \
            >> {log} 2>&1 
        """
