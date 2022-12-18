build = config["build"]

rule genome_build:
    # input: "resources/ref_genomes/build/genome.fa"
    input: report(f"resources/ref_genomes/{build}/genome.fa", category="genome")
    # output: "resources/ref_genomes/build/Bowtie2/genome"
    params: 
        base = f"resources/ref_genomes/{build}/Bowtie2/genome"
    output:
        index_1 = f"resources/ref_genomes/{build}/Bowtie2/genome.1.bt2"
    log: f"logs/rule/analysis/{build}/log/bowtie2_build.log"
    benchmark: f"logs/rule/analysis/{build}/log/bowtie2_build.benchmark.txt"
    conda:
        "../envs/align.yaml"
    resources:
        memory="16GB",
        cpu=1
    shell: 
        """
        (echo "`date -R`: Building indexes..." &&
        bowtie2-build  \
        {input} \
        {params.base} &&
        echo "`date -R`: Success! Indexes are build." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
        """
