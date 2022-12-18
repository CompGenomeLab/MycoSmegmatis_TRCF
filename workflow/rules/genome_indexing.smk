
rule genome_indexing:
    input:
        "resources/ref_genomes/{build}/genome.fa",
    output:
        "resources/ref_genomes/{build}/genome.fa.fai",
    benchmark:
        "logs/rule/analysis/{build}/log/indexing.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    wrapper:
        "0.69.0/bio/samtools/faidx"