
rule genome_index2ron:
    input:
        "resources/ref_genomes/{build}/genome.fa.fai",
    output:
        "resources/ref_genomes/{build}/genome.ron",
    log:
        "logs/rule/analysis/{build}/log/indexing.log",
    benchmark:
        "logs/rule/analysis/{build}/log/indexing.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    shell:
        "python3 workflow/scripts/idx2ron.py -i {input} -o {output} -l {log}"