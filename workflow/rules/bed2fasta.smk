

rule bed2fasta:
    input:
        bed="results/{sample}/{sample}_lengthMode.bed",
        genome="resources/ref_genomes/{build}/genome.fa",
    output:
        temp("results/{sample}/{sample}_lengthMode.fa"),
    log:
        "logs/rule/analysis/{sample}/{sample}_bed2fasta.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_bed2fasta.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        (echo "`date -R`: Converting {input.bed} to fasta format..." &&
        bedtools getfasta \
        -fi {input.genome} \
        -bed {input.bed} \
        -fo {output} \
        -s &&
        echo "`date -R`: Success! {input.bed} is converted." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) > {log} 2>&1
        """

rule bed2fasta_input:
    input:
        bed=lambda w: input4fasta(w, config["input"]["files"], 
            config["srr"]["enabled"], config["srr"]["codes"]),
        genome="resources/ref_genomes/{build}/genome.fa",
    output:
        "results/input/{sample}/{sample}.fasta",
    log:
        "logs/rule/analysis/{sample}/{sample}_bed2fasta_input.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_bed2fasta_input.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        (echo "`date -R`: Converting {input.bed} to fasta format..." &&
        bedtools getfasta \
        -fi {input.genome} \
        -bed {input.bed} \
        -fo {output} \
        -s &&
        echo "`date -R`: Success! {input.bed} is converted." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) > {log} 2>&1
        """