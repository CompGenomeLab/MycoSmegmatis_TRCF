
rule plot_length:
    input:
        rules.fastq2length.output.original,
    output:
        report("results/{sample}/{sample}_length_distribution.pdf", 
                category="Length Distribution"),
    params:
        "{sample}",
    log:
        "logs/rule/figs/{sample}/{sample}_plot_length.log",
    benchmark:
        "logs/rule/figs/{sample}/{sample}_plot_length.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    conda:
        "../envs/plot_nuc.yaml"
    shell:  
        """
        Rscript workflow/scripts/lengthDistPlot.r \
        -i {input} \
        -s {params} \
        -o {output} \
        -l {log}
        """