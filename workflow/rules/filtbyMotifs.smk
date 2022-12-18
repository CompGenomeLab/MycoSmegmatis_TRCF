
rule filtbyMotifs:
    input:
        "results/{sample}/{sample}_sorted_{strand}_10.fa",
    output:
        "results/{sample}/{sample}_sorted_ds_dipyrimidines_{strand}.bed",
    params:
        lambda w: getMotif(w.samples, config["meta"][w.samples]["product"]),
    log:
        "logs/rule/analysis/{sample}/{sample}_filtbyMotifs_{strand}.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_filtbyMotifs_{strand}.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    shell:  
        """
        (echo "`date -R`: Filtering by the given motif..." &&
        workflow/scripts/fa2bedByChoosingReadMotifs.py \
        -i {input} \
        -o {output} \
        -r {params} &&
        echo "`date -R`: Success! Filtering is done." ||
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) > {log} 2>&1
        """