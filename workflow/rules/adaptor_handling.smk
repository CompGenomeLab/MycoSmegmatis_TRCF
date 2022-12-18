
rule cutadapt_se:
    input:
        rules.removeDuplicatesAtFastq.output,
    output:
        fastq="results/{sample}/{sample}_cut.fastq",
        qc=report("results/{sample}/{sample}_cutqc.txt", category="QC"),   
    params:
        adapters=config["adaptor_se"],
        extra='--discard-untrimmed'  
    log:
        "logs/rule/analysis/{sample}/{sample}_cut.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_cut.benchmark.txt",
    resources:
        memory=config["resources"]["memory"],
        cpu=config["resources"]["cpu"],
    wrapper:
        "0.75.0/bio/cutadapt/se"
