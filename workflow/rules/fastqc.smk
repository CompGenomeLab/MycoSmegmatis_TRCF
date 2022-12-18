configfile: "config/config_XR_initial.yaml" 

rule fastqc_se:
    input:
        "resources/samples/{samples}.fastq.gz",
    output:
        html=report("results/{sample}/{sample}.html", category="QC"),
        zip="results/{sample}/{sample}_fastqc.zip",
    params: ""
    log:
        "logs/rule/analysis/{sample}/{sample}_fastqc.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_fastqc.benchmark.txt",
    resources:
        memory=config["resources"]["memory"],
        cpu=config["resources"]["cpu"],
    # threads: 1
    wrapper:
        "0.69.0/bio/fastqc"

rule fastqc_pe:
    input:
        "resources/samples/{samples}_{ext}.fastq.gz", 
    output:
        html=report("results/{sample}/{sample}_{ext}.html", category="QC"), 
        zip="results/{sample}/{sample}_{ext}_fastqc.zip", 
    params: ""
    log:
        "logs/rule/analysis/{sample}/{sample}_fastqc_{ext}.log", 
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_fastqc_{ext}.benchmark.txt",
    resources:
        memory=config["resources"]["memory"],
        cpu=config["resources"]["cpu"],
    # threads: 1
    wrapper:
        "0.69.0/bio/fastqc"