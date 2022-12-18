rule mark_duplicates_se:
    input: rules.bam2bed_se.output.bam
    # optional to specify a list of BAMs; this has the same effect
    # of marking duplicates on separate read groups for a sample
    # and then merging
    output:
        bam="results/{sample}/{sample}_cutadapt_se_dedup.bam",
        metrics="results/{sample}/{sample}_cutadapt_se_dedup.metrics.txt",
    log:
        "logs/rule/analysis/{sample}/{sample}_cutadapt_se_dedup.log",
    params:
        extra="--REMOVE_DUPLICATES true",
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=1024,
        memory="16GB",
        cpu=1
    wrapper:
        "v1.3.2/bio/picard/markduplicates"

# rule mark_duplicates_pe:
#     input: "results/{sample}/{sample}_cutadapt_pe.bam"
#     # optional to specify a list of BAMs; this has the same effect
#     # of marking duplicates on separate read groups for a sample
#     # and then merging
#     output:
#         bam="results/{sample}/{sample}_cutadapt_pe_dedup.bam",
#         metrics="results/{sample}/{sample}_cutadapt_pe_dedup.metrics.txt",
#     log:
#         "logs/rule/analysis/{sample}/{sample}_cutadapt_pe_dedup.log",
#     params:
#         extra="--REMOVE_DUPLICATES true",
#     # optional specification of memory usage of the JVM that snakemake will respect with global
#     # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
#     # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
#     # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
#     resources:
#         mem_mb=1024,
#         memory="16GB",
#         cpu=1
#     wrapper:
#         "v1.3.2/bio/picard/markduplicates"