
rule igv_report:
    input:
        fasta="resources/ref_genomes/{build}/genome_chr{chr}.fa",
        vcf="resources/ref_genomes/{build}/variation_bgzip.vcf.gz",
        tracks=["results/{sample}/{sample}_sorted_plus.bed", "results/{sample}/{sample}_sorted_minus.bed"],
    output:
        report("results/{sample}/{sample}_igv_report_chr{chr}.html", category="IGV"),
    params:
        extra="",  
    log:
        "logs/rule/analysis/{sample}/{sample}_igv_report_chr{chr}.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_igv_report_chr{chr}.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    wrapper:
        "0.69.0/bio/igv-reports"