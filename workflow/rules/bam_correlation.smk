rule bam_correlation:
    input:
        lambda w: input4PCA(config["sample"], config["srr"]["enabled"], 
            config["srr"]["codes"], config["build"], w.duplicate),
    output:
        out="results/readCounts_{duplicate}.npz",
        raw_out="results/readCounts_{duplicate}.tab",
    log:
        "logs/rule/analysis/bam_correlation_{duplicate}.log",
    benchmark:
        "logs/rule/analysis/bam_correlation_{duplicate}.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    conda:
        "../envs/bam_correlation.yaml"
    shell:  
        """
        (echo "`date -R`: MultiBam summary..." &&
        multiBamSummary bins \
        --bamfiles {input} \
        --minMappingQuality 20 \
        -out {output.out} --outRawCounts {output.raw_out} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.out}; exit 1; }} ) > {log} 2>&1
        """
