rule report_stats:
    input:
        originalFastq = "resources/samples/{sample}.fastq",
        nodupFastq = "results/{sample}/{sample}_removedDup.fastq",
        cutFastq = "results/{sample}/{sample}_cut.fastq",
        bam = "results/{sample}/{sample}.bam",
        filteredBed = "results/{sample}/{sample}.bed",
        bedLength = "results/{sample}/{sample}_bedLength.txt",
    output: report("results/{sample}/{sample}_report.txt", category="report"), 
    log:
        "logs/rule/analysis/{sample}/{sample}_reportStats.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_reportStats.benchmark.txt",
    resources:
        memory="4GB",
        cpu=1
    conda:
        "../envs/align.yaml"
    shell:
        """
        (echo "`date -R`: Preparing a report for read statistics..." &&
        A1=$(grep -c "^+" {input.originalFastq}) &&
        A2=$(grep -c "^+" {input.nodupFastq}) &&
        A3=$(grep -c "^+" {input.cutFastq}) &&
        A4=$(samtools view -c -F 260 {input.bam}) &&
        A5=$(grep -c "^" {input.filteredBed}) &&
        printf "Original read number from {input.originalFastq}: $A1\\n" > {output} &&
        printf "Deduplicated read number from {input.nodupFastq}: $A2\\n" >> {output} &&
        printf "Read number after adaptor trimming from {input.cutFastq}: $A3\\n" >> {output} &&
        printf "Aligned read number from {input.bam}: $A4\\n" >> {output} &&
        printf "Filtered read number from {input.filteredBed}: $A5\\n" >> {output} &&
        printf "Length distribution of the bed file:\\n" >> {output} &&
        cat {input.bedLength} >> {output} &&
        printf "\\n" >> {output} &&
        echo "`date -R`: Success! Report for read statistics is generated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
        """