use rule bam2bed_se as bam2bed_dedup with:
    input:
        bam= rules.mark_duplicates_se.output.bam,
        fasta="resources/ref_genomes/" + config["build"] + "/singletons.fa"
    output:
        bed="results/{sample}/{sample}_dedup.bed",
        bam="results/{sample}/{sample}_dedup_sortedbyCoordinates.bam",
        sam="results/{sample}/{sample}_dedup_sortedbyCoordinates.sam",
        idx="results/{sample}/{sample}_dedup_sortedbyCoordinates.bam.bai",
    log:
        "logs/rule/analysis/{sample}/{sample}_dedup_bam2bed.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_dedup_bam2bed.benchmark.txt",