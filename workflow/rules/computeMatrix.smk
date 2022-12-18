rule compute_matrix_for_operon_start:
    input:
         # Please note that the -R and -S options are defined via input files
         bed="resources/ref_genomes/{build}/operons.bed",
         bigwig=lambda w: input4computeMatrix(),
    output:
        # Please note that --outFileName, --outFileNameMatrix and --outFileSortedRegions are exclusively defined via output files.
        # Usable output variables, their extensions and which option they implicitly call are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/deeptools/computematrix.html.
        matrix_gz="results/{sample}/{sample}_sorted_located_matrix.gz",   # required
        # optional output files
        matrix_tab="results/{sample}/{sample}_sorted_located_matrix.tab",
        matrix_bed="results/{sample}/{sample}_sorted_located_matrix.bed",
    log:
        "logs/rule/analysis/{sample}/{sample}_locateDamage_matrix.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_locateDamage_matrix.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1,
    params:
        # required argument, choose "scale-regions" or "reference-point"
        command="reference-point",
        # optional parameters
        extra="--referencePoint TSS -b 3000 -a 3000 --verbose --skipZeros -bs 100"
    wrapper:
        "v1.3.2/bio/deeptools/computematrix"

rule compute_matrix_for_operon_start_sim:
    input:
         # Please note that the -R and -S options are defined via input files
         bed="resources/ref_genomes/{build}/operons.bed",
         bigwig=lambda w: input4computeMatrix_sim(),
    output:
        # Please note that --outFileName, --outFileNameMatrix and --outFileSortedRegions are exclusively defined via output files.
        # Usable output variables, their extensions and which option they implicitly call are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/deeptools/computematrix.html.
        matrix_gz="results/{samples}/simulation/{samples}_sorted_located_matrix.gz",   # required
        # optional output files
        matrix_tab="results/{samples}/simulation/{samples}_sorted_located_matrix.tab",
        matrix_bed="results/{samples}/simulation/{samples}_sorted_located_matrix.bed",
    log:
        "logs/rule/analysis/{samples}/simulation/{samples}_locateDamage_matrix.log",
    benchmark:
        "logs/rule/analysis/{samples}/simulation/{samples}_locateDamage_matrix.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1,
    params:
        # required argument, choose "scale-regions" or "reference-point"
        command="reference-point",
        # optional parameters
        extra="--referencePoint TSS -b 3000 -a 3000 --verbose --skipZeros -bs 100"
    wrapper:
        "v1.3.2/bio/deeptools/computematrix"