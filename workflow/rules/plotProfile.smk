rule plot_profile:
    input:
        # matrix file from deepTools computeMatrix tool
        "results/{sample}/{sample}_sorted_located_matrix.gz",
    output:
        # Please note that --outFileSortedRegions and --outFileNameData are exclusively defined via output files.
        # Usable output variables, their extensions and which option they implicitly call are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/deeptools/plotprofile.html.
        # Through the output variables image file and more output options for plot profile can be selected.
        plot_img=report("results/plots/{samples}_TSSprofile.png",
                    category="TSS Profile"),  # required
        # optional output files
        regions="results/{sample}/{sample}_sorted_located_matrix_regions.bed",  # required
        data="results/{sample}/{sample}_sorted_located_matrix_regions.tab",  # required
    log:
        "logs/rule/analysis/{sample}/{sample}_sorted_located_matrix_plot.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_sorted_located_matrix_plot.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1,
    params:
        # optional parameters
        "--dpi 150 "
    wrapper:
        "v1.3.2/bio/deeptools/plotprofile"

rule plot_heatmap:
    input:
        "results/{sample}/{sample}_sorted_located_matrix.gz",
    output:
        heatmap_img=report("results/plots/plot_heatmap/{samples}_heatmap.png", 
        category="heatmap"),  # required
    log:
        "logs/rule/analysis/{sample}/{sample}_sorted_located_matrix_heatmap.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_sorted_located_matrix_heatmap.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1,
    params:
        # "--missingDataColor 0 ",
        "--colorMap purples ",
    wrapper:
        "v1.4.0/bio/deeptools/plotheatmap"

rule plot_profile_sim:
    input:
        # matrix file from deepTools computeMatrix tool
        "results/{samples}/simulation/{samples}_sorted_located_matrix.gz",
    output:
        # Please note that --outFileSortedRegions and --outFileNameData are exclusively defined via output files.
        # Usable output variables, their extensions and which option they implicitly call are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/deeptools/plotprofile.html.
        # Through the output variables image file and more output options for plot profile can be selected.
        plot_img=report("results/plots/simulation/{samples}_TSSprofile.png",
            category="TSS Profile"),  # required
        # optional output files
        regions="results/{samples}/simulation/{samples}_sorted_located_matrix_regions.bed",  # required
        data="results/{samples}/simulation/{samples}_sorted_located_matrix_regions.tab",  # required
    log:
        "logs/rule/analysis/{samples}/simulation/{samples}_sorted_located_matrix_plot.log",
    benchmark:
        "logs/rule/analysis/{samples}/simulation/{samples}_sorted_located_matrix_plot.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1,
    params:
        # optional parameters
        "--dpi 150 "
    wrapper:
        "v1.3.2/bio/deeptools/plotprofile"
