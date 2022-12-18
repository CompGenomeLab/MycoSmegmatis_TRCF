project_ = config['project']

# rule length_dist:
#     input:
#         rules.cutadapt_se.output,
#     output:
#         "results/{sample}/{sample}_length.txt",
#     log:
#         "logs/rule/analysis/{sample}/{sample}_length.log",
#     benchmark:
#         "logs/rule/analysis/{sample}/{sample}_length.benchmark.txt",
#     resources:
#         memory="16GB",
#         cpu=1
#     shell:  
#         """
#         (echo "`date -R`: Calculating the read length distribution..." &&
#         awk '{{print $3-$2}}' {input} |&
#         sort -k1,1n |& 
#         uniq -c |& 
#         sed 's/\s\s*/ /g' |&
#         awk '{{print $2"\\t"$1}}' > {output} &&
#         echo "`date -R`: Success! Length distribution is calculated." || 
#         {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) > {log} 2>&1
#         """



rule fastq2length:
    input:
        rules.cutadapt_se.output.fastq,
    output:
        original="results/{sample}/{sample}_length.txt",
        nozero="results/{sample}/{sample}_lengthExc0.txt"
    params:
        fields = lambda w: getSampleFields(config["meta"][w.sample])
    log:
        "logs/rule/analysis/{sample}/{sample}_length.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_length.benchmark.txt",
    resources:
        memory="4GB",
        cpu=1
    shell:  
        re.sub(' +', ' ',
        """
        (echo "`date -R`: Calculating the read length distribution..." &&
        python3 workflow/scripts/fastq2lengthDist.py \
        -i {input} \
        -nozero {output.nozero} \
        | awk '{{print $0, "\\t{params.fields}"}}' \
        > {output.original} &&
        echo "`date -R`: Success! Length distribution is calculated." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) > {log} 2>&1
        """
        )



rule merge_length:
    input:
        lambda w: input4mergeLength(config),
    output:
        out=f"results/{project_}/length.tsv",
    log: f"logs/rule/analysis/{project_}/length.log",
    resources:
        memory="16GB",
        cpu=1
    shell:  
        """
        (echo "`date -R`: Merge read length distribution files..." &&
        cat {input} > {output.out} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.out}; exit 1; }} ) > {log} 2>&1
        """

