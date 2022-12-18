
rule prepareOperon:
    input:
        "resources/ref_genomes/{build}/operons.txt",
    output:
        "resources/ref_genomes/{build}/operons.bed",
    params:
        chromosome=chromosome(),
    log:
        "logs/rule/analysis/{build}/operonsBed.log",
    benchmark:
        "logs/rule/analysis/{build}/log/operonBed.benchmark.txt",
    resources:
        memory="2GB",
        cpu=1
    shell: 
        re.sub(' +', ' ', 
        """
        (echo "`date -R`: Converting operons to bed..." &&
        python3 workflow/scripts/operonbed.py -i {input} \
        -c {params.chromosome} \
        | grep "-" \
        > {output} &&
        echo "`date -R`: Success! Operon bed file is generated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
        """
        )