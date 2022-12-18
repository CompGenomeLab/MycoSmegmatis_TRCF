rule simulation:
    input:
        bed= rules.sort_filter.output,
        genome="resources/ref_genomes/{build}/genome.fa",
        bowtie2="resources/ref_genomes/{build}/Bowtie2/genome.1.bt2", 
        inpfile="resources/ref_genomes/{build}/genome.ron",
    output:
        fa=temp("results/{samples}/simulation/{samples}.fa"),
        sim="results/{samples}/simulation/{samples}_sim.fa",
        tembed=temp("results/{samples}/simulation/{samples}.temp.bed"),
        simbed="results/{samples}/simulation/{samples}.bed",
    params:
        ref_genome="resources/ref_genomes/{build}/Bowtie2/genome",
        boquila_seed=config["boquila"]["seed"],
        boquila_sens=config["boquila"]["sens"],
    log:
        "logs/rule/analysis/{samples}/simulation/{samples}.log",
    benchmark:
        "logs/rule/analysis/{samples}/simulation/{samples}.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    conda:
        "../envs/simulation.yaml"
    shell:
        """
        (echo "`date -R`: Converting {input.bed} to fasta format..." &&
        bedtools getfasta \
        -fi {input.genome} \
        -bed {input.bed} \
        -fo {output.fa} \
        -s &&
        echo "`date -R`: Success! {input.bed} is converted." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Simulating reads according to reference genome..." &&
        boquila \
        --fasta {output.fa} \
        --bed {output.tembed} \
        --ref {input.genome} \
        --seed {params.boquila_seed} \
        --sens {params.boquila_sens} \
        --regions {input.inpfile} \
        > {output.sim} &&
        sort -k1,1 -k2,2n -k3,3n {output.tembed} > {output.simbed} &&
        echo "`date -R`: Success! Simulation is done." || 
        {{ echo "`date -R`: Process failed..."; rm {output.simbed}; exit 1; }}  ) >> {log} 2>&1
        """
