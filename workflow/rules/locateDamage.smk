
rule locateDamage:
    input:
        plus=rules.sep_strands.output.plus,
        minus=rules.sep_strands.output.minus,
    output:
        plus="results/{sample}/{sample}_fslp.bed",
        minus="results/{sample}/{sample}_fslm.bed",
    log: 
        "results/{sample}/{sample}_fsl.log",
    benchmark: 
        "results/{sample}/{sample}_fsl.benchmark.txt",
    resources:
        memory="2GB",
        cpu=1,
    shell:  
        re.sub(' +', ' ', 
        """
        (echo "`date -R`: Locating exact damage sites for bacterial repair map..." &&
        cat {input.plus} | python3 workflow/scripts/locateDamage.py > {output.plus} &&
        cat {input.minus} | python3 workflow/scripts/locateDamage.py > {output.minus} &&
        echo "`date -R`: Success! Damage location is done." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) > {log} 2>&1
        """
        )


use rule locateDamage as locateDamage_sim with:
    input:
        plus=rules.sep_strands_sim.output.plus,
        minus=rules.sep_strands_sim.output.minus,
    output:
        plus="results/{samples}/simulation/{samples}_fslp.bed",
        minus="results/{samples}/simulation/{samples}_fslm.bed",
    log: 
        "results/{samples}/simulation/{samples}_fsl.log",
    benchmark: 
        "results/{samples}/simulation/{samples}_fsl.benchmark.txt",