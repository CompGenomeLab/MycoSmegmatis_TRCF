#!/bin/env python

#### configuration file ####
configfile: "config/myco_rna.yaml"

def getSampleFields(sample):
    l = []
    for key in sorted(sample.keys()):
        l.append(str(sample[key]))
    return '\t'.join(l)

def input4merge(sampleList):
    inputList = []
    for sample in sampleList:
        inputList.append(f"results/rna/pe/mapped/{sample}.tsv")
    return inputList

include: "workflow/rules/rna.smk"

output = ["results/rna/readCounts.tsv"]
# for accession in config['sample']:
    # output.append(f"results/rna/pe/mapped/{accession}.tsv")

rule all:
    input: output
        
 