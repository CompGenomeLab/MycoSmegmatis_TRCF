#!/bin/env python

#### configuration file ####
configfile: "config/config_XR_initial.yaml" 
configfile: "config/myco.yaml"
include: "workflow/rules/common.smk"


rule all:
    input:
        lambda w: allInput(config)
 