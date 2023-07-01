#!/bin/bash
sample="${1}"
mapfile="${sample}.mapfile"
genome_dir="/home/db/public/SoftwareIndex/celescope1.11.0/mouse"
thread=10
expected_cell_num=10000
multi_rna\
	--mapfile ${mapfile}\
	--genomeDir ${genome_dir}\
	--thread ${thread}\
	--outdir ${sample}\
	--gtf_type "gene"\
	--expected_cell_num ${expected_cell_num}\
	--steps_run "all"\
	--allowNoPolyT \
	--mod shell
