#!/bin/bash
for sample in SC-Y-4M SC-O-23M
do
	echo $sample
	conda run -n celescope1.11.0 bash ./generate_shell_scripts.sh ${sample}
done
