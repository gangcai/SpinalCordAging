#!/bin/bash
for sample in SC-Y-6M SC-O-21M
do
	echo $sample
	conda run -n celescope1.11.0 bash ./generate_shell_scripts.sh ${sample}
done
