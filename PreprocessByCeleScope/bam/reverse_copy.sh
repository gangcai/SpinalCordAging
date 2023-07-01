#!/bin/bash
for sample in SC-O-21M SC-O-23M SC-Y-4M SC-Y-6M
do
	cp ${sample}.bam.bai ../${sample}/${sample}/04.featureCounts/${sample}_Aligned.sortedByCoord.out.bam.featureCounts.bam.bai
done
