#!/bin/bash
for sample in SC-O-21M SC-O-23M SC-Y-4M SC-Y-6M
do
	samtools index ${sample}.bam ${sample}.bam.bai
done
