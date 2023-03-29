#!/bin/sh
# Takes a bam file [1], a reference [2], a restriction site [3]

# For now only sites at a read end on the plus strand are counted

# Filter plus entries for those whose cigar string starts with a match
# with a digestion overhang at the read end
PLUS_FILTERED=$(samtools view -h -F 16 $1 | awk '$3~/'$2'/ && $6~/^[0-9]+M/ && $10~/^'$3'/ {print $4}')

echo "$PLUS_FILTERED"