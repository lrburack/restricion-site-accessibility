#!/bin/sh
# Takes a bam file with an index file in the same directory [1], a reference [2], a restriction site [3]

# For now only sites at a read end on the plus strand are counted

# Filter only for provided ref on plus strand
PLUS=$(samtools view -h -F 16 $1 $2)

# Filter plus entries for those whose cigar string starts with a match
# with a digestion overhang at the read end
PLUS_FILTERED=$(echo "$PLUS" | awk '$6~/^[0-9]+M/ && $10~/^'$3'/ {print $4}')

echo "$PLUS_FILTERED"