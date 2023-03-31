#!/bin/sh
# Takes a bam file [1], a restriction site [2]

# For now only sites at a read end on the plus strand are counted
OUTFILE_BEFORE="digestions_${1%.*}_"
OUTFILE_AFTER=".txt"

# Filter plus entries for those whose cigar string starts with a match
# with a digestion overhang at the read end
samtools view -h -F 16 $1 | awk '($6~/^[0-9]+M/) && ($10~/^'$3'/) {print $4 >> "'$OUTFILE_BEFORE'"$3"'$OUTFILE_AFTER'"}'