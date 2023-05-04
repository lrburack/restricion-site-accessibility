#!/bin/sh
# Takes a bam file [1]

# Inteded for comparison with digestions
# Only read ends on the plus strand
OUTFILE_BEFORE="readEnds_plus_${1%.*}_"
OUTFILE_AFTER=".txt"

# Filter plus entries for those whose cigar string starts with a match
samtools view -h -F 16 $1 | awk '($6~/^[0-9]+M/) && ($10~/^'$2'/) {print $4 >> "'$OUTFILE_BEFORE'"$3"'$OUTFILE_AFTER'"}'