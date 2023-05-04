#!/bin/sh
# Takes a bam file [1] and a restriction motif [2]

# Inteded for comparison with digestions
# Only read ends on the plus strand
OUTFILE_BEFORE="readEnds_plus_${1%.*}_"
OUTFILE_AFTER=".txt"

MOTIF_LENGTH=$(expr length $2)

# Filter plus entries for: 
#   Cigar string starts with match
#   Sequence starts with restriction motif
samtools view -F 16 $1 | awk '($6~/^[0-9]+M/) && ($10~/^'$2'/) {print $4 >> "'$OUTFILE_BEFORE'"$3"'$OUTFILE_AFTER'"}'
# Filter minus entries for:
#   Cigar string is a full match
#   Sequence ends with restriction motif
# Shift index by seqlength - motif_length
samtools view -h -f 16 $1 | awk '($6~/^[0-9]+M$/) && ($10~/'$2'$/) {print $4 + length($10) - '$MOTIF_LENGTH' >> "'$OUTFILE_BEFORE'"$3"'$OUTFILE_AFTER'"}'