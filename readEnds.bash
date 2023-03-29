#!/bin/sh
# Returns a list of read end indexable entries from a bam file

# For an entry to be read end indexable it must either:
# Be reverse complemented AND end with a match
	# Returned ind = ind + bases after first match -1
# Not be reverse complemented AND start with a match
	# Returned ind = ind

PLUS=$(samtools view -F 16 $1)
MINUS=$(samtools view -f 16 $1)

# Filter plus entries for those whose cigar string starts with a match
PLUS_FILTERED=$(echo "$PLUS" | awk '$6~/^[0-9]+M/{print substr($10, 0, 4)}')
# Filter minus entries for those whose cigar string ends with a match
MINUS_FILTERED=$(echo "$MINUS" | awk '$6~/[M]$/{print substr($10, 92, 4)}')

echo "$PLUS_FILTERED
$MINUS_FILTERED"