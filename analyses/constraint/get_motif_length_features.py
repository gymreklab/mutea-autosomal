#!/usr/bin/env python

import pyfaidx
import sys

try:
    lobref = sys.argv[1]
    REFDB = sys.argv[2]
except:
    sys.stderr.write("Usage: ./get_motif_length_features.py <lobref> <reffa>\n")
    sys.exit(1)

def ReverseComplement(seq):
    """ Get the reverse complement of a nucleotide string """
    newseq = ""
    size = len(seq)
    for i in range(len(seq)):
        char = seq[len(seq)-i-1]
        if char == "A": newseq += "T"
        if char == "G": newseq += "C"
        if char == "C": newseq += "G"
        if char == "T": newseq += "A"
    return newseq

def GetUninterruptedLength_(seq, motif):
    loc = 0
    starts = []
    ends = []
    ind = seq.find(motif, loc) # get first one
    if ind != -1:
        starts.append(ind)
        loc = ind + len(motif)
#    print "starts", starts
    while True:
#        print "looking"
        ind = seq.find(motif, loc)
        if ind == -1: # didn't find anymore, end run
            if len(starts) > 0:
                ends.append(loc-1)
            break
#        print "ind, loc", ind, loc
        if ind > loc: # imperfect, end run and start new one
            ends.append(min([loc-1, len(seq)-1]))
            starts.append(ind)
        loc = ind + len(motif)
#        print "starts, ends",  starts, ends
    intervals = zip(starts, ends)
    intervals = map(list, intervals)
#    print "intervals", intervals
    # For each interval, try to extend
    maxlen = 0
    bestint = [-1, -1]
    for i in range(len(intervals)):
        subseq = seq[intervals[i][0]:intervals[i][1]+1]
#        print "i, subseq", i, subseq, intervals[i]
        prefix = seq[max([intervals[i][0]-len(motif), 0]):intervals[i][0]]
        suffix = seq[intervals[i][1]+1:min([intervals[i][1]+1+len(motif), len(seq)])]
        length = intervals[i][1]-intervals[i][0] # full match
        for j in range(1, len(motif)+1): # match prefix
            if len(prefix) < j: break
#            print "prefix", j, prefix[-1*j], motif[-1*j]
            if prefix[-1*j] == motif[-1*j]:
                intervals[i][0] -= 1
            else:
                break
        for j in range(len(motif)): # match suffix
            if len(suffix) <= j: break
#            print "suffix", j, motif[j], suffix[j]
            if suffix[j] == motif[j]:
                intervals[i][1] += 1
            else:
                break
#        print "extended", seq[intervals[i][0]:intervals[i][1]+1]
        if intervals[i][1]-intervals[i][0]+1 > maxlen:
            maxlen = intervals[i][1]-intervals[i][0]+1
            bestint = intervals[i]
#    print bestint, maxlen, seq[bestint[0]:bestint[1]+1]
    return maxlen

def test():
    seq="TTTAGTTTATTTTTTATTTTTATTTATTTATTTATTTATTT"
    motif="ATTT"
    assert GetUninterruptedLength_(seq, motif) == 23
    seq="ATATATATATATATAT"
    motif="AT"
    assert GetUninterruptedLength_(seq, motif) == 16
    seq="TGTGTGTGTGTG"
    motif="GT"
    assert GetUninterruptedLength_(seq, motif) == 12
    seq="TCTTTCTTTCTCCTTTTCTTCTCTTTTCTTTCTTTCTTTTCTTTCTTTCTTCTTCGTTCTTTC"
    motif=ReverseComplement("AAAG")
    assert GetUninterruptedLength_(seq, motif) == 15
#test()
#sys.exit(1)

def GetUninterruptedLength(chrom, start, end, motif, genome):
    seq = genome[chrom][start-1:end].upper()
    rmotif = ReverseComplement(motif)
    sub1 = GetUninterruptedLength_(seq, motif)
    sub2 = GetUninterruptedLength_(seq, rmotif)
    return max([sub1, sub2])

# Load reference genome
genome = pyfaidx.Fasta(REFDB, as_raw=True)

# Get features for each STR
with open(lobref, "r") as f:
    for line in f:
        items = line.strip().split()
        chrom = items[0]
        start = int(items[1])
        end = int(items[2])
        motif = items[14]
        length = end - start + 1
        unint_length = GetUninterruptedLength(chrom, start, end, motif, genome)
        print "\t".join(map(str, [chrom, start, end, motif, length, unint_length]))
