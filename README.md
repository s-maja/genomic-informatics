# genomic-informatics
Repository for the homework and the project done within the course 'Genomics informatics'

This repo contains an implementation of a Burrows-Wheeler algorithm and Needleman-Wunsch algorithm and then combined using Seed and Extand method to match reads on refernce genome.

* Implementation of Burrows-Wheeler Transformation and FM index
* Algorithm for global alignment aka Needleman-Wunsch algorithm
* Usage of “Seed and Extend” method. Input files format: .fasta, .fastq
* Comparing results with BWA MEM tool and presenting results graphically and tabularly


# Burrows-Wheel Transformation and FM index

This algorithm rearranges a string into runs of similar character and it is very usful for compression. This transformation is reversible. BWT in combination with FM index is off-line exact string matching algorithm.


# Needleman-Wunsch algorithm

This is global alignment technique and applications of dynamic programming for comparing biological sequences. This algorithm uses weighted edit-distance (scoring matrix).

# Seed and Extand

This method takes two steps:

1)find seed (first n characters from read) and map seed to reference
2)take rest characters from read and align it to reference (also check for read reverse complement)
