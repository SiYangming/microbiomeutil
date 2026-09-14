#!/usr/bin/env python

import os, sys, re

class FastaReader:

    def __init__(self, fastaFile):
        self.filename = fastaFile
        self.fh = open(fastaFile, 'r')
        self.curr_line = "\n";

    def nextEntry(self):
        while (self.curr_line and self.curr_line.find(">") != 0):
            self.curr_line = self.fh.readline()
        if not self.curr_line:
            # end of file reached
            return(None)
        header = self.curr_line
        header = header.rstrip()
        header = header[1:] # remove carat
        accession = header.split()[0]

        sequenceLines = []
        self.curr_line = self.fh.readline()
        while (self.curr_line and self.curr_line.find(">") != 0):
            sequenceLines.append(self.curr_line)
            self.curr_line = self.fh.readline()

        sequence = "".join(sequenceLines)

        fastaSequenceObject = FastaSequence(accession, header, sequence)
        return(fastaSequenceObject)


    def close(self):
        self.fh.close()
        

class FastaSequence:

    def __init__(self, accession, header, sequence):
        self.accession = accession
        self.header = header
        self.sequence = sequence

    def getRawSequence(self):
        seq = self.sequence
        seq = re.sub(r"\s", "", seq) # remove all whitespace
        return(seq)


if __name__ == '__main__':
    usage = "\n\nusage: " + sys.argv[0] + " fastaFile\n\n";
    if (len(sys.argv) < 2):
        raise Exception(usage)

    fasta_filename = sys.argv[1]
    fasta_reader = FastaReader(fasta_filename)
    fastaSeqObj = fasta_reader.nextEntry()
    while(fastaSeqObj):
        print "Accession: " + fastaSeqObj.accession
        print "Header: " + fastaSeqObj.header
        print "Sequence: " + fastaSeqObj.getRawSequence() + "\n\n"
        fastaSeqObj = fasta_reader.nextEntry()
    
    sys.exit(0)

