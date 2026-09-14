#!/usr/bin/env python

import os, sys, re

from CdbTools import *

def contigs_to_scaffold_quals (ordered_contigs, contig_quals_file):
    """Builds the scaffold quality file according to the contig layout and
    individual contig quality values"""

    scaff_quals = []

    for contig in ordered_contigs:
        accession = contig['accession']
        match = re.search(r"^(\d+)_", accession)
        if match is None:
            raise Exception("Cannot parse contig number from accession: " + accession)

        contig_id = match.group(1)
        qual_entry = cdbyank(contig_id, contig_quals_file)
        qual_entry = qual_entry.rstrip()
        qual_lines = qual_entry.split("\n")
        qual_lines.pop(0) # rid header
        qual_lines = " ".join(qual_lines)
        qual_lines = qual_lines.lstrip()
        qual_lines = qual_lines.rstrip()
        qual_vals = re.split(r"\s+", qual_lines)
        scaff_quals += qual_vals

        if contig['gap'] > 0 :
            scaff_quals += ["0" for i in range(contig['gap']) ]
        elif contig['gap'] < 0:
            scaff_quals = scaff_quals[:contig['gap']]

    scaff_qual_text = " ".join(scaff_quals)
    
    return(scaff_qual_text)


def contigs_to_scaffold_sequence (ordered_contigs, contig_fasta_file):
    """Builds the scaffold sequence according to the contig layout and individual
    contig fasta sequences"""

    scaff_seq = ""
    
    for contig in ordered_contigs:
        accession = contig['accession']
        match = re.search(r"^(\d+)_", accession)
        if match is None:
            raise Exception("Cannot parse contig number from accession: " + accession)

        contig_id = match.group(1)
        (header, fasta_seq) = cdbyank_linear(contig_id, contig_fasta_file)
        scaff_seq += fasta_seq

        if contig['gap'] > 0 :
            scaff_seq += ("N" * contig['gap'])
        elif contig['gap'] < 0:
            scaff_seq = scaff_seq[:contig['gap']] # strip off negative gap length from contig
        
    return(scaff_seq)
    


def build_ordered_contigs(scaff_file):
    fh = open(scaff_file, 'r')

    contigs = []

    scaffold_accession = ""
    
    for line in fh.readlines():
        if line[0] == '>':
            if len(contigs) != 0:
                raise Exception("Error! more than one scaffold reported.")
            
            header_tokens = line.split()
            scaffold_accession = header_tokens[0]
            scaffold_accession = scaffold_accession[1:]
            continue
        
        line = line.rstrip()
        (accession, orient, length, gap) = line.split()

        contig = { 'accession':accession,
                   'orient':orient,
                   'length':length,
                   'gap': int(gap) }

        contigs.append(contig)


    return(scaffold_accession, contigs)




from optparse import OptionParser


if __name__ == '__main__':
    
    parser = OptionParser()
    parser.add_option("-L", "--scaffold_layout", dest="scaffold_layout_file",
                      help="scaffold layout filename")
    parser.add_option("-C", "--contigs_quals", dest="contigs_quals_file",
                      help="contig quals file")
    parser.add_option("-S", "--contigs_fasta", dest="contigs_fasta_file",
                      help="contigs fasta file")
    
    (options, args) = parser.parse_args()

    if not (options.scaffold_layout_file and options.contigs_quals_file and options.contigs_fasta_file):
        raise Exception("\n\n\n****** use -h for usage info ******* \n\n\n\n")


    (scaffold_accession, ordered_contigs) = build_ordered_contigs(options.scaffold_layout_file)
    

    scaff_seq = contigs_to_scaffold_sequence(ordered_contigs, options.contigs_fasta_file)
    print ">" + scaffold_accession + "\n" + scaff_seq + "\n"

    scaff_quals = contigs_to_scaffold_quals(ordered_contigs, options.contigs_quals_file)
    print ">" + scaffold_accession + "\n" + scaff_quals + "\n"
    
    sys.exit(0)



    
