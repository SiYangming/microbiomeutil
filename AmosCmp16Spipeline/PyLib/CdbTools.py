#!/usr/bin/env python

import os, sys
import CmdRunner


class Cdbyank_exception(Exception): pass

def cdbyank (acc, database):

    index_file = database + ".cidx"
    if not os.path.isfile(index_file):
        # build cdbfasta index file
        CmdRunner.process_cmd("cdbfasta " + database)
        
    cmd = "cdbyank -a '" + acc + "' " + index_file
    output = CmdRunner.process_cmd(cmd)
    return(output)

def cdbyank_linear (acc, database):
    output = cdbyank(acc, database)
    output = output.split("\n")

    header = output.pop(0)
    header = header[1:] # remove carat

    sequence = "".join(output)

    return(header, sequence)



if __name__ == '__main__':
    usage = "\n\n\nusage: " + sys.argv[0] + " accession database\n\n\n"
    if (len(sys.argv) < 3):
        raise Exception(usage)

    (acc, db) = sys.argv[1:3]
    
    # full fasta record retrieval
    output = cdbyank(acc, db)
    print output

    # just the sequence retrieval
    (header, sequence) = cdbyank_linear(acc, db)
    print "Header: " + header
    print "Sequence: " + sequence
    
    sys.exit(0)

