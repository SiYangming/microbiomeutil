#!/usr/bin/env python


import os, sys, re, tempfile
from FastaReader import *
from CmdRunner import *
from CdbTools import *



def homologyTrimLucy(seqFile, qualFile, referenceSeqFile, log=sys.stderr, minPercentIdentity=0.0):

    log.write("homologyTrimLucy(seqFile = " + seqFile \
              + " qualFile=" + qualFile \
              + " referenceSeqFile=" + referenceSeqFile + "\n")
                 
    
    trimSeqFH = open(seqFile + ".htrim.seq", 'w')
    trimQualFH = open(seqFile + ".htrim.qual", 'w')
    
    seqsProcessed = []
    fr = FastaReader(seqFile)
    
    while(True):
        seqObj = fr.nextEntry()
        if not seqObj:
            break

        seqAccession = seqObj.accession
        
        (fd, tmpfilename) = tempfile.mkstemp()
                
        os.write(fd, ">" + seqObj.header + "\n" + seqObj.sequence)
        os.close(fd)


        if not os.path.isfile(referenceSeqFile + ".nin"):
            cmd = "formatdb -i " + referenceSeqFile + " -p F"
            CmdRunner.process_cmd(cmd)
        
        results = None
        
        try:
            #cmd = "megablast -i " + tmpfilename + " -d " + referenceSeqFile \
            cmd = "blastall -p blastn -i " + tmpfilename + " -d " + referenceSeqFile \
                  + " -v 1 -b 1 -m 8 -e 1e-20 "
            log.write("RefHomologyCheck: " + cmd + "\n")
            results = CmdRunner.process_cmd(cmd)
        
        except (CmdRunnerException):
            # nonzero exit of megablast should be treated as having no hits (unfortunately)
            log.write("No matches detected for " + seqAccession + " according to nonzero megablast exit value.\n")
            
        
        os.remove(tmpfilename)
        
                

        if results is None or re.search(r"\w", results) is None:
            log.write("Warning, no match found for " + seqAccession + "\n")
            continue

       
        log.write(results)
        
        results = results.rstrip()
        matches = results.split("\n")
        match_lends = []
        match_rends = []

        for match in matches:
            print "Match: " + match
            match = match.rstrip()
            x = match.split("\t")
            perID = float(x[2])
            if perID < minPercentIdentity:
                continue
            
            query_lend = int(x[6])
            query_rend = int(x[7])
            match_lends.append(query_lend)
            match_rends.append(query_rend)

            print match_lends
            print match_rends

        if match_lends:
            log.write("Match lends: " + ",".join([str(x) for x in match_lends]) + "\n")
            log.write("Match rends: " + ",".join([str(x) for x in match_rends]) + "\n")
        
            query_lend = min(match_lends)
            query_rend = max(match_rends)
            
            log.write("Match coordinates:" + str(query_lend) + "-" + str(query_rend) + "\n")
            
            # check the Lucy trim coordinates:
            header_entries = re.split(r"\s+", seqObj.header)
            lucy_rend = int(header_entries.pop())
            lucy_lend = int(header_entries.pop())
            log.write("Lucy trim coordinates: " + str(lucy_lend) + "-" + str(lucy_rend) + "\n")
                    
            trim_lend = max(query_lend, lucy_lend)
            trim_rend = min(query_rend, lucy_rend)
            log.write("Final trim: " + str(trim_lend) + "-" + str(trim_rend) + "\n")
                                        
            htrim_seq = seqObj.getRawSequence()[trim_lend-1:trim_rend]
            trimSeqFH.write(">" + seqObj.accession + " trimmed to " + str(trim_lend) + " " + str(trim_rend) \
                            + "\n" + htrim_seq + "\n")
                                        
            write_trimmed_qual_file(seqObj.accession, qualFile, trim_lend, trim_rend, trimQualFH)
            seqsProcessed.append(seqAccession)

        else:
            log.write("Warning, no maqtch found for " + seqAccession + " after applying minPerID threshold\n")

    
    return(seqsProcessed)



def write_trimmed_qual_file(accession, qual_file, lend, rend, FH):
    qualEntry = cdbyank(accession, qual_file)
    lines = qualEntry.split("\n")

    lines.pop(0) # remove header
    val_string = " ".join(lines)
    vals = re.split(r"\s+", val_string)

    vals = vals[lend-1:rend]

    FH.write(">" + accession + " trimmed to " + str(lend) + " " + str(rend) \
             + "\n" + " ".join(vals) + "\n")

    return





if __name__ == '__main__':
    
    usage = "\n\n\nusage: " + sys.argv[0] + " reads.seq reads.qual reference.seq\n\n\n"

    if (len(sys.argv) != 4):
        raise Exception(usage)

    (seqFile, qualFile, referenceSeqFile) = sys.argv[1:]

    homologyTrimLucy(seqFile, qualFile, referenceSeqFile)

    sys.exit(0)


