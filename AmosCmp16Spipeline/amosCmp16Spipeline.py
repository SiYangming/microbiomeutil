#!/usr/bin/env python

import os, sys, re, shutil, traceback

sys.path.append(sys.path[0] + "/PyLib")

from CdbTools import *
from CmdRunner import *
from FastaReader import *

import homology_trim_lucy_output
import contig_layout_to_scaffold


DEBUG = False #True
minPercentIdentity = 90.0  # require at least this percent identity between sequence read and reference genome


#############################################################
## Must set path to the following programs in conf.txt file:

conf = {}
for line in open(sys.path[0] + "/conf.txt"):
    if not re.search(r"\w", line):
        continue
    
    match = re.search(r"^(\S+)\s+(.*)\s*$", line)
    if match:
        prog_dir_token = match.group(1)
        prog_dir = match.group(2)
        conf[prog_dir_token] = prog_dir
    else:
        print "Cannot parse conf.txt entry: " + line

#print conf


# update PATH environmental variable so tools will be found
required_dirnames = ['LUCY_DIR', 'MEGABLAST_DIR', 'AMOS_DIR']
for dirname in required_dirnames:
    if not conf.has_key(dirname):
        raise Exception("Please specify path corresponding to (" + dirname + ") in file conf.txt")
    os.environ['PATH'] += ":" + conf[dirname]


def main():

    usage = "\n\n\nusage: " + sys.argv[0] + " reads.fasta reads.qual reads.pairs reference_db.fasta [minPercentIdentity=90]\n\n"
    
    if (len(sys.argv) < 5):
        raise Exception(usage)

    (reads_file, qual_file, pairs_file, reference_db) = sys.argv[1:5]
    if (len(sys.argv) == 6):
        minPercentIdentity = float(sys.argv[5])
        
    # require that complete paths are provided for each of the inputs

    currdir = os.getcwd();

    if reads_file[0] != '/':
        reads_file = currdir + "/" + reads_file

    if qual_file[0] != '/':
        qual_file = currdir + "/" + qual_file
        
    if pairs_file[0] != '/':
        pairs_file = currdir + "/" + pairs_file
        
    if reference_db[0] != '/':
        reference_db = currdir + "/" + reference_db
    
    read_pairs = retrieve_read_pairs(pairs_file)

    base = os.path.basename(reads_file)
    
    
    asm_seqs_FH = open(base + ".asm.seqs", 'w')
    asm_quals_FH = open(base + ".asm.quals", 'w')
    asm_log_FH = open(base + ".asm.log", 'w')
    
    scaffold_counter = 0
        
    for read_cluster in read_pairs:
        try:
            scaffold_counter += 1
            reads_text_line = "## Ref16S Assisted Assembly: " + " ".join(read_cluster) + "\n"
            asm_log_FH.write("\n\n\n----------------------------------------------------------\n" \
                             + reads_text_line)
            
            assemble_reads(scaffold_counter, read_cluster, reads_file, qual_file, reference_db,
                           asm_seqs_FH, asm_quals_FH, asm_log_FH)
            
        except (CmdRunnerException), e:
            err_msg = "Error processing read cluster: " + ",".join(read_cluster) + "\n" + str(e)
            sys.stderr.write(err_msg)
            asm_log_FH.write(err_msg)
            traceback.print_exc()

        except (AssembleReadsException):
            asm_log_FH.write("Error, Read assembly failed.")
            traceback.print_exc()
        

    # clean up 
    asm_seqs_FH.close()
    asm_quals_FH.close()
    asm_log_FH.close()

    sys.exit(0)


class AssembleReadsException(Exception): pass


def assemble_reads(scaffold_counter, read_cluster, reads_file, qual_file, reference_db,
                   asm_seqs_FH, asm_quals_FH, asm_log_FH):

    cwd = os.getcwd()
    
    pid = os.getpid()
    tmpdir = cwd + "/tmp_" + str(pid)
    os.mkdir(tmpdir)
    
    os.chdir(tmpdir)

    try:

        # pull out the reads and quals for read pairs
        local_reads_file = "reads.seq"
        local_qual_file = "reads.qual"
        
        readsFH = open(local_reads_file, 'w')
        qualsFH = open(local_qual_file, 'w')
        
        for read in read_cluster:
            read_fasta = cdbyank(read, reads_file)
            read_fasta = read_fasta.replace('\t', ' ')
            readsFH.write(read_fasta)
            
            read_qual = cdbyank(read, qual_file)
            qualsFH.write(read_qual)

        readsFH.close()
        qualsFH.close()

        reference_accession = find_best_matching_reference_sequence(local_reads_file, reference_db, asm_log_FH)
        if reference_accession is None:
            asm_log_FH.write("No reference match for reads: " + ",".join(read_cluster) + "\n")
            raise AssembleReadsException()
        
        # write the reference sequence to a file.
        local_refseq_file = "ref.seq"
        refFH = open(local_refseq_file, 'w')
        ref_fasta = cdbyank(reference_accession, reference_db)
        refFH.write(ref_fasta)
        refFH.close()

        # Run Lucy
        lucy_cmd = " ".join(["lucy -output lucy.seq lucy.qual", local_reads_file, local_qual_file])
        asm_log_FH.write("Running Lucy: " + lucy_cmd + "\n")
        process_cmd(lucy_cmd)
        
        # Run Homology Trimmer
        seq_accessions_processed = homology_trim_lucy_output.homologyTrimLucy("lucy.seq", "lucy.qual", local_refseq_file, asm_log_FH, minPercentIdentity)
    
        if len(seq_accessions_processed) == 0:
            errmsg = "\n\n\tWARNING: ** No sequences passed the homology trim filter. **\n\n"
            sys.stderr.write(errmsg)
            asm_log_FH.write(errmsg)
            raise AssembleReadsException()
        else:
            asm_log_FH.write("Accessions homology and lucy trimmed: " + ",".join(seq_accessions_processed) + "\n")

        
        # Run AMOScmp
        cmd = "tarchive2amos -o raa lucy.seq.htrim.seq"
        asm_log_FH.write("AMOScmp: " + cmd + "\n")
        process_cmd(cmd)

        cmd = "AMOScmp -D TGT=raa.afg -D REF=ref.seq raa"
        asm_log_FH.write("AMOScmp: " + cmd + "\n")
        process_cmd(cmd)
        
        
        # write outputs
        
        # process the contig info
        fh = open("raa.contig_layout", 'r')
        contig_layout = fh.read()
        asm_log_FH.write("ASMcmp contig layout\n" + contig_layout + "\n")
        fh.close()

        placed_reads = _parse_placed_reads_from_contig_layout(contig_layout)

        if len(placed_reads) == 0:
            asm_log_FH.write("* no reads were placed by AMOScmp\n\n")
            raise AssembleReadsException()
        else:
            asm_log_FH.write("The following reads were placed or assembled by AMOScmp: " + ",".join(placed_reads) + "\n")


        # build the scaffold sequence and quals
        (parsed_scaff_acc, ordered_contigs) = contig_layout_to_scaffold.build_ordered_contigs("raa.scaffold_layout")

        scaff_seq = contig_layout_to_scaffold.contigs_to_scaffold_sequence(ordered_contigs, "raa.fasta")
        scaff_quals = contig_layout_to_scaffold.contigs_to_scaffold_quals(ordered_contigs, "raa.qual")

        # check for matching lengths between scaff sequence and lines
        scaffold_seq_length = len(scaff_seq)
        qual_vals = scaff_quals.split(" ")

        if (scaffold_seq_length != len(qual_vals)):
            err_msg = "Error, qual and sequence length do not match:" + \
            "Saffold length: " + str(scaffold_seq_length) + \
            " Number of qual values: " + str(len(qual_vals)) + "\n\n"
            sys.stderr.write(err_msg)
            asm_log_FH.write(err_msg)
            raise AssembleReadsException()
                
        
        asm_seqs_FH.write(">scaff_" + str(scaffold_counter) + " " + ",".join(placed_reads) + "\n" + scaff_seq + "\n")

        # capture the scaffold layout information
        fh = open("raa.scaffold_layout", 'r')
        asm_log_FH.write("AMOScmp scaffolding information:\n" + fh.read() + "\n")
        fh.close()
        
        asm_quals_FH.write(">scaff_" + str(scaffold_counter) + "\n" + scaff_quals + "\n")
        

    finally:
        # go back to CWD and delete sandbox
        os.chdir(cwd); #sys.exit(1) # set exit here to debug (yes, will make option at some point)
        try:
            if not DEBUG:
                shutil.rmtree(tmpdir)
                
        except Exception, e:
            print e
    return



def find_best_matching_reference_sequence(query_file, reference_database, asm_log_FH):

    # check to see that reference database is properly formatted for blast
    if not os.path.isfile(reference_database + ".nin"):
        cmd = "formatdb -i " + reference_database + " -p F"
        process_cmd(cmd)
    

    # run megablast
    cmd = " ".join(["megablast -i", query_file, "-d", reference_database, "-e 1e-20 -v 1 -b 1 -m 8"])
    asm_log_FH.write("Running megablast to find best matching reference sequence\n" + "CMD: " + cmd + "\n")
    results = process_cmd(cmd)
    asm_log_FH.write(results + "\n")
    
    if len(results) == 0:
        asm_log_FH.write("** no blast hits\n")
        return None
    
    # sum up the bit score for each hit
    reference_matches = {}
    results = results.split("\n")
    for result in results:
        #print result
        if not re.search(r"\w", result):
            continue
        x = result.split("\t")
        hit = x[1]
        score = int(x[3])

        if reference_matches.has_key(hit):
            reference_matches[hit] += score
        else:
            reference_matches[hit] = score
        

    hits = reference_matches.items()
    
    hits.sort(hit_sorter)

    # dump hit scores to log
    asm_log_FH.write("Sorted scores for top hits:\n");
    for hit in hits:
        asm_log_FH.write( "\t".join([str(x) for x in hit])  + "\n");
    

    top_hit = hits.pop()

    top_hit_acc = top_hit[0]
    top_hit_acc = re.sub("\#.*$", "", top_hit_acc)  # megablast concatenates the accession with the second field when tab-delimited!

    top_hit_report = "\nTop hit: " + top_hit_acc + " with cum-score: " + `top_hit[1]` + "\n"
    asm_log_FH.write(top_hit_report)
    print top_hit_report
    
    return(top_hit_acc)



def hit_sorter(a,b):
    if(a[1] < b[1]):
        return(-1)
    elif(a[1] == b[1]):
        return(0)
    else:
        return(1)


    

class BlastException(Exception): pass



def retrieve_read_pairs(file):
    fh = open(file, 'r')
    data = fh.read()

    clusters = []
    lines = data.split("\n")
    for line in lines:
        if not re.search(r"\w", line):
            continue
        reads = re.split(r"\s+", line)
        clusters.append(reads)

    return(clusters)



def _parse_placed_reads_from_contig_layout (contig_layout):

    placed_reads = []
    
    for line in contig_layout.split("\n"):
        match = re.search(r"^\#([^\#][^\(\s]+)", line)
        if match:
            acc = match.group(1)
            placed_reads.append(acc)

    return(placed_reads)




main()
