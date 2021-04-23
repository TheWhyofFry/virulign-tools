"""

The purpose of this script is to run virulign on input sequences and
deal with any pesky FrameShift failures that occur.

Procedure

1. Run Virulign on query sequences using max of 0 FrameShifts
2. If some make it through
   - recreate a reference from those sequences (assuming they do not start/end with gaps)
   - run virulign on these references 
3. If none make it through/procedure in (2) fails on some:
   - Run sequences through blast database
   - For each match, run virulign on a reference made from the target sequence
   - Do until matches are exhausted/all sequences have been aligned
4. If 3 fails - last ditch effort:
   - Run the sequences through blastx
   - Get all the HSPs/target
   - Highest bit score wins
   - Fill in positionss not covered with "X"
5. If 2-3 fails, for the nucleotide alignment, two options:
   - Just do a normal alignment
   - If some of the sequences made it through, add
     the subsequent sequences to the existing alignment



virulign ref.xml/fasta target --exportWithReference yes --exportAlphabet Nucleotides/Aminoacids --exportKind GlobalAlignment


"""

import os
import shutil
import subprocess
import io
import tempfile
import pandas as pd

def read_fasta(f,upper=True):
    #fasta_dict = {}
    fasta_list = []
    curseq = ""
    currentkey = ""
    fastafile = open(f, "r") if type(f) is str else f
    for line in fastafile:
        if line.startswith(">"):
            if curseq != "":
                #fasta_dict[currentkey] = curseq
                try:
                    fasta_list.append((currentkey,curseq.upper()))
                    currentkey = line[1:].strip()
                    curseq = ""
                except:
                    print(f)
                    print(curseq)
                    raise
            else:
                
                currentkey = line[1:].strip()
        else:
            curseq += line.strip().replace(".","-")

    fastafile.close()


    fasta_list.append((currentkey, curseq))
    return fasta_list


def write_fasta(fasta_list, output_file="string"):

    output = "\n".join([">%s\n%s\n"%(title,seq) for title, seq in fasta_list])

    if output_file == "string":
        return output
    else:
        with open(output_file,"w") as outfile:
            outfile.write(output)
        return True

def run_output(command):
    print(command)
    command_split = list(map(lambda x:x.replace("#", " "),command.split(" ")))
    proc = subprocess.Popen(command_split, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    stdout = stdout.decode('utf-8').replace("\\n","\n")
    stderr = stderr.decode('utf-8').replace("\\n","\n")
    returncode = proc.returncode
    
    return stdout, stderr, returncode

def parse_virulign(stdout):
    
    with io.StringIO(stdout) as virulign_output:
        fasta_entries = read_fasta(virulign_output)

    return fasta_entries

 
# Fasta is a fasta list produced by read_fasta
def find_missing_entries(query_fasta, virulign_fasta):

    query_entries = [title for title, seq in query_fasta]
    virulign_entries = [title for title, seq in virulign_fasta]

    
    missing_entries =  set(query_entries).difference(virulign_entries)
    print("MISSING ENTRIES:", missing_entries, len(query_fasta), len(virulign_fasta))
    return [(title, seq) for title,seq in query_fasta if title in missing_entries]


def blast(missing_entries, blast_db, blast_prog="blastn", other_options=""):


    blast_command = "{blast_prog} -query {query} -db {blast_db} -outfmt '6 qseqid sseqid pident length mismatch qstart qend sstart ssend bitscore qseq sseq' {other_options}"



    stdout, stderr, returncode = run_output(blast_command)



def initial_alignment(fasta_file, virulign_ref_file, virulign_command="virulign",
                      alphabet="AminoAcids", exportkind="GlobalAlignment",
                      exportwithref="yes", maxFrameShifts=3):
    
    virulign_command_template = "{v} {virulign_ref_file} {fasta_file} --exportKind {exportkind} --exportAlphabet {alphabet} --exportReferenceSequence {ref}" + \
                                   " --maxFrameShifts {maxFrameShifts}"

        
    stdout, stderr, returncode = run_output(virulign_command_template.format(v=virulign_command,
                                                                     virulign_ref_file=virulign_ref_file,
                                                                     fasta_file=fasta_file,
                                                                     exportkind=exportkind,
                                                                     alphabet=alphabet,
                                                                     ref=exportwithref,
                                                                     maxFrameShifts=maxFrameShifts))
    # Fix that quirky case when arbitrary characters are written out by virulignmp
    if ">" in stdout:
        #stdout = stdout[stdout.index(">"):]

        with io.StringIO(stdout) as fasta_output:
            fasta_list = read_fasta(fasta_output)
    else:
        fasta_list = []
    return fasta_list


        
def do_secondary(missing_entries, ref_entry, virulign_params={}):
    
    code, temp_ref_file = tempfile.mkstemp(prefix="virulign-tools", suffix=".fasta", dir="/tmp")
    code, temp_missing_file = tempfile.mkstemp(prefix="virulign-tools", suffix=".fasta", dir="/tmp")
    print("TEMP MISSING FILE:", temp_missing_file)
    write_fasta(missing_entries, temp_missing_file)

    with open(temp_ref_file, 'w') as temp_ref:
        temp_ref.write(">%s\n%s\n"%(ref_entry))

    
    # We'll just directly assume frameshifts are OK
    secondary_alignment = initial_alignment(temp_missing_file, temp_ref_file,  **virulign_params)

    secondary_missing_entries = find_missing_entries(missing_entries, secondary_alignment)

    if len(secondary_missing_entries) > 0:
        missing_entries = secondary_missing_entries
    

    #os.unlink(temp_ref_file)
    #os.unlink(temp_missing_file)


    return secondary_alignment, secondary_missing_entries








if __name__ == "__main__":
    import argparse
    import textwrap

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", type=str, dest="input_fasta", help="Input FASTA file")
    parser.add_argument("-n", type=str, dest="output_fasta", help="Output FASTA file (nuc)")
    parser.add_argument("-a", type=str, dest="output_fasta", help="Output FASTA file (aa)")
    parser.add_argument("-r", type=str, dest="virulign_ref_file", help="Initial Virulign reference file")
    parser.add_argument("-b", type=str, dest="blastn", help="BLASTN DB", default="")
    parser.add_argument("-B", type=str, dest="blastp", help="BLASTP DB", default="")
    parser.add_argument("-X", action='store_true', dest="do_blastx", help="If all else fails, align with BLASTX")
    parser.add_argument("-R", action='store_true', dest="do_selfref", help="Align failed sequences to sequences that passed first")
    parser.add_argument("-N", action='store_true', dest="do_blastn", help="Repeat alignment for failed sequences using the specified BLASTN DB")
    parser.add_argument("-s", type=int, dest="num_alignments", help="Number of blast alignments", default=20)

    
    args = parser.parse_args()

    # Intial FASTA

    fasta_initial = read_fasta(args.input_fasta)
    n_query = len(fasta_initial)

    # Initial alignment
    fasta_initial_alignment = initial_alignment(args.input_fasta, args.virulign_ref_file, maxFrameShifts=0, alphabet="Nucleotides")
    fasta_initial_alignment_aa = initial_alignment(args.input_fasta, args.virulign_ref_file, maxFrameShifts=0, alphabet="AminoAcids")
    reference = fasta_initial_alignment[0]
    # Missing alignemnts

    missing_entries = find_missing_entries(fasta_initial, fasta_initial_alignment)
    print([title for title,seq in fasta_initial])
    print([title for title,seq in fasta_initial_alignment])
    temp_alignments = []
    n_missing_entries = len(missing_entries)
    print("Missing entries:", n_missing_entries)

    if n_missing_entries > 0:
        code, missing_entries_file = tempfile.mkstemp(prefix="virulign-tools", suffix=".fasta",dir="/tmp")
        print("MISSING ENTRY FILE:", missing_entries_file, code)
        write_fasta(missing_entries, missing_entries_file)

        blastn_fasta_entries = []

        if (args.do_blastn):


            # Get N hits from a blast db to try and use these seqs as reference
            # Note that these are _valid_ references, i.e. seqs where virulign could resolve
            # A fundamental reference sequence (e.g. HXB2 in HIV-1) with a query
            # and WITHOUT any frame shift compensation

            BLASTN = "blastn -db {db} -query {query} -num_alignments {n} -outfmt 6#qseqid#sseqid#bitscore".format(db=args.blastn,
                                                                                                                    query=missing_entries_file,
                                                                                                                    n=args.num_alignments)
            print(BLASTN)
            blastn_output, stderr, returncode = run_output(BLASTN)
            with io.StringIO(blastn_output) as blastn_output_f:
                blastn_df = pd.read_table(blastn_output_f,names=["qseqid", "sseqid", "bitscore"])
                blast_hits = blastn_df.value_counts("sseqid").sort_values(ascending=False).reset_index().sseqid.values
                print("BLAST HITS:")
                blastn_entries, stderr, returncode = run_output('blastdbcmd -entry "{entry}" -db {db}'.format(entry=','.join(blast_hits),db=args.blastn))
                
                with io.StringIO(blastn_entries) as blast_file:
                    blastn_fasta_entries = read_fasta(blast_file)
                    print("BLASTN Entries:", blastn_fasta_entries)

            



        if (args.do_selfref) and (len(fasta_initial_alignment) > 1):
            fasta_initial_alignment_ref = fasta_initial_alignment[:]
        else:
            fasta_initial_alignment_ref = []

        
        # Placeholder
        rep = True
        if rep:
            print("Finding appropriate references ... ")
            alternative_alignments = fasta_initial_alignment_ref + blastn_fasta_entries
            print([title for title,seq in alternative_alignments])
            for fasta_entry in alternative_alignments:
                title, seq = fasta_entry
                print(title)
                secondary_alignment, missing_entries = do_secondary(missing_entries, fasta_entry, virulign_params=dict(alphabet="AminoAcids", exportwithref="no"))

                if len(secondary_alignment) == 0:
                    continue
                
                temp_alignments.extend(secondary_alignment)

                n_missing_entries -= len(secondary_alignment)

                if n_missing_entries == 0:
                    break

    final_alignments_aa = fasta_initial_alignment_aa + temp_alignments

    print(len(final_alignments_aa))
    write_fasta(final_alignments_aa, "/dev/stderr")
        

                




            



















    


