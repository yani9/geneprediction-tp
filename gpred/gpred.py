import argparse
import sys
import os
import csv
import re


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def isdir(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isdir(path):
        if os.path.isfile(path):
            msg = "{0} is a file".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='genome_file', type=isfile, required=True, 
                        help="Complete genome file in fasta format")
    parser.add_argument('-g', dest='min_gene_len', type=int, 
                        default=50, help="Minimum gene length to consider")
    parser.add_argument('-s', dest='max_shine_dalgarno_distance', type=int, 
                        default=16, help="Maximum distance from start codon "
                        "where to look for a Shine-Dalgarno motif")
    parser.add_argument('-d', dest='min_gap', type=int, default=40,
                        help="Minimum gap between two genes (shine box not included).")
    parser.add_argument('-p', dest='predicted_genes_file', type=str, 
                        default=os.curdir + os.sep +"predict_genes.csv",
                        help="Tabular file giving position of predicted genes")
    parser.add_argument('-o', dest='fasta_file', type=str,
                        default=os.curdir + os.sep + "genes.fna",
                        help="Fasta file giving sequence of predicted genes")
    return parser.parse_args()


def read_fasta(fasta_file):
    """Extract the complete genome sequence as a single string
    """
    seq = []
    with open(fasta_file, 'r') as filin: 
        lines = filin.readlines()
        for i in range(1,len(lines)):
            if not lines[i].startswith(">"):  
                seq.append(lines[i].rstrip().upper())
    return "".join(seq)


def find_start(start_regex, sequence, start, stop):
    """Find the start codon
    """
    sequence = sequence.replace("U","T")
    res_start = start_regex.search(sequence, start, stop)
    if res_start!= None: 
        return res_start.start(0)
    return None


def find_stop(stop_regex, sequence, start):
    """Find the stop codon
    """
    sequence = sequence.replace("U","T")
    iter_stop = stop_regex.finditer(sequence, start, len(sequence))
    for span in iter_stop:
        if span!=None: 
            for stop in span.span(): 
                # Check if start and stop are in same reading frame 
                if ((stop - start) % 3) == 0: 
                    return stop
    return None


def has_shine_dalgarno(shine_regex, sequence, start, max_shine_dalgarno_distance):
    """Find a shine dalgarno motif before the start codon
    """
    shine_pos = shine_regex.search(sequence, start- max_shine_dalgarno_distance)
    if shine_pos!=None and start-shine_pos.end()>6: 
        return True
    return False

def predict_genes(sequence, start_regex, stop_regex, shine_regex, 
                  min_gene_len, max_shine_dalgarno_distance, min_gap):
    """Predict most probable genes
    """
    probable_genes = []
    current_pos = 0
    while (len(sequence)-current_pos >=min_gap):
        current_pos = find_start(start_regex, sequence, current_pos, len(sequence))
        if current_pos!=None: 
            stop = find_stop(stop_regex, sequence, current_pos)
            if stop!=None:
                if stop-current_pos+1>=min_gene_len:
                    has_sd = has_shine_dalgarno(shine_regex, sequence, current_pos , max_shine_dalgarno_distance)
                    if has_sd:
                        probable_genes.append([current_pos+1,stop+3])
                        current_pos = stop + 3 + min_gap
                    else: 
                        current_pos = current_pos + 1
                else: 
                    current_pos = current_pos + 1
            else:
                current_pos = current_pos + 1
    return probable_genes

def write_genes_pos(predicted_genes_file, probable_genes):
    """Write list of gene positions
    """
    try:
        with open(predicted_genes_file, "wt") as predict_genes:
            predict_genes_writer = csv.writer(predict_genes, delimiter=",")
            predict_genes_writer.writerow(["Start", "Stop"])
            predict_genes_writer.writerows(probable_genes)
    except IOError:
        sys.exit("Error cannot open {}".format(predicted_genes_file))


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def write_genes(fasta_file, sequence, probable_genes, sequence_rc, probable_genes_comp):
    """Write gene sequence in fasta format, 
    """
    try:
        with open(fasta_file, "wt") as fasta:
            for i,gene_pos in enumerate(probable_genes):
                fasta.write(">gene_{0}{1}{2}{1}".format(
                    i+1, os.linesep, 
                    fill(sequence[gene_pos[0]-1:gene_pos[1]])))
            i = i+1
            for j,gene_pos in enumerate(probable_genes_comp):
                fasta.write(">gene_{0}{1}{2}{1}".format(
                            i+1+j, os.linesep,
                            fill(sequence_rc[gene_pos[0]-1:gene_pos[1]])))
    except IOError:
        sys.exit("Error cannot open {}".format(fasta_file))


def reverse_complement(kmer):
    """Get the reverse complement"""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in kmer[::-1]])


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Gene detection over genome involves to consider a thymine instead of
    # an uracile that we would find on the expressed RNA
    #start_codons = ['TTG', 'CTG', 'ATT', 'ATG', 'GTG']
    #stop_codons = ['TAA', 'TAG', 'TGA']
    start_regex = re.compile('AT[TG]|[ATCG]TG')
    stop_regex = re.compile('TA[GA]|TGA')
    # Shine AGGAGGUAA
    #AGGA ou GGAGG 
    shine_regex = re.compile('A?G?GAGG|GGAG|GG.{1}GG')
    # Arguments
    args = get_arguments()

    # Let us do magic in 5' to 3'
    sequence = read_fasta(args.genome_file)
    start_regex = re.compile('AT[TG]|[ATCG]TG')
    stop_regex = re.compile('TA[GA]|TGA')
    shine_regex = re.compile('A?G?GAGG|GGAG|GG.{1}GG')
    probable_genes = predict_genes(sequence, start_regex, stop_regex, shine_regex, 50, 16, 40)

    # We reverse and complement
    sequence_rc = reverse_complement(sequence)

    # Call to output functions
    print("reverse sequence : ",sequence_rc)
    res_rc = predict_genes(sequence_rc, start_regex, stop_regex, shine_regex, 50,16,40)
    probable_genes_comp = []
    for i in range(len(res_rc)):
        pos1_53 = len(sequence)-res_rc[i][1]
        pos2_53 = len(sequence)-res_rc[i][0]
        probable_genes_comp.append([pos1_53, pos2_53])
        print(probable_genes_comp)
    print("len predict genes : ", len(probable_genes)+len(probable_genes_comp))

    # Write positions of predicts genes of the two strands 
    write_genes_pos(args.predicted_genes_file, probable_genes+sorted(probable_genes_comp))
    write_genes(args.fasta_file, sequence, probable_genes, sequence_rc, probable_genes_comp)
  

if __name__ == '__main__':
    main()
