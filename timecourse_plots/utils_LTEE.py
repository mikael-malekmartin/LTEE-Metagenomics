from Bio import SeqIO
from Bio.Seq import Seq

base_table = {'A':'T',
              'T':'A',
              'G':'C',
              'C':'G'
             }

genetic_code = {'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
                'AGA': 'R', 'AGG': 'R', 'AAT': 'N', 'AAC': 'N', 'GAT': 'D', 'GAC': 'D', 'TGT': 'C', 'TGC': 'D',
                'CAA': 'Q', 'CAG': 'Q', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
                'CAT': 'H', 'CAC': 'H', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L',
                'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'AAA': 'K', 'AAG': 'K', 'ATG': 'M', 'TTT': 'F', 'TTC': 'F',
                'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
                'AGT': 'S', 'AGC': 'S', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'TGG': 'W', 'TAT': 'Y',
                'TAC': 'Y', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'TAA': '!', 'TGA': '!', 'TAG': '!' 
               }

def genbank(pathToGenbankFile: str):
    '''
    Convert Genbank file to object
    
    :param pathToGenbankFile: Path to Genbank file
    :type pathToGenbankFile: str
    :return: Object BioPython
    :rtype: <class 'Bio.SeqRecord.SeqRecord'>
    '''

    with open(pathToGenbankFile, mode='r') as f:
        genbank = SeqIO.read(f, 'gb')
    return genbank

def cds_features(genbank) -> list:
    '''
    Fetch and store CDS information from a genbank file object in a list
    
    :param genbank: A BioPython object from genbank file
    :type genbank: <class 'Bio.SeqRecord.SeqRecord'>
    :return: List of CDS
    :rtype: list
    '''

    cds = []
    for feature in genbank.features:
        if feature.type == 'CDS':
            cds.append(feature)

    return cds

def geneProfile(genbankRef, cds_list: list, geneOfInterest:str):
    '''
    Fetch information of Gene of interest and return sequence, codons and strand in genome
    
    :param genbankRef:  BioPython object from genbank file
    :type genbankRef: <class 'Bio.SeqRecord.SeqRecord'>
    :param cds_list: List of CDS previously generate
    :type cds_list: list
    :param geneOfInterest: Gene of interest name
    :type geneOfInterest: str
    :return: return sequence, codons and strand in genome
    :rtype: 
    '''

    geneList = []
    for gene in cds_list :
        if 'gene' in gene.qualifiers.keys():
            if gene.qualifiers['gene'][0].casefold() == geneOfInterest.casefold() :
                geneList.append(gene)
                strand = gene.location.strand
    
    gene_start = list(geneList[0].location)[0]
    gene_end = list(geneList[0].location)[-1] + 1
    
    seq_ref = genbankRef.seq
    
    n = 3
    
    if strand < 0:
        gene_seq_ref = seq_ref[gene_end - 1:gene_start + 1]
    else:
        gene_seq_ref = seq_ref[gene_start:gene_end]
    
    ref_codons = [(gene_seq_ref[i:i+n]) for i in range(0, len(gene_seq_ref), n)]
        
    return geneList, seq_ref, ref_codons, strand

def geneModification(geneList: list, seq_ref, ref_codons: list, mutation: str, position: int, strand: int):
    '''
    Use Genome reference and mutation information (nucleotide, position, strand) to predict amino acid changing
    
    :param geneList: Gene details
    :type geneList: list
    :param seq_ref: Sequence from ancestor or reference
    :type seq_ref: BioPthon sequence object
    :param ref_codons: Codons list from ancestor or reference
    :type ref_codons: list
    :param mutation: The mutated nucleotide (A, T, G, C)
    :type mutation: str
    :param position: Mutation position in genome
    :type position: int
    :param strand: The strand where gene is located in geneome (-/+)
    :type strand: int
    :return: Amino acid in ancestor/reference and in mutated genome
    :rtype 

    '''

    gene_start = list(geneList[0].location)[0]
    gene_end = list(geneList[0].location)[-1] + 1

    seq_mut = list(seq_ref)

    if strand < 0:
        seq_mut[position - 2] = base_table[mutation][0]
    else:
        seq_mut[position - 1] = mutation
    
    n = 3
    
    if strand < 0:
        gene_seq_mut = Seq(''.join(seq_mut))[gene_end - 1:gene_start + 1]
    else: 
        gene_seq_mut = Seq(''.join(seq_mut))[gene_start:gene_end]
    
    mut_codons = [(gene_seq_mut[i:i+n]) for i in range(0, len(gene_seq_mut), n)]
    
    if strand < 0:
        aa_ref = genetic_code[ref_codons[(position - gene_end - 1) // 3].reverse_complement()]
        aa_mut = genetic_code[mut_codons[(position - gene_end - 1) // 3].reverse_complement()]
    else:
        aa_ref = genetic_code[ref_codons[(position - gene_start - 1) // 3]]
        aa_mut = genetic_code[mut_codons[(position - gene_start - 1) // 3]]

    return aa_ref, aa_mut

if __name__ == "__main__":
    help(genbank)