import vcf
from pyfaidx import Fasta
from Bio import Seq
from Bio.SeqUtils import nt_search
from pathlib import Path
import re

class Genome():
    def __init__(self, genome:str):
        self.genome = genome
        self.fasta = Fasta(genome)
    def subseq(self, chrom, start, end, reverse=False):
        if not reverse:
            return Seq.Seq(self.fasta[chrom][start-1:end].seq).upper()
        else:
            return Seq.Seq(self.fasta[chrom][start-1:end].seq).reverse_complement().upper()

def nt_match(seq, motif):
    if len(seq) == len(motif) and len(nt_search(str(seq), str(motif))) == 2:
        return True
    else:
        return False

def if_targetable(record, edit_from, edit_to):
    """_summary_

    Args:
        record (vcf.model._Record): _description_
        edit_from (str): _description_
        edit_to (str): _description_
    Returns:
        int: return 1, 0 or -1.
    """
    ref = str(record.REF)
    alt = str(record.ALT[0])
    if nt_match(alt, edit_from) and nt_match(ref, edit_to):
        return 1
    elif nt_match(Seq.Seq(alt).complement(), edit_from) and nt_match(Seq.Seq(ref).complement(), edit_to):
        return -1
    else:
        return 0

def if_multiple_target(Genome, record, window_size, strand, edit_from):
    """_summary_

    Args:
        Genome (Genome): _description_
        record (vcf.model._Record): _description_
        window_size (int): _description_
        strand (int): 1 mean plus strand, -1 means minus strand
    Return:
        sequence(str):
    """
    chrom = record.CHROM
    start = record.POS - window_size + 1
    end = record.POS + window_size - 1
    reverse = True if strand==-1 else False
    flanking_seq = Genome.subseq(chrom, start, end, reverse)
    flanking_seq = f'{flanking_seq[:window_size-1]}{edit_from.lower()}{flanking_seq[window_size:]}'
    multiple_targets = False
    ## 
    for i in range(window_size):
        if flanking_seq[i:i+window_size].count(edit_from) < 2:
            multiple_targets = True
    return multiple_targets, flanking_seq

def get_motif(Genome, record, motif_window, strand, edit_from):
    """_summary_

    Args:
        Genome (Genome): _description_
        record (vcf.model._Record): _description_
        motif (list): _description_
    Return:
        motif_seq (str): seuqnece length = motif[0] + motif[1] + 1
    """
    reverse = True if strand==-1 else False
    chrom = record.CHROM
    start = record.POS - motif_window[0] if not reverse else record.POS - motif_window[1]
    end = record.POS + motif_window[1] if not reverse else record.POS + motif_window[0]
    motif_seq = Genome.subseq(chrom, start, end, reverse)
    return f'{motif_seq[:motif_window[0]]}{edit_from.lower()}{motif_seq[motif_window[0]+1:]}'

def search_sgrnas(Genome, record, strand, sgrna_len, pam_len, edit_from):
    """_summary_

    Args:
        Genome (Genome): _description_
        record (vcf.model._Record): _description_
        sgrna_len (int): _description_
        pam_len (int): _description_
    Return:
        sgrna_list (list), sgRNA (with PAM sequence) list
    """
    reverse = True if strand==-1 else False
    chrom = record.CHROM
    start = record.POS-sgrna_len+1 if not reverse else record.POS-sgrna_len-pam_len+1
    end = record.POS+sgrna_len+pam_len-1 if not reverse else record.POS+sgrna_len-1
    flnaking_seq = Genome.subseq(chrom, start, end, reverse)
    sgrna_list = []
    for i in range(sgrna_len):
        sgrna = f'{flnaking_seq[i:sgrna_len-1]}{edit_from.lower()}{flnaking_seq[sgrna_len:i+sgrna_len+pam_len]}'
        sgrna_list.append(sgrna)
    return sgrna_list

def modify_sgrna_pam(sgrna_pam, sgrna_len, pam_len):
    """_summary_

    Args:
        sgrna_pam (str): guide rna with PAM sequence, the length is sgrna_len + pam_len
        sgrna_len (int): 
        pam_len (int): 
    """
    if len(sgrna_pam) != sgrna_len + pam_len:
        return sgrna_pam
    new_sgrna_pam = sgrna_pam[:sgrna_len] + "|" + sgrna_pam[sgrna_len:]
    return new_sgrna_pam


def if_base_editor_suitable(sgrna_pam, target_position, sgrna_len, pam_motif, edit_from, edit_to, editing_window, 
                            motif_window=None, motif_seq=None):
    labels = []
    # split sgrna as sgrna, pam and window_seq
    sgrna = sgrna_pam[:sgrna_len]
    pam = sgrna_pam[sgrna_len:]
    window_seq = sgrna[editing_window[0] - 1:editing_window[1]]
    # if target nucleotide is in editing window
    if target_position >= editing_window[0] and target_position <= editing_window[1]:
        labels.append("TARGET")
    else:
        labels.append("target")
    # if mutiple target in editing window
    if window_seq.upper().count(edit_from.upper()) > 1:
        labels.append("MULTIPLE")
    else:
        labels.append("multiple")
    # if pam mathch
    if nt_match(pam, pam_motif):
        labels.append("PAM")
    else:
        labels.append("pam")
    # if motif match
    if not motif_window or not motif_seq:
        pass
    else:
        # if motif match
        sg_motif = sgrna_pam[target_position-motif_window[0]-1:target_position+motif_window[1]]
        if nt_match(sg_motif.upper(), motif_seq):
            labels.append("MOTIF")
        else:
            labels.append("motif")
    return labels
    

def if_record_editable(sgrna_pam_list, sgrna_len, pam_motif, edit_from, edit_to, editing_window, 
                       motif_window=None, motif_seq=None):
    is_editable_list = []    
    target_position =  sgrna_len
    labels_list = []
    if not motif_window or not motif_seq:
        for sgrna_pam in sgrna_pam_list:
            labels = if_base_editor_suitable(sgrna_pam, target_position, sgrna_len, pam_motif, edit_from, 
                                             edit_to, editing_window)
            if "TARGET" in labels and "MULTIPLE" not in labels and "PAM" in labels:
                is_editable_list.append(True)
                labels_list.append("True:" + "|".join(labels))
            else:
                is_editable_list.append(False)
                labels_list.append("False:" + "|".join(labels))
            target_position -= 1
    else:
        for sgrna_pam in sgrna_pam_list:
            labels = if_base_editor_suitable(sgrna_pam, target_position, sgrna_len, pam_motif, edit_from, 
                                             edit_to, editing_window, motif_window, motif_seq)
            if "TARGET" in labels and  "PAM" in labels and "MOTIF" in labels:
                is_editable_list.append(True)
                labels_list.append("True:" + "|".join(labels))
            else:
                is_editable_list.append(False)
                labels_list.append("False:" + "|".join(labels))
            target_position -= 1
    is_editable = True if is_editable_list.count(True)>0 else False
    return is_editable, labels_list

def if_pam_match(sgrna_pam_list, sgrna_len, pam_motif, edit_from, edit_to, editing_window):
    is_pam_match_list = []
    target_position =  sgrna_len
    for sgrna_pam in sgrna_pam_list:
        labels = if_base_editor_suitable(sgrna_pam, target_position, sgrna_len, pam_motif, edit_from, 
                                             edit_to, editing_window)

        if "PAM" in labels and "TARGET" in labels:
            is_pam_match_list.append(True)
        target_position -= 1
    is_pam_match = True if is_pam_match_list.count(True)>0 else False
    return is_pam_match
