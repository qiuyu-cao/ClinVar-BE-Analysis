from filter_anno_vcf import feilter_anno_vcf
from stats_rules import nt_match, Genome, if_targetable, \
    search_sgrnas, if_record_editable, get_motif, if_pam_match
import vcf
from pathlib import Path
from loguru import logger



genome = "hg38_ucsc/hg38.fa"
edit_from = "A"
edit_to = "G"
editing_window = [3, 9]
out_prefix = "results/clinvar_20240716"

motif_window=[1,0]
sgrna_len=20
pam_motif="NGG"

# editor paramaters
## ABE8e
abe8e_editing_window = editing_window
abe8e_edit_from = edit_from
abe8e_edit_to = edit_to
## ABE10
abe10_editing_window = editing_window
abe10_edit_from = edit_from
abe10_edit_to = edit_to
abe10_motif_window = [1, 0]
abe10_motif_seq = "YA"
## ABE9
abe9_editing_window = [5,6]
abe9_edit_from = edit_from
abe9_edit_to = edit_to
## pcABE
pcabe_editing_window = editing_window
pcabe_edit_from = edit_from
pcabe_edit_to = edit_to

from filter_anno_vcf import feilter_anno_vcf
clinvar_vcf = "clinvar_20240716.vcf"
out_prefix = "results/clinvar_20240716"
vcf_file = "clinvar_20240716_filtered_ann_chr.vcf"
if not Path(vcf_file).exists():
    feilter_anno_vcf(clinvar_vcf, vcf_file)
    

def main(genome, vcf_file, edit_from, edit_to, editing_window, out_prefix, sgrna_len, pam_motif):

    genome_reader = Genome(genome)
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    window_size = editing_window[1] - editing_window[0] + 1
    pam_len = len(pam_motif)
    out_title = ["CHROM", "POS", "REF", "ALT", "CLNHGVS", "GeneName", "FeatureType", "HGVS", "Strand"]
    out_title2 = ["Motif", "ABE8e_NGG", "ABE8e_NG", "ABE10_NGG", "ABE10_NG", "ABE9_NGG", "ABE9_NG",
                  "ABE8e_NGG_labels", "ABE10_NGG_labels", "ABE9_NGG_labels", "IsNGGMatch", "sgrna_ngg_list",
                  "ABE8e_NG_labels", "ABE10_NG_labels", "ABE9_NG_labels", "IsNGMatch", "sgrna_ng_list"]
    # all list
    out_path = out_prefix + "_all_targetable.xls"
    out_open = open(out_path, 'w')
    out_open.write('\t'.join(out_title + out_title2) + '\n')
    logger.info('Start to filter the VCF file')
    for record in vcf_reader:
        # if targetable
        strand = if_targetable(record, edit_from, edit_to)
        if strand == 0:
            continue
        # get sgrna
        sgrna_ngg_list = search_sgrnas(genome_reader, record, strand,sgrna_len, 3, edit_from)
        sgrna_ng_list = search_sgrnas(genome_reader, record, strand, sgrna_len, 2, edit_from)
        #
        clnhgvs = ','.join(record.INFO['CLNHGVS'])
        #
        ann_list = record.INFO['ANN'][0].split('|')
        gene_name = ann_list[3]
        feature_type = ann_list[5]
        hgvs = '{}({}):{}({})'.format(ann_list[6], ann_list[3], ann_list[9], ann_list[10]) if ann_list[10] \
                            else '{}({}):{}'.format(ann_list[6], ann_list[3], ann_list[9])
        # ABE8e NGG
        if_abe8e_ngg_editable, abe8e_ngg_labels = if_record_editable(sgrna_ngg_list, sgrna_len, "NGG", edit_from, 
                                                             edit_to, editing_window, motif_window=None, motif_seq=None)
        # ABE8e NG
        if_abe8e_ng_editable, abe8e_ng_labels = if_record_editable(sgrna_ng_list, sgrna_len, "NG", edit_from, 
                                                             edit_to, editing_window, motif_window=None, motif_seq=None)
        # ABE10 NGG
        if_abe10_ngg_editable, abe10_ngg_labels = if_record_editable(sgrna_ngg_list, sgrna_len, "NGG", edit_from, 
                                                             edit_to, editing_window, motif_window=abe10_motif_window, motif_seq=abe10_motif_seq)
        # ABE10 NG 
        if_abe10_ng_editable, abe10_ng_labels = if_record_editable(sgrna_ng_list, sgrna_len, "NG", edit_from, 
                                                             edit_to, editing_window, motif_window=abe10_motif_window, motif_seq=abe10_motif_seq)
        # ABE10 motif
        abe10_sgrna_motif_seq = get_motif(genome_reader, record, abe10_motif_window, strand, edit_from)
        # ABE9 NGG
        if_abe9_ngg_editable, abe9_ngg_labels = if_record_editable(sgrna_ngg_list, sgrna_len, "NGG", edit_from, 
                                                             edit_to, [5, 6], motif_window=None, motif_seq=None)
        # ABE9 NG
        if_abe9_ng_editable, abe9_ng_labels = if_record_editable(sgrna_ng_list, sgrna_len, "NG", edit_from, 
                                                             edit_to, [5, 6], motif_window=None, motif_seq=None)
        # if pam match
        ## if ngg mathc
        is_ngg_match = if_pam_match(sgrna_ngg_list, sgrna_len, "NGG", edit_from, edit_to, editing_window)
        ## if ng match
        is_ng_match = if_pam_match(sgrna_ng_list, sgrna_len, "NG", edit_from, edit_to, editing_window)
        
        line = [record.CHROM, record.POS, record.REF, record.ALT[0], clnhgvs, gene_name, feature_type, hgvs, strand]

        line2 = [abe10_sgrna_motif_seq, 
                 if_abe8e_ngg_editable, if_abe8e_ng_editable, if_abe10_ngg_editable, 
                 if_abe10_ng_editable, if_abe9_ngg_editable, if_abe9_ng_editable,
                 ",".join(abe8e_ngg_labels), ",".join(abe10_ngg_labels), ",".join(abe9_ngg_labels), is_ngg_match,",".join(sgrna_ngg_list),
                 ",".join(abe8e_ng_labels), ",".join(abe10_ng_labels), ",".join(abe9_ng_labels), is_ng_match,",".join(sgrna_ng_list),]
        out_open.write('\t'.join(map(str, line + line2)) + '\n')
    out_open.close()

main(genome, vcf_file, edit_from, edit_to, editing_window, out_prefix, sgrna_len, pam_motif)