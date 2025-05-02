import vcf
from pathlib import Path
import os
import time
import re


def filter_vcf(vcf_file, vcf_filtered_file):
    vcf_reader = vcf.Reader(open(vcf_file))
    vcf_writer = vcf.Writer(open(vcf_filtered_file, 'w'), vcf_reader)
    # vcf_writer._write_header()
    chroms = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']
    for record_ in vcf_reader:
        # if the var type is not snp, then skip
        if not record_.is_snp :
            continue
        # if the ALT is  N or more than one, then skip
        if record_.ALT[0] == 'N' or len(record_.ALT) != 1:
            continue
        # if the record is not in the chroms, then skip
        if record_.CHROM not in chroms:
            continue
        # if the redord.INFO not have CLNSIG, then skip
        if 'CLNSIG' not in record_.INFO:
            continue
        # if the CLNSIG contain the 'Pathogenic', then write the record
        clnsig_ = record_.INFO['CLNSIG']
        clnsig2_ = []
        [clnsig2_.extend(re.split('[|/,]', n)) for n in clnsig_]
        if 'Pathogenic' in clnsig2_:
            vcf_writer.write_record(record_)

def snpeff_ann(vcf_filtered_file, vcf_filtered_ann_file):
    cmd = 'snpEff ann  -canon hg38 {} > {}'.format(vcf_filtered_file, vcf_filtered_ann_file)
    os.system(cmd)
    
def add_chr(vcf_filtered_ann_file, vcf_filtered_ann_chr_file):
    cmd = 'sed \'/^[^#]/s/^/chr/\' {} >  {}'.format(vcf_filtered_ann_file, vcf_filtered_ann_chr_file)
    os.system(cmd)


def feilter_anno_vcf(input_vcf, output_prefix):
    input_vcf = Path(input_vcf)
    # filter the ClinVar VCF
    filtered_vcf = output_prefix + '_filtered.vcf'
    filter_vcf(input_vcf, filtered_vcf)
    # annotate the filtered VCF
    filtered_vcf_ann = output_prefix + '_filtered_ann.vcf'
    snpeff_ann(filtered_vcf, filtered_vcf_ann)
    # add chr to the filtered VCF
    output_vcf = output_prefix + '_filtered_ann_chr.vcf'
    add_chr(filtered_vcf_ann, output_vcf)
    # remove the intermediate files
    # os.remove(filtered_vcf)
    # os.remove(filtered_vcf_ann)
