#!/usr/bin/env python3
import os
import argparse
import re
from collections import defaultdict

# Adding description and example usage to the help message
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description="""Convert gene count to RPKM or TPM.
    
    By default, the input file should have the first column with Gene Name/ID, 
    the last column with gene length, and other columns with counts, separated by tabs.
    If -l is specified, the input file shouldn't have gene length in the last column,
    instead lengths from the specified file will be used. Only genes contained in the 
    length file will be considered.
    
    Sample names in the first row should contain _R1/_R2/_R3, etc, identifying its replicate number.
    By default, both sample-based and group-based results are generated. If sample name includes 
    _R1/_R2 etc. at the end, samples will be merged to generate group-based results.
    
    Example usage:
    python count2rpkm_tpm.py -i input_file.txt -l length_file.txt -m RPKM
    python count2rpkm_tpm.py -i input_file.txt -l length_file.txt -m TPM
    """
)

parser.add_argument("-i", "--input", required=True, help="File including counts, with a header indicating sample names.", metavar="input file")
parser.add_argument("-l", "--length_file", default=None, help="File including gene ID and their lengths.", metavar="length file")
parser.add_argument("-m", "--mode", required=True, choices=['RPKM', 'TPM'], help="Mode of conversion: RPKM or TPM.")
parser.add_argument("-s", "--sample_mode", action='store_true', help="If specified, only sample-based result will be generated.")
parser.add_argument("-g", "--group_mode", action='store_true', help="If specified, only group-based result will be generated.")
parser.add_argument("-o", "--output_prefix", default="", help="Output file prefix. Default output file names are TPM.sample.txt and TPM.group.txt.")
parser.add_argument("-d", "--output_directory", default=None, help="Output directory. Defaults to the directory of the input file.")

args = parser.parse_args()
input_file = args.input
length_file = args.length_file
mode = args.mode
sample_mode = args.sample_mode
group_mode = args.group_mode
output_prefix = args.output_prefix
output_directory = args.output_directory

if output_directory:
    output_directory = output_directory
else:
    output_directory = os.path.dirname(input_file)
if output_directory == "":
    output_directory = "."
output_directory += "/"

if output_prefix:
    output_prefix += "."

assert (sample_mode and group_mode) == False, "-s and -g cannot be both specified!"

if sample_mode:
    print("Run on sample mode, only generating sample-based result")
if group_mode:
    print("Run on group mode, only generating group-based result")

names = locals()
d_len = dict()
input_genes = list()
length_genes = set()
final_genes = list()

if length_file:
    with open(length_file) as f:
        lines = f.readlines()
        if lines[0].split("\t")[1].strip().isnumeric():
            lines = lines
        else:
            lines = lines[1:]
        for line in lines:
            gene, length = line.strip().split("\t")
            length_genes.add(gene)
            d_len[gene] = float(length)

    with open(input_file) as f:
        header = f.readline()
        samples = header.strip().split("\t")[1:]
        groups = set()
        for sample in samples:
            if re.search(r"(.*?)_R\d+$", sample):
                groups.add(re.search(r"(.*?)_R\d+$", sample).group(1))
            else:
                groups.add(sample)
        for sample in samples:
            names["d_%s" % sample] = defaultdict(float)
        for group in groups:
            names["d_%s" % group] = defaultdict(float)
        lines = f.readlines()
        for line in lines:
            line_split = line.strip().split("\t")
            gene = line_split[0]
            input_genes.append(gene)
            res = line_split[1:]
            for i, v in zip(samples, res):
                try:
                    names["d_%s" % i][gene] = float(v)
                except ValueError:
                    print(line)
                if re.search(r"(.*?)_R\d+$", i):
                    names["d_%s" % re.search(r"(.*?)_R\d+$", i).group(1)][gene] += float(v)

    final_genes = [g for g in input_genes if g in length_genes]
    final_genes_set = set(final_genes)

else:
    with open(input_file) as f:
        header = f.readline()
        samples = header.strip().split("\t")[1:-1]
        groups = set()
        for sample in samples:
            if re.search(r"(.*?)_R\d+$", sample):
                groups.add(re.search(r"(.*?)_R\d+$", sample).group(1))
            else:
                groups.add(sample)
        for sample in samples:
            names["d_%s" % sample] = defaultdict(float)
        for group in groups:
            names["d_%s" % group] = defaultdict(float)
        lines = f.readlines()
        for line in lines:
            line_split = line.strip().split("\t")
            gene = line_split[0]
            input_genes.append(gene)
            res = line.strip().split("\t")[1:-1]
            for i, v in zip(samples, res):
                try:
                    names["d_%s" % i][gene] = float(v)
                except ValueError:
                    print(line)
                if re.search(r"(.*?)_R\d+$", i):
                    names["d_%s" % re.search(r"(.*?)_R\d+$", i).group(1)][gene] += float(v)
            d_len[gene] = float(line_split[-1])

    final_genes = input_genes

def sample_prepare():
    for sample in samples:
        names["sample_READS_%s" % sample] = defaultdict(float)
        names["sample_all_READS_%s" % sample] = 0
        for gene in final_genes:
            names["sample_READS_%s" % sample][gene] = names["d_%s" % sample][gene]
            names["sample_all_READS_%s" % sample] += names["d_%s" % sample][gene]
    for sample in samples:
        print(sample + "\t" + "all_READS\t%s" % names["sample_all_READS_%s" % sample])

def group_prepare():
    for group in groups:
        names["group_READS_%s" % group] = defaultdict(float)
        names["group_all_READS_%s" % group] = 0
        for gene in final_genes:
            names["group_READS_%s" % group][gene] = names["d_%s" % group][gene]
            names["group_all_READS_%s" % group] += names["d_%s" % group][gene]
    for group in groups:
        print(group + "\t" + "all_READS\t%s" % names["group_all_READS_%s" % group])

def sample_out_rpkm():
    with open("%sRPKM.sample.txt" % (output_directory + output_prefix), "w") as f:
        f.write("GID\t" + "\t".join(sorted(samples)) + "\n")
        for gene in final_genes:
            f.write(gene)
            for sample in sorted(samples):
                f.write("\t" + str(names["sample_READS_%s" % sample][gene] / (d_len[gene] / 1000) / names["sample_all_READS_%s" % sample] * 1000000))
            f.write("\n")

def group_out_rpkm():
    with open("%sRPKM.group.txt" % (output_directory + output_prefix), "w") as f:
        f.write("GID\t" + "\t".join(sorted(groups)) + "\n")
        for gene in final_genes:
            f.write(gene)
            for group in sorted(groups):
                f.write("\t" + str(names["group_READS_%s" % group][gene] / (d_len[gene] / 1000) / names["group_all_READS_%s" % group] * 1000000))
            f.write("\n")

def sample_out_tpm():
    for sample in samples:
        names["sample_RPK_%s" % sample] = defaultdict(float)
        names["sample_all_RPK_%s" % sample] = 0
        for gene in final_genes:
            names["sample_RPK_%s" % sample][gene] = names["d_%s" % sample][gene] / (d_len[gene] / 1000)
            names["sample_all_RPK_%s" % sample] += names["sample_RPK_%s" % sample][gene]
    with open("%sTPM.sample.txt" % (output_directory + output_prefix), "w") as f:
        f.write("GID\t" + "\t".join(sorted(samples)) + "\n")
        for gene in final_genes:
            f.write(gene)
            for sample in sorted(samples):
                f.write("\t" + str(names["sample_RPK_%s" % sample][gene] / names["sample_all_RPK_%s" % sample] * 1000000))
            f.write("\n")

def group_out_tpm():
    for group in groups:
        names["group_RPK_%s" % group] = defaultdict(float)
        names["group_all_RPK_%s" % group] = 0
        for gene in final_genes:
            names["group_RPK_%s" % group][gene] = names["d_%s" % group][gene] / (d_len[gene] / 1000)
            names["group_all_RPK_%s" % group] += names["group_RPK_%s" % group][gene]
    with open("%sTPM.group.txt" % (output_directory + output_prefix), "w") as f:
        f.write("GID\t" + "\t".join(sorted(groups)) + "\n")
        for gene in final_genes:
            f.write(gene)
            for group in sorted(groups):
                f.write("\t" + str(names["group_RPK_%s" % group][gene] / names["group_all_RPK_%s" % group] * 1000000))
            f.write("\n")


if mode == 'RPKM':
    if sample_mode:
        sample_prepare()
        sample_out_rpkm()
    elif group_mode:
        group_prepare()
        group_out_rpkm()
    else:
        sample_prepare()
        group_prepare()
        sample_out_rpkm()
        group_out_rpkm()
elif mode == 'TPM':
    if sample_mode:
        sample_prepare()
        sample_out_tpm()
    elif group_mode:
        group_prepare()
        group_out_tpm()
    else:
        sample_prepare()
        group_prepare()
        sample_out_tpm()
        group_out_tpm()
