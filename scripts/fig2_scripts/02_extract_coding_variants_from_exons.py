#!/gpfs/commons/home/jeinson/python3_env/bin/python

import argparse
import pysam
import pandas as pd
import gzip
import sys

def column_sample_map(path, start_col=9, line_key="#CHR"):
    stream_in = gzip.open(path, "r")
    out_map = {}
    for line in stream_in:
        if line_key in line.decode('ASCII'):
            line = line.decode('ASCII').rstrip().split("\t")
            for i in range(len(line)):
                out_map[i] = line[i]
            break
    stream_in.close()
    return out_map

def fun_dict_from_vcf_info(info):
    info = info.split(";")
    dict_info = {}
    for item in info:
        fields = item.split("=")
        if len(fields) == 2:
            tag = fields[0]
            value = fields[1]
            dict_info[tag] = value
    return dict_info

def flush_print(text):
    print(text)
    sys.stdout.flush()

### MAIN SCRIPT
def main():

    # Arguments passed
    parser = argparse.ArgumentParser()
    # required
    parser.add_argument("--vcf_all", help="Phased genotype VCF with all variants")
    parser.add_argument("--cadd", help="CADD annotation file")
    parser.add_argument("--cadd_indel", help="CADD indel annotation file")
    parser.add_argument("--gnomad", help="gnomad file to check rare variant status")
    parser.add_argument("--chrom", help="Chromosome to use")
    parser.add_argument("--exon_bed", help = "A bed file containing coordinates of exons to extract variants from")

    global args
    args = parser.parse_args()

    # First print the header for the output file
    flush_print("\t".join(['indv', 'gene', 'exon_id', 'exon_coord', 'genotype', 'var', 'chr', 'pos', 'ref', 'alt', 'gnomad_af', 'consequence', 'cadd_phred']))

    # Start by importing all exons to a data frame
    exon_coords = pd.read_csv(args.exon_bed, sep = "\t", index_col = False)
    exon_coords_chr = exon_coords[exon_coords['#Chr'] == args.chrom]

    # Make parser objects for the VCF, CADD, and Gnomad file
    gtex_vcf_tabix = pysam.Tabixfile(args.vcf_all, "r")
    gnomad_tabix = pysam.Tabixfile(args.gnomad, "r")
    cadd_tabix = pysam.Tabixfile(args.cadd, "r")
    cadd_indel_tabix = pysam.Tabixfile(args.cadd_indel, "r")

    # Get the column names for the VCF
    gtex_sample_map = column_sample_map(args.vcf_all)

    # Now loop through these  exons
    for index, row in exon_coords_chr.iterrows():

        gene = row['gene_id']
        exon_id = row['exon_id']
        exon_coord = row['ID']

        variants = gtex_vcf_tabix.fetch(row['#Chr'], row['start'], row['end'])
        for var in variants:
            var = var.split("\t")

            chrom = var[0]
            pos = var[1]
            ref = var[3]
            alt = var[4]

            info_dict = fun_dict_from_vcf_info(var[7])

            # Check the gnomad allele frequency. If it's rare enough, continue
            gnom = gnomad_tabix.fetch(chrom, int(pos) - 1, int(pos))
            gnomad_af = 'NA'
            for record in gnom:
                record = record.split("\t")
                
                if(ref == record[3] and alt == record[4]):
                    gnomad_af = fun_dict_from_vcf_info(record[7])['AF']
                               
            # If the gnomad AF is still NA, it means this variant is missing from gnomad
            # (which I don't think would ever happen since GTEx is part of gnomad)
            if(gnomad_af != 'NA'):
                if(float(gnomad_af) > .01):
                    continue

            # Grab the consequence and CADD score for the variant
            if(len(ref) > 1 or len(alt) > 1): # If the variant is an inddel
                # This will only stay when the indel isn't in CADD
                cadd_phred = "NA"
                consequence = "indel_nos"

                cadd = cadd_indel_tabix.fetch(chrom.replace("chr", ""), int(pos) - 1, int(pos))
                for record in cadd:
                    record = record.split("\t")
                    if(ref == record[2] and alt == record[3] and gene == record[18]):
                        cadd_phred = record[-1]
                        consequence = record[9]
                        break
               
            
            else: # If the variant is a SNP
                cadd = cadd_tabix.fetch(chrom.replace("chr", ""), int(pos) - 1, int(pos))
                for record in cadd:
                    record = record.split("\t")
                    if(ref == record[2] and alt == record[3] and gene == record[18]):
                        cadd_phred = record[-1]
                        consequence = record[9]
                        continue

            # Now loop through the variants to find the individuals who carry it
            for i, gt in enumerate(var):
                if(i < 9): continue
                if (gt.split(":")[0] not in {'0/0', "./."}):

                    # Finally print the output!!
                    indv = gtex_sample_map[i]
                    genotype = gt.split(":")[0]
                    full_variant = var[2]

                    flush_print("\t".join(map(str, [
                        indv, gene, exon_id, exon_coord, genotype, full_variant, chrom, pos, ref, alt, gnomad_af, consequence, cadd_phred
                    ])))

                #import code; code.interact(local=locals()); sys.exit()


if __name__ == "__main__":
    main()
