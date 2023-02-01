#!/gpfs/commons/home/jeinson/python3_env/bin/python
import argparse
import gzip
import subprocess
import sys

import random
import pandas
import csv
import pysam
import time
import numpy as np
from scipy import stats

# This is a fork of Dafni's haplotype_configuration_gtex.py script. It is adapted to work 
# with PSI splicing QTLs. 

# One major difference is that it designates the HIGHER included haplotype as 'a' in the haplotype calls

# To test:
# python scripts/haplotype_configuration_gtex_JE.py --vcf_coding $GTEX_CODING_VCF --vcf_all $GTEX_ALL_VCF --cadd $CADD --genes sQTL/sGenes_by_covariates/fifteen_PEERS_gt_tech_PSI_sQTLs.tsv --chrom 4

### FUNCTIONS
def sample_column_map(path, start_col=9, line_key="#CHR"):
    stream_in = gzip.open(path, "r")
    out_map = {}
    for line in stream_in:
        if line_key in line.decode('ASCII'):
            line = line.decode('ASCII').rstrip().split("\t")
            for i in range(start_col, len(line)):
                out_map[line[i]] = i
            break
    stream_in.close()
    return out_map


def flush_print(text):
    print(text)
    sys.stdout.flush()


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

def is_in_between(mid, lower, upper):
    return mid in range(lower, upper)

### MAIN SCRIPT
def main():
    vep_cons = np.loadtxt("/gpfs/commons/groups/lappalainen_lab/jeinson/data/vep_severity.txt", dtype = "str", encoding = None)[:,0]

    # Arguments passed
    parser = argparse.ArgumentParser()
    # required
    parser.add_argument("--vcf_coding", help="Phased genotype VCF with coding variants")
    parser.add_argument("--vcf_all", help="Phased genotype VCF with all variants")
    parser.add_argument("--cadd", help="CADD annotation file")
    parser.add_argument("--gnomad", help="gnomad file to check rare variant status")
    parser.add_argument("--anc_allele", help="File containing the ancestral alleles")
    parser.add_argument("--genes", help="File containing the sQTL results. It should have a column called 'group' for the genes, and a column called 'top_pid' for the top exon from permuations")
    parser.add_argument("--chrom", help="Chromosome to use")
    parser.add_argument("--indvs", help="A list of individuals to used in the QTL mapping, to recompute the minor allele frequency")
    parser.add_argument("--recompute", default = True, help="Should I requery the VCF? Speeds up runtime for debugging")

    global args
    args = parser.parse_args()

    df_eqtl = pandas.read_csv(args.genes, sep="\t", index_col=False)
    # sGenes used in our analysis
    set_egenes = set(df_eqtl['group']) # changed from gene

    vcf_phased_path = args.vcf_all
    tabix_vcf_phased = pysam.Tabixfile(vcf_phased_path, "r")
    phased_sample_map = sample_column_map(vcf_phased_path)

    vcf_coding_path = args.vcf_coding

    flush_print("\t".join(map(str,
                              ['indv', 'gene', 'exon_coord', 'top_exon_coord', 'chr', 
                              'csnp', 'effect', 'cadd_effect', 'CADD', 'csnp_af',
                              'csnp_af_gnomad', 'csnp_type', 'csnp_ref_allele', 
                              'csnp_alt_allele', 'csnp_anc_allele', 'dummy',
                              'esnp', 'esnp_hi_inc_allele', 'esnp_low_inc_allele', 'esnp_anc_allele', 
                              'slope', 'esnp_af', 'esnp_gnomad_af', 'haplotype'])))
    
    # Open a tabix object for the CADD file
    tabix_cadd = pysam.Tabixfile(args.cadd, "r")

    # Open a tabix object for the gnomad file
    gnomad_cadd = pysam.Tabixfile(args.gnomad, "r")

    # Open a tabix object for the ancestral allele file
    tabix_anc_allele = pysam.Tabixfile(args.anc_allele, "r")
    
    # Make a dictionary that maps genes names to gene IDs
    map_path = "/gpfs/commons/groups/lappalainen_lab/jeinson/data/gene_id_symbol_map.tsv"
    gene_id_map={}
    with open(map_path) as f:
        next(f)
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            gene_id_map[row[1]] = row[0]        

    # Read in the list of individuals used for QTL mapping
    #indvs = [line.rstrip('\n') for line in open(args.indvs)]
    indvs = list(phased_sample_map.keys())
    # Make sure to only use the intersection


    # retrieve specific chromosome
    # Get the name of the tissue being worked on, for saving a unique tmp file
    tiss_long = args.genes.split("/")[-1]

    # Make a new empty vcf file with coding variant coordinates
    recompute = args.recompute
    
    # Define tmp file names outside of the if statement, so we still get them if recompute==False
    random_number = str(random.randint(1,1e15)) # to prevent weirdness when multiple runs happening simultaneously
    coord_fp = "tmp/" + random_number + tiss_long + args.chrom + "_exons_tabix_coords.tsv"
    vcf_fp = tiss_long + args.chrom + random_number + ".vcf"
    unique_vcf_fp = vcf_fp + ".unique"

    if recompute == True:
        subprocess.call("> tmp/" + vcf_fp, shell = True)
        
        # Save only regions that cover top exons! This is a difference from the other script
        # Also make sure we're only using significant sQTLs
        
        # Try an alternative way of doing this
        coord_file = open(coord_fp, "w")
        
        exons = df_eqtl[(df_eqtl.chr == "chr" + args.chrom) & (df_eqtl.padj < .05)].pid
               
        for exon in exons:
            exon_coords = exon.split("_")
            # Get rid of the "chr" when using SSC
            coord_file.write("\t".join([exon_coords[0], exon_coords[1], exon_coords[2]]) + "\n")
        coord_file.close()
        
        # Sort the coord file using subprocess
        sorted_coord_fp = coord_fp + "sort"
        subprocess.call("sort -V " + coord_fp + " | uniq > " + sorted_coord_fp, shell = True)
        
        # Build a VCF by querying based on the coordinates
        subprocess.call("tabix -R " + sorted_coord_fp + " " + vcf_coding_path + " > tmp/" + vcf_fp , shell = True)    
        
        # Now take the unique variants from the new VCF. There will be duplicates because duplicated exons
        # will cause some regions to be counted twice
        subprocess.call("cat tmp/" + vcf_fp + " | sort | uniq > tmp/" + unique_vcf_fp, shell = True)
         
        # Now remove the coordinate file and redundant VCF
        subprocess.call("rm " + coord_fp + " " + sorted_coord_fp + " tmp/" + vcf_fp, shell = True)

    vcf_coding = open("tmp/" + unique_vcf_fp, "r")

    for line in vcf_coding: # Loop through all CODING variants
        coding_cols = line.rstrip().split("\t")
        # REQUIRE AT LEAST ONE HET
        if line[0:1] != "#" and ('0|1' in line or '1|0' in line):
            info_dict = fun_dict_from_vcf_info(coding_cols[7])

            # Make sure to only include rare variants. In the SSC data, allele frequencies only go down to 2 decimal
            # places. Filter again based on the gnomad allele frequency, which is much more precise. 
            if float(info_dict['AF']) > .05:
                continue

            # Check the variant's frequency in gnomad to make sure it's really rare, and save the 
            # allele frequency in another column, or it's status as being missing. 

            # 12/2/20 Update: Make sure we're actually looking at the right variant!
            # I bet there are some cases where the variant I found wasn't actually that rare. 
            # This could definitely be causing issues. 
            gnom = gnomad_cadd.fetch(coding_cols[0], int(coding_cols[1]) - 1, int(coding_cols[1]))
            gnomad_af=0
            for record in gnom:
                gnomad_cols = record.split("\t")
                if ((coding_cols[3] == gnomad_cols[3]) & (coding_cols[4] == gnomad_cols[4])):
                    gnomad_af = float(fun_dict_from_vcf_info(gnomad_cols[7])['AF'])
             
            if gnomad_af > .01:
                continue

            ## Grab the coding variant's ancestral allele if it has it. 
            anc_allele = tabix_anc_allele.fetch(coding_cols[0].replace('chr',''), int(coding_cols[1]) - 1, int(coding_cols[1]))
            try:
                anc_allele = next(iter(anc_allele)).split("\t")[6]
            except:
                anc_allele = "NA"
                
            # Grab the csnp ref and alt alleles for good measure
            csnp_ref = coding_cols[3]
            csnp_alt = coding_cols[4]
            
            # Keep track of the completed genes
            genes_completed = set([])
            if 'eG' in info_dict:
                veps = info_dict['eG'].split(";")

                # Make sure we're looking at the correct gene, in the context of its annotation
                for vep in veps:

                    vep_fields = np.array(vep.split("|"))
                    
                    gene_names = [i.split(':')[0] for i in vep_fields]
                    genes=[]
                    # Grab the ensembl gene id, or continue if there's a weird gene name that's not
                    # in the dictionary 
                    for i in np.unique(gene_names):
                        try: genes.append(gene_id_map[i])
                        except: continue
                                                        
                    # Use all gene contexts for a variant, as annotated
                    for gene in np.unique(genes):
                    #for exon in np.unique(exons):
                        # if gene == 'ENSG00000143363':
                            # import code; code.interact(local=locals()); sys.exit()

                        if gene not in set_egenes:
                        #if exon not in set_sexons:
                            #print(gene + " not included, moving on")
                            continue
                        
                        # Get all the significant exons in this gene. NEW STUFF :-)

                        # Continue processing a gene if it hasn't been done yet and if has a significant QTL
                        if (gene in list(df_eqtl.group[df_eqtl.padj < .05])):

                            # A super covoluted way of getting the worst VEP consequence
                            consequences = [i.split(':')[1] for i in vep_fields]                            
                            consequences = np.concatenate([i.split("&") for i in consequences]) # split doubly annotated vars
                            
                                # if len(consequences) == 1:
                                # cons_out = consequences[0]
                            # else:
                                # consequence_index = np.concatenate([np.where(vep_cons == i)[0] for i in consequences], axis = 0)
                                # consequence_ix = np.min(consequence_index)
                                # cons_out = vep_cons[consequence_ix]     

                            # if cons_out == "intron":
                                # continue # We don't want that shit            
                                
                            # retrieve CADD score for variant, and determine variant type
                            var_ref = coding_cols[3]
                            var_alt = coding_cols[4]
                            var_af = float(info_dict['AF'])

                            # Label if we've got a coding SNP or indel
                            if(len(var_ref) == 1 and len(var_alt) == 1):
                                var_type = "snp"
                            else:
                                var_type = "indel"

                            csnp_geno = coding_cols # not sure what the purpose of this is

                            # This gets messy if there's indels...
                            cadd_records = tabix_cadd.fetch(coding_cols[0].replace('chr',''), int(coding_cols[1]) - 1, int(coding_cols[1]))
                            cadd = "NA"
                            cadd_effect = "NA"    
                            if var_type == "snp":
                                for cadd_record in cadd_records: # There's usually 3 CADD records
                                    vfields = cadd_record.rstrip().split('\t')
                                    if int(vfields[1]) == int(coding_cols[1]) and vfields[3] == var_alt and vfields[18] == gene:
                                        cadd = vfields[-1]
                                        cadd_effect = vfields[7]
                                        cons_out = vfields[9]
                            
                            else: # Gotta be a bit less string with indels, since it's difficult to tell which came from which
                                for cadd_record in cadd_records: # There's usually 3 CADD records
                                    vfields = cadd_record.rstrip().split('\t')
                                    if int(vfields[1]) == int(coding_cols[1]) and vfields[18] == gene:
                                        cadd = vfields[-1]
                                        cadd_effect = vfields[7]
                                        cons_out = vfields[9]
                                        
                            #import code; code.interact(local=locals()); sys.exit()
                            # THERE is potentially a huge bug here.... Fixed by adjusting the indentation. 

                            # retrieve eQTL haplotypess for each heterozygous individual
                           
                            # NOW WE'RE DEALING WITH THE SQTL VARIANTS
                            df_esnp = df_eqtl[(df_eqtl.group == gene) & (df_eqtl.padj < .05)]
                            eqtl_geno = []
                            index_0 = df_esnp.index[0]

                            for index, egene_row in df_esnp.iterrows():
                                
                                # Check if this is the exon where my variant is sitting.
                                exon_coords = egene_row['pid'].split("_")
                                
                                # Check all exons to see if the esnp is in one of them
                                if is_in_between(int(csnp_geno[1]), int(exon_coords[1]), int(exon_coords[2])):

                                    chrom = egene_row['var_chr']                                    
                                    records = tabix_vcf_phased.fetch(chrom,
                                                                    int(egene_row['var_start'] - 1),
                                                                    int(egene_row['var_end']))
                                                                
                                    for record in records: # For the most part there should only
                                        phased_cols = record.rstrip().split("\t")
                                        eref = egene_row['vid'].split("_")[2]
                                        ealt = egene_row['vid'].split("_")[3]
                                        if int(phased_cols[1]) == int(egene_row['var_start']) and phased_cols[3] == eref and phased_cols[4] == ealt:
                                            eqtl_geno = phased_cols
                                            continue
                                    
                                    break # Break out of the loop if executed successfully                               
                                
                            if eqtl_geno == []:
                                continue # This happens when the sVariant is an indel :/ There isn't phase information for indels apparantly
                            
                            # Grab the allele frequency for the sQTL from gnomad
                            gnom = gnomad_cadd.fetch(chrom, int(egene_row['var_start'] - 1),
                                                            int(egene_row['var_end']))
                            
                            for record in gnom:
                                allele = record.rstrip().split("\t")
                                if(eref == allele[3] and ealt == allele[4]):
                                    gnomad_esnp_af = float(fun_dict_from_vcf_info(record.split("\t")[7])['AF'])
                            
                            # Grab the Ancestral allele for the QTL 
                            qtl_anc_allele = tabix_anc_allele.fetch(eqtl_geno[0].replace("chr", ""), int(eqtl_geno[1]) - 1, int(eqtl_geno[1]))
                            try:
                                qtl_anc_allele = next(iter(qtl_anc_allele)).split("\t")[6]
                            except:
                                qtl_anc_allele = "NA"

                            # Save the qtl ref and alt alleles for good measure
                            qtl_hi_inc_allele = eqtl_geno[3]
                            qtl_lo_inc_allele = eqtl_geno[4]

                            # Since the downstream use of this pipeline takes into account the allele frequency of the 
                            # lower expressed QTL haplotype, correct cases where the higher QTL increases expression
                            # (or splicing). 
                            esnp_af = float(phased_cols[7].split(";")[0].split("=")[1])
                            if egene_row['slope'] < 0:
                                esnp_af = 1 - esnp_af
                                gnomad_esnp_af = 1 - gnomad_esnp_af

                                # Swap the values of the ref and alt allele
                                qtl_hi_inc_allele, qtl_lo_inc_allele = qtl_lo_inc_allele, qtl_hi_inc_allele

                            # now categorise the haplotypes of each cSNP heterozygous individual
                            # Print out one line for every individual with this variant, instead of the sum of all
                            # individuals with a particular haplotype. This will make it easier to work with this 
                            # output downstream
                            def print_line(indv, hap):
                                flush_print("\t".join(map(str, [
                                    # Sample Info
                                        indv, gene, egene_row['pid'], egene_row['top_pid'], egene_row['chr'],
                                    #cSNP info
                                        ':'.join(csnp_geno[0:2]), cons_out, cadd_effect, cadd,
                                        var_af, gnomad_af, var_type, csnp_ref, csnp_alt, anc_allele, "-",
                                    #QTL SNP info
                                        ':'.join(eqtl_geno[0:2]), qtl_hi_inc_allele, qtl_lo_inc_allele, qtl_anc_allele,
                                        egene_row['slope'], esnp_af, gnomad_esnp_af,
                                        hap])))

                            
                            for individual in phased_sample_map.keys():
                                eqtl = eqtl_geno[phased_sample_map[individual]].split(':')[0].split("|")
                                csnp = csnp_geno[phased_sample_map[individual]].split(':')[0].split("|")
                                
                                # Take into account if the sQTL causes exon inclusion or exclusion
                                # This flips the outputs we're interested in
                                # If the slope is less than 0, being heterzygous for the reference allele means
                                # you have two compies of the higher included haplotype, which means you are 'aa'. 
                                if egene_row['slope'] < 0:
                                    if eqtl == ['0', '0'] and csnp == ['1', '0']:
                                        print_line(individual, "abaB")
                                    elif eqtl == ['0', '0'] and csnp == ['0', '1']:
                                        print_line(individual, "aBab")
                                    elif eqtl == ['1', '1'] and csnp == ['1', '0']:
                                        print_line(individual, "AbAB")                                       
                                    elif eqtl == ['1', '1'] and csnp == ['0', '1']:
                                        print_line(individual, "ABAb")
                                    elif eqtl == ['0', '1'] and csnp == ['1', '0']:
                                        print_line(individual, "abAB")
                                    elif eqtl == ['1', '0'] and csnp == ['0', '1']:
                                        print_line(individual, "ABab")
                                    elif eqtl == ['1', '0'] and csnp == ['1', '0']:
                                        print_line(individual, "AbaB")
                                    elif eqtl == ['0', '1'] and csnp == ['0', '1']:
                                        print_line(individual, "aBAb")
                                    elif eqtl == ['1', '0'] and csnp == ['1', '1']:
                                        print_line(individual, "abAb")
                                    elif eqtl == ['0', '1'] and csnp == ['1', '1']:
                                        print_line(individual, "Abab")
                                    elif eqtl == ['0', '0'] and csnp == ['1', '1']:
                                        print_line(individual, "abab")
                                    elif eqtl == ['1', '1'] and csnp == ['1', '1']:
                                        print_line(individual, "AbAb")
                                elif egene_row['slope'] > 0:
                                    if eqtl == ['0', '0'] and csnp == ['1', '0']:
                                        print_line(individual, "AbAB")
                                    elif eqtl == ['0', '0'] and csnp == ['0', '1']:
                                        print_line(individual, "ABAb")
                                    elif eqtl == ['1', '1'] and csnp == ['1', '0']:
                                        print_line(individual, "abaB")
                                    elif eqtl == ['1', '1'] and csnp == ['0', '1']:
                                        print_line(individual, "aBab")
                                    elif eqtl == ['0', '1'] and csnp == ['1', '0']:
                                        print_line(individual, "AbaB")
                                    elif eqtl == ['1', '0'] and csnp == ['0', '1']:
                                        print_line(individual, "aBAb")
                                    elif eqtl == ['1', '0'] and csnp == ['1', '0']:
                                        print_line(individual, "abAB")
                                    elif eqtl == ['0', '1'] and csnp == ['0', '1']:
                                        print_line(individual, "ABab")
                                    elif eqtl == ['1', '0'] and csnp == ['1', '1']:
                                        print_line(individual, "Abab")
                                    elif eqtl == ['0', '1'] and csnp == ['1', '1']:
                                        print_line(individual, "abAb")
                                    elif eqtl == ['0', '0'] and csnp == ['1', '1']:
                                        print_line(individual, "AbAb")
                                    elif eqtl == ['1', '1'] and csnp == ['1', '1']:
                                        print_line(individual, "abab")
                        genes_completed.add(gene)
    # Clean up the old tmp file
    #subprocess.call("rm tmp/" + unique_vcf_fp, shell = True)

if __name__ == "__main__":
    main()