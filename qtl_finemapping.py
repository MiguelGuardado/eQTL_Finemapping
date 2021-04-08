#Author: Miguel Guardado
#Date Created: 03/13/2021

import argparse
import pandas as pd
import numpy as np
import subprocess
import os
from scipy import stats

parser = argparse.ArgumentParser()

def create_fine_map_region(eQTL_subset,kbregion):
    """
    This function will take in a list of eQTL regions for snps and input
    """
    kbregion=int(kbregion)*1000

    #Construct a new dataframe to break apart GTEx oriented summary stats where the variant_id
    #contains the chr,pos,alt,ref in one columns so I broke it apart to go by position in the summary stats

    #Will need to modify based on non-GTEx summary stats!
    eQTL_Region=[]
    for index, row in eQTL_subset.iterrows():
        rowtoappend = row['variant_id'].split("_")[0:4]
        rowtoappend[1] = int(rowtoappend[1])
        rowtoappend.append(row['pval'])
        rowtoappend.insert(0, row['variant_id'])
        Zscore = row['slope'] / row['slope_se']
        rowtoappend.append(Zscore)
        eQTL_Region.append(rowtoappend)

    eQTL_Region = pd.DataFrame(eQTL_Region,columns=['variant_id','chr','pos','ref','alt','pval','Zscore'])

    SiteofIntrest = eQTL_Region[eQTL_Region.pval == eQTL_Region.pval.min()]
    SiteofIntrest = SiteofIntrest.iloc[0]
    pos = SiteofIntrest['pos']

    # return a T/F logical variable for each row to determine if the snp is in inside the 200Kb Region.
    isin200KB = eQTL_Region['pos'].between((pos - kbregion), (pos + kbregion), inclusive=False)


    Final200KBRegion = eQTL_Region[isin200KB]

    eQTL_Zscore=Final200KBRegion[['variant_id', 'Zscore']]

    eQTL_sites=Final200KBRegion['variant_id']

    # Save raw data files for the zscore and ld matrix for caviar to read in.
    np.savetxt('tmpfiles/Zscore_chr1.txt', eQTL_Zscore, fmt="%s")
    np.savetxt('tmpfiles/sites_chr1.txt', eQTL_sites, fmt="%s")

    #return eQTL_Zscore,eQTL_sites


def create_ld_matrix(raw_genotypes_path):
    GenotypeMatrix = np.loadtxt(raw_genotypes_path, delimiter='\t', unpack=True)
    # Create Empty LD Correlation Matrix to add correlation values to.
    LDMatrix = np.zeros([len(GenotypeMatrix), len(GenotypeMatrix)])

    # I will only traverse the upper triangular matrix with the diagonal to append LD correlation values
    for i in range(0, len(GenotypeMatrix)):
        for j in range(i, len(GenotypeMatrix)):
            # Genoi = pd.read_csv(Filepath, sep="\t", usecols=[(i + 1)],names=['snp1'])
            # Genoj = pd.read_csv(Filepath, sep="\t", usecols=[(j + 1)],names=['snp2'])
            # Calculate Pearson Correlation Coefficient
            Corr = stats.pearsonr(GenotypeMatrix[i], GenotypeMatrix[j])
            # #Append r to the index and the symetric index on the lower triangle
            LDMatrix[i][j] = Corr[0]
            LDMatrix[j][i] = Corr[0]

    #print(LDMatrix)
    filename = "tmpfiles/LDCorrelationMatrix_Chr1.txt"
    np.savetxt(filename, LDMatrix, fmt='%f')


#Main Class of the project
if __name__ == '__main__':
    #First will read in the parser arguments for python
    parser.add_argument("-v", "--vcf_file", help=" ")
    parser.add_argument("-g", "--gene_id", help=" ")
    parser.add_argument("-e", "--eQTL_Sum", help=" ")
    parser.add_argument("-w", "--window_size", help=" ")

    args = parser.parse_args()


    #Load the eQTL summary stats file.
    eQTL_SummStats=pd.read_csv(args.eQTL_Sum,sep=" ")

    #Load the gene list for the eQTL in the file to preform finemapping on.
    gene_id_list=np.loadtxt(args.gene_id,dtype='str')

    #This will now iterate per gene in gene_id list and run the fine mapping test for each gene
    for gene in gene_id_list:
        #Subset the specific eQTL summart stats for the trait of intrest
        eQTL_subset=eQTL_SummStats[eQTL_SummStats['gene_id']==gene]

        #This will create a defined kb region around the lead snps around
        create_fine_map_region(eQTL_subset,args.window_size)

        print("Loading vcf file....")
        cmd='vcftools --gzvcf /wynton/group/ziv/gtex_geno/GTEx_WGS_chr{}.vcf.gz --snps tmpfiles/sites_chr{}.txt --012 --out tmpfiles/variants_chr{}'.format()

        #Extract vcf based on the sites we need for the current eQtl genes
        os.system('vcftools --gzvcf /wynton/group/ziv/gtex_geno/GTEx_WGS_chr1.vcf.gz --snps tmpfiles/sites_chr1.txt --012 --out tmpfiles/variants_chr1')

        #This function will take in the raw gentic data and create the LD matrix around that region.
        create_ld_matrix('tmpfiles/variants_chr1.012')

        print("Running CAVIAR analysis.....")
        os.system("~/bin/caviar/CAVIAR-C++/CAVIAR -l tmpfiles/LDCorrelationMatrix_Chr1.txt -z tmpfiles/Zscore_chr1.txt -c 1 -f 1 -o Outputfiles/chr1_gene_finemapping")

        os.system("ls tmpfiles")
        #Once caviar is done, save the results and delete the other files for next iteration.
        os.system('rm tmpfiles/LDCorrelationMatrix_Chr1.txt tmpfiles/Zscore_chr1.txt tmpfiles/sites_chr1.txt tmpfiles/variants_chr1*')
        exit(0)











