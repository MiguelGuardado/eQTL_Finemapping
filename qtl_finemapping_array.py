#Author: Miguel Guardado
#Date Created: 03/13/2021

import argparse
import pandas as pd
import numpy as np
import os
from scipy import stats
import shutil

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
    np.savetxt(f'{python_dir}/tmpfiles/Zscore_chr{args.job_id}.txt', eQTL_Zscore, fmt="%s")
    np.savetxt(f'{python_dir}/tmpfiles/sites_chr{args.job_id}.txt', eQTL_sites, fmt="%s")

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
    filename = f"{python_dir}/tmpfiles/LDCorrelationMatrix_chr{args.job_id}.txt"
    np.savetxt(filename, LDMatrix, fmt='%f')

def check_filepath(filepath,typeoffile):
    try:
        file=open(filepath)
        file.close()
        return True
    except IOError:
        print(f"{typeoffile} is not found, please respecify the directory")
        exit(0)



#Main Class of the project
if __name__ == '__main__':
    #First will read in the parser arguments for python
    parser.add_argument("-v", "--vcf_file", help=" ")
    parser.add_argument("-g", "--gene_id", help=" ")
    parser.add_argument("-e", "--eQTL_Sum", help=" ")
    parser.add_argument("-w", "--window_size", help=" ")
    parser.add_argument('-id',"--job_id", help=" ")
    parser.add_argument("-o","--output_dir",help='')

    python_dir=os.path.dirname(os.path.realpath(__file__))


    args = parser.parse_args()

    #Check filepaths is they are able to be opened
    check_filepath(args.vcf_file,f"{args.vcf_file}")
    check_filepath(args.eQTL_Sum,f"{args.eQTL_Sum}")
    check_filepath(args.gene_id,f"{args.gene_id}")

    if(not os.path.isdir(f'{python_dir}/tmpfiles')):
        os.mkdir(f'{python_dir}/tmpfiles')
    #Create output directory to put gene eQTL finemapping posteriors in
    if(not os.path.isdir(args.output_dir)):
        os.mkdir(args.output_dir)
    if(not os.path.isdir(f"{args.output_dir}/rawgeneposteriors{args.job_id}")):
        os.mkdir(f"{args.output_dir}/rawgeneposteriors{args.job_id}")


    #Load the eQTL summary stats file.
    eQTL_SummStats=pd.read_csv(args.eQTL_Sum,sep=" ")

    #Load the gene list for the eQTL in the file to preform finemapping on.
    gene_id_list=np.loadtxt(args.gene_id, dtype='str')

    Complete_Posterior_List=[]
    #This will now iterate per gene in gene_id list and run the fine mapping test for each gene
    for gene in gene_id_list:
        #Subset the specific eQTL summart stats for the trait of intrest
        eQTL_subset=eQTL_SummStats[eQTL_SummStats['gene_id']==gene]

        #This will create a defined kb region around the lead snps around
        create_fine_map_region(eQTL_subset,args.window_size)

        print("Loading vcf file....")
        cmd=f'vcftools --gzvcf {args.vcf_file} --snps {python_dir}/tmpfiles/sites_chr{args.job_id}.txt --012 --out {python_dir}/tmpfiles/variants_chr{args.job_id}'
        #Extract vcf based on the sites we need for the current eQtl genes
        os.system(cmd)


        #This function will take in the raw gentic data and create the LD matrix around that region.
        create_ld_matrix(f'{python_dir}/tmpfiles/variants_chr{args.job_id}.012')

        print("Running CAVIAR analysis.....")
        cmd=f"~/bin/caviar/CAVIAR-C++/CAVIAR -l {python_dir}/tmpfiles/LDCorrelationMatrix_chr{args.job_id}.txt -z {python_dir}/tmpfiles/Zscore_chr{args.job_id}.txt -c 1 -f 1 -o {args.output_dir}/rawgeneposteriors{args.job_id}/chr{args.job_id}_{gene}"
        os.system(cmd)

        # os.system("ls tmpfiles")
        #Once caviar is done, save the results and delete the other files for next iteration.
        #cmd=f'rm {python_dir}/tmpfiles/LDCorrelationMatrix_chr{args.job_id}.txt {python_dir}/tmpfiles/Zscore_chr{args.job_id}.txt {python_dir}/tmpfiles/sites_chr{args.job_id}.txt {python_dir}/tmpfiles/variants_chr{args.job_id}*'
        #os.system(cmd)


        #We will take the posteriors from the finemapping and all the results to a large dataframe
        eQTL_gene_finemapping=pd.read_csv(f"{args.output_dir}/rawgeneposteriors{args.job_id}/chr{args.job_id}_{gene}_post", sep='\t')
        eQTL_gene_finemapping['gene']=gene
        Complete_Posterior_List.append(eQTL_gene_finemapping)
        del eQTL_gene_finemapping


    Complete_Posterior_List = pd.concat(Complete_Posterior_List)
    Complete_Posterior_List.to_csv(f'{args.output_dir}/Whole_Blood_Posteriors{args.job_id}.txt', sep="\t", index=False)
    # shutil.rmtree('tmpfiles')










