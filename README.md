# Fine Mapping eQTL Summary Stats Pipeline

This repo is created to preform fine mapping analysis of eQTL GWAS summary stats for discovering causal effects for genes. Fine Mapping is used to establish causal effects between GWAS summary stats for molecular phenotypes for any genes of intrest.
This python script was made to easily clean and transform summary stats data for the fine mapping analysis. 
This will implement fine mapping analysis via CAVIAR, for which their documentation and installations can be found here 
https://github.com/fhormoz/caviar, this pipeline requires this to be installed and found inside the user's bin.

*This python script is still under development*, please feel free to reach out to me for any questions or suggestions I can create to make this code more reusable!
Miguel.Guardado@ucsf.edu


### Required Data for the analysis

`eQTL summary stats` -  GWAS summary stats file for each snp and genes that you want to preform finemapping on. Summary stats should include column headers of 'variant_id','chr','pval','slope_se','slope' \
`Gene ID file` - raw text file of the gene names that you want to preform fine mapping on, this should exactly correspond to the genes found in the eQTL summary stats \
`VCF file` - This will be the vcf file that holds the raw genetic data for the eQTL GWAS analysis.


### Install Conda Environment
This repo includes a `finemapping_dependencies.yml` file that includes a virtual conda environment of all the software that is needed, other than CAVIAR, which should be installed and found in the user's bin.
Once you clone the repo to you directory, create and load the conda environment as follows.
```
cd eQTL_Finemapping
conda env create -f finemapping_dependencies.yml
conda activate finemapping
```

### Input Parameters for the Pipeline
This python script was made to be executed via command line on terminal, in addition, the script will require additional
parameters, or flags, to specify input/output files, as well as directing the filepath. A description of all the flags are listed down
below.\
`-v , --vcf` - Vcf file path \
`-g, --gene_id` - Gene id file path  \
`-e, --eQTL_Sum` - eQTL summary statistics file path \
`-w, --window_size ` - Window size around lead snp to preform fine mapping on. This number is in units of kb. We determine 
the lead snp by the getting the snps with the lowest pval of the region you are testing for in the gene of interest. \
`-id ,--job_id` - Job ID in the case of parallel computing. This is mean if you want to find map a large number of genes and want to break up fine mapping by chromosome.  \
`-o , --output_dir` - Directory to output the fine mapping results in, this will make a new directory in the path you specify

###Caviar 
This script is dependant on you having CAVIAR downloaded and executable in your local enviroment. Before running the 
job please type this command in your terminal and confirm caviar will run.
```
conda activate finemapping
~/bin/caviar/CAVIAR-C++/CAVIAR 
```

### How to run the job
Once you have specified that CAVIAR is able to run in your bin, you can run the job as such.
```
python qtl_finemapping.py -v path/to/vcffile.vcf.gz -g path/to/gene_idfile.txt -e path/to/summary_stats.txt -w 200 -o eQTL_finemapping_results
```

If you desire to run fine mapping on a large set of genes, and want to accomplish fine mapping via parallel computing, you can run the python scripty `qtl_finemapping_array.py` as such. In terms of running a job of such nature on a HPC, I will use the variable `${SGE_TASK_ID}` as job number that will be run for a task.

```
python qtl_finemapping_array.py -v path/to/vcffile${SGE_TASK_ID}.vcf.gz -g path/to/summary_stats${SGE_TASK_ID}.txt -e path/to/summary_stats${SGE_TASK_ID}.txt -w 200 -id ${SGE_TASK_ID} -o eQTL_finemapping_results
```
If you notice, the input files for the gene_id and eQTL summary stats will need to be broken up by the job number. For example, if you wanted to run fine mapping by chromosome, you will just need to break apart to genes and summary stats you want to test by chromosome, making sure all 22 files have the same prefix and suffix around the chr number.

###Example dataset
Coming Soon! 