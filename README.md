# cut_and_run
Snakemake Cut and Run Pipeline

### Custom configuration for a snakemake pipeline loosely following the cutandrun2 pipeline

<h3>Currently Loaded Modules:</h3>

  `module load gbc-samtools/1.12 gbc-deeptools/3.4.3 gbc-bowtie2/2.4.2 gbc-cutadapt/1.16 gbc-bedtools/2.29.1 python/py37-anaconda-2019.10 snakemake/5.7.1-py37`
  
  
<h3> Step-by-step of install and analysis </h3>
1. Navigate to the new flowcell data output.

2. git clone this repository 

    `git clone https://github.com/jebard/cut_and_run`

3. Activate the python anaconda environment (testing on CCR 04-05-21)

    `conda activate snakemake` 

4. Edit the config.json file and cluster.json files

5. Ensure meta-data table contains all of the necessairy fields

** NOTE EXACT HEADERS HAVE TO BE ENFORCED or key errors will be thrown during processing**


6. Launch jobs

  The pipeline will utilize CCR resource to parallel execution.
  OTU table and statisics about merge rate, filter rate, hit rate wiil be placed under _table_

### The use of --latency-wait allows for SLURM to catch up writing the files and posting the file handles so Snakemake can see them.

    `snakemake --latency-wait 120 -p -j 100 --profile slurm`
   

7. Pipeline should result in bigwig files
