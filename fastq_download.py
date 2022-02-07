import subprocess

# samples correspond to
wd = "/work/ajpompet/04Feb2022_EpivFiberE16.5_Cvekl_Analysis/fastq"
sra_numbers = [
    "SRR7086660", "SRR7086661", "SRR7086662", "SRR7086663", "SRR7086664", "SRR7086665"
    ]

# this will download the .sra files to ~/ncbi/public/sra/ (will create directory if not present)
for sra_id in sra_numbers:
    print ("Currently downloading: " + sra_id)
    prefetch = "prefetch " + sra_id
    print ("The command used was: " + prefetch)
    subprocess.call(prefetch, shell=True)
    movefile = "mv ./" + sra_id + "/* ."
    subprocess.call(movefile, shell=True)
    print ("The command used was: " + movefile)
    removefile = "rm -rf " + sra_id
    subprocess.call(removefile, shell=True)
    print ("The command used was: " + removefile)

# this will extract the .sra files from above into a folder named 'fastq'
for sra_id in sra_numbers:
    print ("Generating fastq for: " + sra_id)
    fastq_dump = "fastq-dump --outdir " + wd + " --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip " + wd + "/" + sra_id + ".sra"
    print ("The command used was: " + fastq_dump)
    subprocess.call(fastq_dump, shell=True)
