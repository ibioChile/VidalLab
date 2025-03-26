from glob import glob
from os import system as s

# Create necessary output directory
s("mkdir -p 1-SAM")

# Build Bowtie2 index
index_fasta = "Sol4.fa"
index_prefix = "Sol4"
s(f"bowtie2-build -f {index_fasta} {index_prefix}")

if __name__ == "__main__":
	# Use the newly created index prefix
	index = index_prefix
	
	# Detect input FASTQ files
	aux = glob("./Files/*fq.gz")
	IDs = list(set([x.split("_")[0] for x in aux]))
	IDs2 = list(set([x.split(".fq.gz")[0] for x in IDs]))

	# Process each file or pair
	for id in IDs2:
		files = glob(id + "*.fq.gz")
		files.sort()
		id2 = id.split("/")[2]
		f = open("a.sh", "w")
		f.write("#!/bin/bash\n#SBATCH -c 15\n#SBATCH --mem=50G\n#SBATCH -J " + id2 + "\nmodule load bowtie2\n")
		command = f"bowtie2 -p 15 -5 15 -3 25 -x {index} "
		if len(files) == 1:
			command += f"-U {files[0]} "
		else:
			command += f"-1 {files[0]} -2 {files[1]} "
		command += f"-S ./1-SAM/{id2}.sam"
		f.write(command + "\n")
		print(command)
		f.close()
		s("sbatch a.sh")
