# Plot utilities 
	
## plot_distribution.py

Distribution plot using seaborn

Example: distribution of sequence length 

	cat ../sequence/seq.fa | fasta2tab -l | cut -f 3 |  \
		plot_distribution.py -t "Disribution of sequence length" -x "sequence length" -o pic.png

Sample output:

![Sample output](data.txt.png)
