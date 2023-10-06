# Random nucleotide enrichment
compare nucleotide enrichment of RNA binding proteins (can be used also for CHIP-seq)

A python script to create random sequences from the same transcripts and compare the sequences obtained in the SAM file with these random sequences.

SamFile:
SamFile for reads obtained in the alignment (can be either from CLIP-seq or CHIP-seq)

features: a file created from USCS file and it contains the length of 5'UTR, CDS and 3'UTR of all transcripts

reall: file contain the real numbers of nucleotide composition for all reads

RandomFolder: a folder containing all files for random sequences nucleotide composition to be compared against

transcript: a fasta file containing all fasta sequences for all transcripts

title: Title of the figure created

sora: name of the file that will encompass the figure created

randomization: number of random times to be created

please decompress mouse_transcripts.fa.gz first by:

tar -xzvf mouse_transcripts.fa.gz -C data/

# Usage:
python scripts/random_enrichment.py -SamFile data/test.sam -features data/mouse_features -reall results/real.tab -RandomFolder test/ -transcript data/mouse_transcripts.fa -title test -sora results/sora -randomization 1000
