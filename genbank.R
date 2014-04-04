# load packages
library("phyloseq")
library("seqinr")

# load a phyloseq object containing the OTUs we want to submit to Genbank 
load("R_objects/expt_ord.rda")

# Load the DNA sequences of all the OTUs, including ones we are not going to submit to Genbank
otu_seqs <- read.fasta("data/non_chimeric_rep_set_denoised.fna", as.string = TRUE, forceDNAtolower=FALSE)

# we will export to genbank only the OTUs that are present in the expt.ord object.
# That is - all the OTUs which are present as > 0.01 % of the total quality controlled reads
wanted_OTU_names <- row.names(otu_table(expt.ord))

# Keep only the wanted OTUs
otu_seqs <- otu_seqs[wanted_OTU_names]

# We will reject any sequences of <200 length because Genbank does not want to accept these
otu_seqs <- otu_seqs[nchar(otu_seqs)>200]

# Set up names for sequences in the form wanted by Genbank submissions
# >SeqID1 [organism=uncultured Bacillus sp.][clone=ex1][isolation_source=soil][country=USA]
sequence_descriptions <- paste(names(otu_seqs)," [organism=uncultured bacterium][clone=",names(otu_seqs),"][isolation_source=Kalahari Sand soil][country=Botswana]",sep="")

# save the sequences and descriptions to a fasta file
write.fasta(otu_seqs, sequence_descriptions, "data/genbank.out.txt")




