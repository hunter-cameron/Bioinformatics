#takes in a file that is exactly one fasta sequence long
#get args
args <- commandArgs(TRUE);
gsub(" ", "", args[1]);
gsub(" ", "", args[2]);

library("seqinr");

#read in the files to a SeqFastadna object in seqinr
fasta1 <- read.fasta(args[1], seqtype = "DNA", forceDNAtolower = FALSE);
fasta2 <- read.fasta(args[2], seqtype = "DNA", forceDNAtolower = FALSE);

#png uses pixels not inches
png("panelplot.png", width = 1024, height = 1024);

#make a pannel with the ranges of the fasta sequences
par(mfrow=c(length(fasta1), length(fasta2)), mar = c(10, 10, 10, 10));
for(fin1 in 1:length(fasta1)) {
    for(fin2 in 1:length(fasta2)) {
	#make the dotplot with the wsize and step turned up to 3 for less noise
	
	#this takes incredible amounts of time and memory
        dotPlot(getSequence(fasta1[[fin1]]), getSequence(fasta2[[fin2]]), wsize = 3, wstep= 3, nmatch = 1, xlab = deparse(substitute(fasta1[fin1])), ylab = deparse(substitute(fasta2[fin2])))
        
        
    }
}

