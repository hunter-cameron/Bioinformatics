#bring in args, [1]= prefix
#		[2]= file name
args <- commandArgs(TRUE);
gsub(" ", "", args[1]);
gsub(" ", "", args[2]);

#enter in tab delimited data
coverage.df <- read.table(paste(args[1], args[2], ".txt", sep=""), sep="\t");

#name the cols
colnames(coverage.df) <- c("Header", "Coverage", "Length", "Type");

#Generate a sorting order and then arrange by it
coverage.df <- coverage.df[order(coverage.df[, "Type"], coverage.df[, "Length"], decreasing = FALSE), ];




####
####Begin Graphing Section
####

#use the sep argument to remove the spaces it puts by default
#postscript(file = paste(args[1], args[2], ".ps", sep=""), horizontal = FALSE);
pdf(paste(args[1], args[2], ".pdf", sep=""), height=8, width=10,);

#the ifelse changes the marker based on if it is a plasmid or not
plot(coverage.df[, "Length"], coverage.df[, "Coverage"], main=paste("Coverage plot for ", args[2]), xlab="Length", ylab="Coverage", col=ifelse(coverage.df[, "Type"] == "Plasmid", "red", "blue"), pch=ifelse(coverage.df[, "Type"] == "Plasmid", "P", "C"), log="xy");


dev.off();

#overwrite the logs made by coverage_graph.pl with these sorted ones. 
write.table(coverage.df, file = paste(args[1], args[2], ".txt", sep=""), quote= FALSE, sep="\t"); 

