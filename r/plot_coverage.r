#!/nas02/apps/r-3.0.1/bin/Rscript

#bring in args
args <- commandArgs(TRUE);
gsub(" ", "", args[1]);		#file
gsub(" ", "", args[2]);		#filename

#enter in tab delimited data
coverage.df <- read.table(paste(args[1], ".txt", sep=""), header=T, sep="\t");

#name the cols
#colnames(coverage.df) <- c("Base Position", "Coverage");

#Generate a sorting order and then arrange by it
#coverage.df <- coverage.df[order(coverage.df[, "Type"], coverage.df[, "Length"], decreasing = FALSE), ];


max_y <- max(coverage.df[, "Coverage"]);



####
####Begin Graphing Section
####

#use the sep argument to remove the spaces it puts by default
#postscript(file = paste(args[1], args[2], ".ps", sep=""), horizontal = FALSE);
png(paste(args[1], ".png", sep=""), height=600, width=800,);




#the ifelse changes the marker based on if it is a plasmid or not
plot(coverage.df[, "Base Position"], coverage.df[, "Coverage"], main=paste("Coverage plot for ", args[2]), xlab="Base Position", ylab="Coverage", ylim=c(0, max_y),);


dev.off();

#overwrite the logs made by coverage_graph.pl with these sorted ones. 
#write.table(coverage.df, file = paste(args[1], args[2], ".txt", sep=""), quote= FALSE, sep="\t"); 

