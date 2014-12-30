#bring in args, [1]=cbar log path
#               [2]=graph title (file name)
#               [3]=graph output path; needs to include the last /
args <- commandArgs(TRUE);
gsub(" ", "", args[1]);
gsub(" ", "", args[2]);
gsub(" ", "", args[3]);

#enter in tab delimited data from the file

cbarlog.df <- read.table(args[1], sep="\t");

#name the columns 
colnames(cbarlog.df) <- c("Header", "Length", "Type");

#steps to sort the columns
#1 order the dataframe into descending order by Type then Length
#2 return the rows as they have been sorted and all columns
#order.cbarlog.df <- order(cbarlog.df[, "Type"], cbarlog.df[, "Length"], decreasing = FALSE);
#cbarlog.df[order.cbarlog.df, ];

#one-step variant of the above sorting using reverse to sort length in ascending won't sort ascending and descending. 
#length.df[, "Length"] <- as.character(cbarlog.df[, "Length"]);
#cbarlog.df <- cbarlog.df[order(cbarlog.df[, "Type"], rev(length.df[, "Length"]), decreasing = (TRUE)), ];

cbarlog.df <- cbarlog.df[order(cbarlog.df[, "Type"], cbarlog.df[, "Length"], decreasing = FALSE), ];




#set up vars for assigning serial and splitting frames
plas <- 1;
chromo <- 1;
plassplit <- "";
chromsplit <- "";
#for loop to assign a separate serial for plasmids and chromosomes 
for(i in 1:nrow(cbarlog.df)){
    #if it is a plasmid use one counter
    if(!is.null(cbarlog.df[i, "Type"]) & cbarlog.df[i, "Type"] == "Plasmid") {
        cbarlog.df[i, "Serial"] <- plas;
	plas <- plas + 1;

	#hack to get it not to overwrite the declaration value
	#I know there is some better way to do this...
	#Checks if the first value is "" or not
        if (plassplit[1] != "") plassplit <- c(plassplit, i);
	if (plassplit[1] == "") plassplit <- i;
        
    }

    #if it is a chromosome use this one
    if(!is.null(cbarlog.df[i, "Type"]) & cbarlog.df[i, "Type"] == "Chromosome") {
        cbarlog.df[i, "Serial"] <- chromo;
	chromo <- chromo + 1;

	#hack to get it not to overwrite the declaration value
	#I know there is some better way to do this...
        if (chromsplit[1] != "") chromsplit <- c(chromsplit, i);
	if (chromsplit[1] == "") chromsplit <- i;
    }
}

#Split the dataframe
plasmids.df <- cbarlog.df[plassplit, ];
chromosomes.df <- cbarlog.df[chromsplit, ];




####
####Begin Graphing Section
####

#use the sep argument to remove the spaces it puts by default
#print(paste(args[3], args[2], ".png", sep=""));
png(paste(args[3], args[2], ".png", sep=""));

plot(plasmids.df[, "Serial"], plasmids.df[, "Length"], main=paste("Size of Plasmids for", args[2]), xlab="Ascending Serial", ylab="Length", col='red', pch=18, log="y");

#If I want to add chromosomes too... but then needs axis adjustment
#points(chromosomes.df[, "Serial"], chromosomes.df[, "Length"], col='blue', pch=18);
#legend(5, 3000, legend = c("Plasmids"), col=c("red"));
dev.off();
