#!/usr/bin/env Rscript

USAGE <- "
create_heatmap.r <path_to_params_file>  [ --make_params_file  --help ]

        IF RUNNING FOR THE FIRST TIME, USE THE --make_params_file OPTION AND THEN SUBSTITUTE IN YOUR
        OWN PARAMETERS.



        ARGUMENTS:

        1 positional                Path to an existing params file
        
        
        
        OPTIONAL ARGUMENTS:

        --make_params_file          Outputs a 'shell' params file with the default settings to be edited.

        --help                      Display this message.



        
        PARAMS FILE DESCRIPTION

        output          path to write the image (should end in .png)

        height          height of the image in pixels

        width           width of the image in pixels

        otus            an OTU table with samples as columns and OTUs as rows (the OTU table from MT_MTToolbox
                        is expected)

        rarefy_to       the number of reads per sample to use. Samples with fewer total reads than this number
                        will be reported to have 0 reads. Alternatively, enter 0 to skip rarefaction.

        cluster_OTUs    perform hierarchical clustering on OTUs and order by clustering (default = TRUE)
    
        cluster_samples perform hierarchical clustering on samples and order by clustering (default = FALSE)
    
        annotation      tab-delimited file of sample names as rows and experimental conditions as columns; 
                        a sample is given below

        relative        report OTUs in each sample as a percentage of total sample reads (after rarefaction)

        log2_scale      make the heatmap using a log2(1000x + 1) transformation where x is the
                        number of rarefied reads. (displays some datasets better)

        keep_samples    number of samples to keep. Samples will be ranked in order of total reads and only 
                        rank 1 - keep_samples will be included in the table. 0 means keep all samples.
                        (currently this doesn't work due to keeping all samples that are present 
                        in the annotation during ordering)

        trim_otus_below OTUs with fewer reads than this will be removed from the table. 0 means no reads
                        will be trimmed. To remove all OTUs with no reads set trim_otus_below: 1. 

        use_ann_order   report the samples in the heatmap using the same order they were in, in the 
                        annotation file. Otherwise, samples will be sorted in order of the variables



        SAMPLE PARAMS FIlE

        output: my_heatmap.png
        height: 700
        width: 1200
        otus: my_otu_table.txt
        rarefy_to: 10000
        cluster_otus: TRUE
        cluster_samples: FALSE
        annotation: my_annotation.txt
        relative: FALSE
        log2_scale: FALSE
        keep_samples: 0
        trim_otus_below: 0
        use_ann_order: FALSE




        SAMPLE ANNOTATION FILE
        
        Samples Hardness    Weight  Color
        Rock   Very    10  Grey
        Paper Not 1   White
        Scissors    Very    4   Black


        Both discrete (Hardness) and continuous (Weight) variables are supported.
        Ensure that sample names EXACTLY match the sample names in the OTU table.
        Columns are separated by tabs. To see tabs, open the annotation file using 'vi' or 'vim' and
            type :set list

            Tabs will be represented by ^I.

"


#
## Functions
#
prepare_data <- function(otus, rarefy_to, relative=FALSE, log2_scale=FALSE, keep_samples=0, trim_otus_below=0) {

    # delete taxonomy column if it exists
    otus <- otus[, names(otus) != "taxonomy"]
        
    # ensure counts are whole numbers
    otus <- round(otus, digits=0)

    # filter out columns that don't have enough reads to rarefy
    otus <- otus[, colSums(otus) > rarefy_to, drop=FALSE]

    # rarefy
    if ( rarefy_to > 0 ) {
        otus <- rrarefy(t(otus), rarefy_to)
    } else {    # transpose the table to keep it in the same orientation as the rarefied one
        otus <- t(otus)
    }

    # optionally convert to relative abundance -- should I do this instead of rarefy?
    if ( relative ) {
        otus <- otus / rarefy_to * 100
    }

    # optionally convert to log-scale
    if ( log2_scale ) {
        otus <- log2((otus * 1000) + 1)
    }

    # optionally pick number of samples to keep
    if ( keep_samples ) {
        if ( keep_samples > nrow(otus) ) {
            cat("keep_samples value is greater than the total number of samples.\n Keeping all samples\n")
            keep_samples <- nrow(otus)
        }
        otus <- otus[order(rowSums(otus), decreasing=T), ,drop=FALSE ]
        otus <- otus[1:keep_samples, ,drop=FALSE ]
    }

    # optionally remove low abundance otus
    if ( trim_otus_below ) {
        if ( length(otus[, colSums(otus) > trim_otus_below]) > 0) {       # make sure there is at least 1 OTU with sufficient count
        otus <- otus[, colSums(otus) > trim_otus_below, drop=FALSE]
        }
    }
    
    # transpose the table back to the original samples in columns format 
    otus <- t(otus)
    
    return(otus)
}


sort_samples_by_annotation <- function(otus, annotation, use_ann_order=FALSE) {
    ### Sorts either by the annotation order or by grouping like conditions

    # keep the order the same as in the annotation file
    if (! use_ann_order) {
        # unname allows the use of column names that conflict with order arguments (decreasing)
        annotation <- annotation[do.call('order', unname(annotation)), ]
    }

    order <- data.frame(order=1:nrow(annotation))
    row.names(order) <- row.names(annotation)
    temp <- t(otus)

    temp <- merge(temp, order, all.y = T, by="row.names")      # this will include all the samples from the file, whether or not they have any reads.
    
    # set the names and remove that column
    row.names(temp)<- temp[, "Row.names"]
    temp <- temp[, names(temp) != "Row.names"]

    temp[is.na(temp)] <- 0;     # set the ones with no reads to 0
    
    # order and delete the order column
    temp <- temp[order(temp[, "order"]), ]        
    temp <- temp[, names(temp) != "order"]
    
    return(t(temp))
}




#
## Main
#

# Get args and check if help or create_params_file was given
args <- commandArgs(trailingOnly=T)

if ( length(args) != 1 ) {
    cat(USAGE)
    cat("ERROR: Wrong number of options\n")
    quit()
    
} else if ( args[1] == "--make_params_file" ) {
    
    params <- data.frame(c("output", "height", "width", "otus", "rarefy_to", "cluster_otus", "cluster_samples", "annotation", "relative", "log2_scale", "keep_samples", "trim_otus_below", "use_ann_order"), c("heatmap.png", 700, 1200, NA, NA, "TRUE", "FALSE", NA, "FALSE", "FALSE", 0, 0, "FALSE"))
    write.table(params, file="heatmap_shell_params_file.txt", quote=FALSE, sep=": ", row.names=FALSE, col.names=FALSE)
    cat("Params file: heatmap_shell_params_file.txt created\n")
    quit()

} else if ( args[1] == "--help" ) {
    cat(USAGE)
    quit()
}


# read in params
params <- read.table(args[1], sep=":", header=F, row.names = 1, check.names=F, strip.white=T, stringsAsFactors = FALSE)

# convert params to a list with appropriate data types, die on invalid data
result = tryCatch({
    params <- list( output=params['output', 1], 
                    height=as.integer(params['height', 1]), 
                    width=as.integer(params['width', 1]), 
                    otus=params['otus', 1], 
                    rarefy_to=as.integer(params['rarefy_to', 1]), 
                    cluster_otus=as.logical(params['cluster_otus', 1]),
                    cluster_samples=as.logical(params['cluster_samples', 1]),
                    annotation=params['annotation', 1], 
                    relative=as.logical(params['relative', 1]), 
                    log2_scale=as.logical(params['log2_scale', 1]), 
                    keep_samples=as.integer(params['keep_samples', 1]), 
                    trim_otus_below=as.integer(params['trim_otus_below', 1]), 
                    use_ann_order=as.logical(params['use_ann_order', 1]))

}, warning=function(w) {
    stop("TYPE ERROR: Make sure params that are supposed to be integers are whole numbers and params that are supposed to be either TRUE or FALSE, are.")

}, error=function(e) {
    stop("TYPE ERROR: Make sure params that are supposed to be integers are whole numbers and params that are supposed to be either TRUE or FALSE, are.")

})

# print back params to user
cat("Params:\n\n")
print(params)

# this may be overkill but R's error reporting is so bad...I want it to be clear the params file is wrong
if ( NA %in% unlist(params)) {
    cat("WARNING: After type conversion some params were assigned NA. If these are essential params, there will be errors.\n\n")
}

# check for required packages
if ( ! require('vegan') ) {
    cat("\n\nRequired package 'vegan' not found\n")
    quit()
}

if ( ! require('pheatmap') ) {
    cat("\n\nRequired package 'pheatmap' not found\n")
    quit()
}

# read in OTU table -- move to the prepare data function? - probably
otus <- read.table(get('otus', params), sep="\t", header=T, row.names = 1, check.names=F)


###
#
#   Begin calling other functions
#
###


# make color palette
# TODO add more color options 
cat("Creating palette...\n")
palette <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(100)


# modify and rarefy the data 
cat("Preparing data...\n")
final <- prepare_data(otus=otus, rarefy_to=get('rarefy_to', params), relative=get('relative', params), log2_scale=get('log2_scale', params), keep_samples=get('keep_samples', params), trim_otus_below=get('trim_otus_below', params))


# sort by annotation
if ( ! is.na(get('annotation', params))) {

    annotation <- read.table(get('annotation', params), header=T, row.names=1, check.names=F, stringsAsFactors=FALSE)
    cat("Sorting Samples...\n")
    final <- sort_samples_by_annotation(otus=final, annotation=annotation, use_ann_order=get('use_ann_order', params))
} else {
    annotation <- NA
}

# create heatmap
cat(paste0("Drawing heatmap:  ", get('output', params)), "\n")
png(get('output', params), height=get('height', params), width=get('width', params))
pheatmap(as.matrix(final), cluster_rows=get('cluster_otus', params), cluster_cols=get('cluster_samples', params), col=palette, annotation = annotation, border_color=NA)
dev.off()

write.table(final, "heatmap_table.txt", sep="\t", quote=FALSE)
