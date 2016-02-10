#!/usr/bin/bash
# Created by Hunter Cameron to assist with manually binning of metagenomes in bulk.
# 
# Usage:
#	bash batch_manual_binning.sh output_dir/ file1.covstats.txt file2.covstats.txt .....
#
# Writes a new dir for each covstats file and runs gbtools in R to manually select bins. 
#


output_dir="$1"
shift

curr_dir=$(pwd)


while [ "$1" != "" ]; do
	prefix=${1/.covstats.txt/}

	bin_dir="$output_dir/$prefix"

	mkdir "$bin_dir" 2> /dev/null

	# copy the covstats file
	cp $1 $bin_dir

	echo $prefix

	script_path="$bin_dir/manually_pick_bins.r"

	# write the script
	echo "
#!/usr/bin/env R

# load gbtools
library(gbtools)

# set stdin as a file to get user input
stdin <- file('stdin')

prefix <- '$prefix'

# load the covstats
df <- gbt(covstats='$1')

# set up plotting enviro
X11()
par(bg='white')
colors <- rainbow(6)

# plot
plot(df)
		
# pick bins in a loop -- setup
bins <- list()
bin_names <- list()
more_bins <- 'y'
	
# begin loop
while (more_bins == 'y') {

	# pick one bin
	bin <- choosebin(df, slice=1, save=TRUE, file=paste0(prefix, '.manual_curation.bin', length(bins) + 1, '.txt'), num.points=4)

	bins[[length(bins) + 1]] <- bin
	bin_names <- c(bin_names, paste0('bin', length(bins)))

	points(bins[[length(bins)]], col=colors[length(bins)], slice=1)

	message('Select another bin [y/n]:')
	more_bins <- scan(file='stdin', nlines=1, quiet=TRUE, what='character')

}
legend('topleft', legend=bin_names, fill=colors[1:length(bin_names)])
dev.copy(png, paste0(prefix, '.manual_curation.bins.png'))
dev.off()
graphics.off()
save.image(paste0(prefix, '.manual_curation.Rsession'))
q(save='no')

" > $script_path

	echo $script_path
	
	# change into the directory
	cd $bin_dir

	echo $(pwd)

	# run R
	R -f manually_pick_bins.r

	cd $curr_dir

	shift
done
