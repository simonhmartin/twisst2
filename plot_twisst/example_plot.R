
################################# overview #####################################

# The main data produced by Twisst is a weights file which has columns for each
# topology and their number of observations of that topology within each
# genealogy. Weights files produced by Twisst also contain initial comment lines
# speficying the topologies.

# The other data file that may be of interest is the window data. That is, the
# chromosome/scaffold and start and end positions for each of the regions or
# windows represented in the weights file.

# Both of the above files can be read into R, manipulated and plotted however
# you like, but I have written some functions to make these tasks easier.
# These functions are provided in the script plot_twisst.R

################### load helpful plotting functions #############################

source("plot_twisst.R")

############################## input files ######################################

# Each genomic region (chromosome, contig etc) shgould have two files@
# A topocounts file giving the count of each topology (one per column)
# And an intervals file, giving the chromosome name, start and end position

# Here we just import one genomic region

intervals_file <- "examples/admix_hiILS_123_l5e6_r1e8_mu1e8.chr1.intervals.tsv.gz"

topocounts_file <- "examples/admix_hiILS_123_l5e6_r1e8_mu1e8.chr1.topocounts.tsv.gz"

################################# import data ##################################

# The function import.twisst reads the topology counts and intervals data files into a list object
# If there are multiple weights files, or a single file with different chromosomes/scaffolds/contigs
# in the window data file, these will be separated when importing.

twisst_data <- import.twisst(intervals_files=intervals_file,
                             topocounts_files=topocounts_file)

#Some additional arguments:

# topos_file=
# If you prefer to provide the topologies as a separate file (rather than as comment lines in the topocounts file)

# ignore_extra_columns=TRUE
# This option will ignore subtrees that did not match any of the defined topologies (usually due to polytomies).
# If you use this option, the weightings for every tree will sum to 1, but you may be throwing away some information
# If you use this option, you might want to set the min_subtrees option too.

# min_subtrees=100
# This option can be used in conjunction with ignore_extra_columns=TRUE.
# For trees with polytomies there may be few subtrees that match any of the defined topologies.
# These weightings would be less reliable (more noisy). 
# Any tree with too fewer subtrees considered than the defined number will be ignored.

# max_interval=10000
# In parts of the genome with bad data, the tree interval can be very large.
# You can simply exclude this information by setting the maximum interval.

# names=
# If you have multiple regions they will be named according to the chromosome by default,
# but you can give your own names if you prefer

# reorder_by_start=TRUE
# If the tree intervals are out of order, this option will reorder them by the start position


############################## combined plots ##################################
# there are a functions available to plot both the weightings and the topologies

#for all plots we will use the 15 colours that come with plot_twisst.R
# But we will reorder them according to the most abundant topology
topo_cols <- topo_cols[order(order(twisst_data$weights_overall_mean[1:length(twisst_data$topos)], decreasing=TRUE))]

#a summary plot shows all the topologies and a bar plot of their relative weightings
plot.twisst.summary(twisst_data, lwd=3, cex=0.7)

plot.twisst.summary.boxplot(twisst_data)


#or plot ALL the data across the chromosome(s)
# Note, this is not recommended if there are large numbers of windows.
# instead, it is recommended to first smooth the weghtings and plot the smoothed values
# There are three plotting modes to try
plot.twisst(twisst_data, mode=1, show_topos=TRUE, ncol_topos=15)
plot.twisst(twisst_data, mode=2, show_topos=TRUE, ncol_topos=15)
plot.twisst(twisst_data, mode=3, show_topos=TRUE, ncol_topos=15)


# make smooth weightings and plot those across chromosomes
twisst_data_smooth <- smooth.twisst(twisst_data, span_bp = 20000, spacing = 1000)
plot.twisst(twisst_data_smooth, mode=2, ncol_topos=15) #mode 2 overlays polygons, mode 3 would stack them


##################### individual plots: raw weights ############################

#plot raw data in "stepped" style, with polygons stacked.
#specify stepped style by providing a matrix of starts and ends for positions
par(mfrow = c(1,1), mar = c(4,4,1,1))
plot.weights(weights_dataframe=twisst_data$weights[[1]], positions=twisst_data$interval_data[[1]][,c("start","end")],
             line_cols=topo_cols, fill_cols=topo_cols, stacked=TRUE)

#plot raw data in stepped style, with polygons unstacked (stacked =FLASE)
#use semi-transparent colours for fill
plot.weights(weights_dataframe=twisst_data$weights[[1]], positions=twisst_data$interval_data[[1]][,c("start","end")],
             line_cols=topo_cols, fill_cols=paste0(topo_cols,80), stacked=FALSE)


#################### individual plots: smoothed weights ########################

#plot smoothed data with polygons stacked
plot.weights(weights_dataframe=twisst_data_smooth$weights[[1]], positions=twisst_data_smooth$pos[[1]],
             line_cols=topo_cols, fill_cols=topo_cols, stacked=TRUE)

#plot smoothed data with polygons unstacked
plot.weights(weights_dataframe=twisst_data_smooth$weights[[1]], positions=twisst_data_smooth$pos[[1]],
             line_cols=topo_cols, fill_cols=paste0(topo_cols,80), stacked=FALSE)



#################### subset to only the most abundant topologies #################

#get list of the most abundant topologies (top 2 in this case)
top2_topos <- order(twisst_data$weights_overall_mean, decreasing=T)[1:2]

#subset twisst object for these
twisst_data_top2topos <- subset.twisst.by.topos(twisst_data, top2_topos)
#this can then be used in all the same plotting functions above.

######################## subset to only specific regions #########################

#regions to keep (more than one can be specified)
regions <- c("chr1")

#subset twisst object for these
twisst_data_chr1 <- subset.twisst.by.regions(twisst_data, regions)
#this can then be used in all the same plotting functions above.


########################### plot topologies using Ape ##########################
#unroot trees if you want to
# for (i in 1:length(twisst_data$topos)) twisst_data$topos[[i]] <- ladderize(unroot(twisst_data$topos[[i]]))

par(mfrow = c(3,length(twisst_data$topos)/3), mar = c(1,1,2,1), xpd=NA)
for (n in 1:length(twisst_data$topos)){
  plot.phylo(twisst_data$topos[[n]], type = "cladogram", edge.color=topo_cols[n], edge.width=5, rotate.tree = 90, cex = 1, adj = .5, label.offset=.2)
  mtext(side=3,text=paste0("topo",n))
  }


