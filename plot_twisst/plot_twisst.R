simple.loess.predict <- function(x, y, span, new_x=NULL, weights = NULL, max = NULL, min = NULL, family=NULL){
    y.loess <- loess(y ~ x, span = span, weights = weights, family=family)
    if (is.null(new_x)) {y.predict <- predict(y.loess,x)}
    else {y.predict <- predict(y.loess,new_x)}
    if (is.null(min) == FALSE) {y.predict = ifelse(y.predict > min, y.predict, min)}
    if (is.null(max) == FALSE) {y.predict = ifelse(y.predict < max, y.predict, max)}
    y.predict
    }

smooth.df <- function(x, df, span, new_x = NULL, col.names=NULL, weights=NULL, min=NULL, max=NULL, family=NULL){
    if (is.null(new_x)) {smoothed <- df}
    else smoothed = df[1:length(new_x),]
    if (is.null(col.names)){col.names=colnames(df)}
    for (col.name in col.names){
        print(paste("smoothing",col.name))
        smoothed[,col.name] <- simple.loess.predict(x,df[,col.name],span = span, new_x = new_x, max = max, min = min, weights = weights, family=family)
        }
    smoothed
    }

smooth.weights <- function(interval_positions, weights_dataframe, span, new_positions=NULL, interval_lengths=NULL){
    weights_smooth <- smooth.df(x=interval_positions,df=weights_dataframe, weights = interval_lengths,
                                span=span, new_x=new_positions, min=0, max=1)
    
    #return rescaled to sum to 1
    weights_smooth <- weights_smooth / apply(weights_smooth, 1, sum)
    
    weights_smooth[is.na(weights_smooth)] <- 0
    
    weights_smooth
    }


stack <- function(mat){
    upper <- t(apply(mat, 1, cumsum))
    lower <- upper - mat
    list(upper=upper,lower=lower)
    }

interleave <- function(x1,x2){
    output <- vector(length= length(x1) + length(x2))
    output[seq(1,length(output),2)] <- x1
    output[seq(2,length(output),2)] <- x2
    output
    }


sum_df_columns <- function(df, columns_list){
    new_df <- df[,0]
    for (x in 1:length(columns_list)){
        if (length(columns_list[[x]]) > 1) new_df[,x] <- apply(df[,columns_list[[x]]], 1, sum, na.rm=T)
        else new_df[,x] <- df[,columns_list[[x]]]
        if (is.null(names(columns_list)[x]) == FALSE) names(new_df)[x] <- names(columns_list)[x]
        }
    new_df
    }


plot.weights <- function(weights_dataframe,positions=NULL,line_cols=NULL,fill_cols=NULL,density=NULL,lwd=1,xlim=NULL,ylim=c(0,1),stacked=FALSE,
                                        ylab="Weighting", xlab = "Position", main="",xaxt=NULL,yaxt=NULL,bty="n", add=FALSE){
    #get x axis
    x = positions
    #if a two-column matrix is given - plot step-like weights with start and end of each interval    
    if (dim(as.matrix(x))[2]==2) {
        x = interleave(positions[,1],positions[,2])
        yreps=2
        }
    else {
        if (is.null(x)==FALSE) x = positions
        else x = 1:nrow(weights_dataframe)
        yreps=1
        }
    
    #set x limits
    if(is.null(xlim)) xlim = c(min(x), max(x))
    
    #if not adding to an old plot, make a new plot
    if (add==FALSE) plot(0, pch = "", xlim = xlim, ylim=ylim, ylab=ylab, xlab=xlab, main=main,xaxt=xaxt,yaxt=yaxt,bty=bty)
    
    if (stacked == TRUE){
        y_stacked <- stack(weights_dataframe)
        for (n in 1:ncol(weights_dataframe)){
            y_upper = rep(y_stacked[["upper"]][,n],each=yreps)
            y_lower = rep(y_stacked[["lower"]][,n],each = yreps)
            polygon(c(x,rev(x)),c(y_upper, rev(y_lower)), col = fill_cols[n], density=density[n], border=NA)
            }
        }
    else{
        for (n in 1:ncol(weights_dataframe)){
            y = rep(weights_dataframe[,n],each=yreps)
            polygon(c(x,rev(x)),c(y, rep(0,length(y))), col=fill_cols[n], border=NA,density=density[n])
            lines(x,y, type = "l", col = line_cols[n],lwd=lwd)
            }
        }
    }

options(scipen = 7)

#Heres a set of 15 colourful colours from https://en.wikipedia.org/wiki/Help:Distinguishable_colors
topo_cols <- c(
"#0075DC", #Blue
"#2BCE48", #Green
"#FFA405", #Orpiment
"#5EF1F2", #Sky
"#FF5005", #Zinnia
"#005C31", #Forest
"#00998F", #Turquoise
"#FF0010", #Red
"#9DCC00", #Lime
"#003380", #Navy
"#F0A3FF", #Amethyst
"#740AFF", #Violet
"#426600", #Quagmire
"#C20088", #Mallow
"#94FFB5") #Jade




########### Below are some more object-oriented tools for working with standard twisst output files

library(ape)
library(tools)

#a function that imports topocounts and computes weights 
import.twisst <- function(topocounts_files, intervals_files=NULL, split_by_chrom=TRUE, reorder_by_start=FALSE, na.rm=TRUE, max_interval=Inf,
                          lengths=NULL, topos_file=NULL, ignore_extra_columns=FALSE, min_subtrees=1, recalculate_mid=FALSE, names=NULL){
    l = list()
    
    if (length(intervals_files) > 1){
        print("Reading topocounts and interval data")
        l$interval_data <- lapply(intervals_files, read.table ,header=TRUE)
        l$topocounts <- lapply(topocounts_files, read.table, header=TRUE)
        if (is.null(names) == FALSE) names(l$interval_data) <- names(l$topocounts) <- names
        }
    
    if (length(intervals_files) == 1){
        print("Reading topocounts and interval data")
        l$interval_data <- list(read.table(intervals_files, header=TRUE))
        l$topocounts <- list(read.table(topocounts_files, header=TRUE))
        if (split_by_chrom == TRUE){
            l$topocounts <- split(l$topocounts[[1]], l$interval_data[[1]][,1])
            l$interval_data <- split(l$interval_data[[1]], l$interval_data[[1]][,1])
            }
        }
    
    if (is.null(intervals_files) == TRUE) {
        print("Reading topocounts")
        l$topocounts <- lapply(topocounts_files, read.table, header=TRUE)
        n <- nrow(l$topocounts[[1]])
        l$interval_data <- list(data.frame(chrom=rep(0,n), start=1:n, end=1:n))
        if (is.null(names) == FALSE) names(l$interval_data) <- names
        }
    
    l$n_regions <- length(l$topocounts)
    
    if (is.null(names(l$interval_data)) == TRUE) {
        names(l$interval_data) <- names(l$topocounts) <- paste0("region", 1:l$n_regions)
        }
    
    print(paste("Number of regions:", l$n_regions))
    
    if (reorder_by_start==TRUE & is.null(intervals_files) == FALSE){
        print("Reordering")
        orders = sapply(l$interval_data, function(df) order(df[,2]), simplify=FALSE)
        l$interval_data <- sapply(names(orders), function(x) l$interval_data[[x]][orders[[x]],], simplify=F)
        l$topocounts <- sapply(names(orders), function(x) l$topocounts[[x]][orders[[x]],], simplify=F)
        }
    
    print("Getting topologies")
    
    #attempt to retrieve topologies
    l$topos=NULL
    #first, check if a topologies file is provided
    if (is.null(topos_file) == FALSE) {
        l$topos <- read.tree(file=topos_file)
        if (is.null(names(l$topos)) == TRUE) names(l$topos) <- names(l$topocounts[[1]])
        }
    else{
        #otherwise we try to retrieve topologies from the (first) topocounts file
        n_topos = ncol(l$topocounts[[1]]) - 1
        topos_text <- read.table(topocounts_files[1], nrow=n_topos, comment.char="", sep="\t", as.is=T)[,1]
        try(l$topos <- read.tree(text = topos_text))
        try(names(l$topos) <- sapply(names(l$topos), substring, 2))
        }
    
    print("Cleaning data")
    
    if (ignore_extra_columns == TRUE & is.null(l$topos)==FALSE) {
        for (i in 1:l$n_regions){
            #columns that are unwanted
            l$topocounts[[i]] <- l$topocounts[[i]][,1:length(l$topos)]
            }
        }
    
    if (na.rm==TRUE){
        for (i in 1:l$n_regions){
            #remove rows containing NA values
            row_sums = apply(l$topocounts[[i]],1,sum)
            good_rows = which(is.na(row_sums) == F & row_sums >= min_subtrees &
                              l$interval_data[[i]]$end - l$interval_data[[i]]$start + 1 <= max_interval)
            l$topocounts[[i]] <- l$topocounts[[i]][good_rows,]
            l$interval_data[[i]] = l$interval_data[[i]][good_rows,]
            }
        }
    
    print("Computing summaries")
    
    l$weights <- sapply(l$topocounts, function(raw) raw/apply(raw, 1, sum), simplify=FALSE)
    
    l$weights_mean <- t(sapply(l$weights, apply, 2, mean, na.rm=T))
    
    #weighting per region as a total. This will be used for getting the overall mean
    weights_totals <- apply(t(sapply(l$weights, apply, 2, sum, na.rm=T)), 2, sum)
    
    l$weights_overall_mean <- weights_totals / sum(weights_totals)
    
    if (is.null(lengths) == TRUE) l$lengths <- sapply(l$interval_data, function(df) tail(df$end,1), simplify=TRUE)
    else l$lengths = lengths
    
    
    for (i in 1:length(l$interval_data)) {
        if (is.null(l$interval_data[[i]]$mid) == TRUE | recalculate_mid == TRUE) {
            l$interval_data[[i]]$mid <- (l$interval_data[[i]]$start + l$interval_data[[i]]$end)/2
            }
        }
    
    l$pos=sapply(l$interval_data, function(df) df$mid, simplify=FALSE)
    
    l
    }


smooth.twisst <- function(twisst_object, span=0.05, span_bp=NULL, new_positions = NULL, spacing=NULL) {
    l=list()
    
    l$topos <- twisst_object$topos
    
    l$n_regions <- twisst_object$n_regions
    
    l$weights <- list()
    
    l$lengths = twisst_object$lengths

    l$pos <- list()
    
    for (i in 1:l$n_regions){
        if (is.null(span_bp) == FALSE) span <- span_bp/twisst_object$length[[i]]
        
        if (is.null(new_positions) == TRUE){
            if (is.null(spacing) == TRUE) spacing <- twisst_object$length[[i]]*span*.1
            new_positions <- seq(twisst_object$pos[[i]][1], tail(twisst_object$pos[[i]],1), spacing)
            }
        
        l$pos[[i]] <- new_positions
        
        interval_lengths <- twisst_object$interval_data[[i]][,3] - twisst_object$interval_data[[i]][,2] + 1
        
        l$weights[[i]] <- smooth.weights(twisst_object$pos[[i]], twisst_object$weights[[i]], new_positions = new_positions, span = span, interval_lengths = interval_lengths)
        }
    
    names(l$weights) <- names(twisst_object$weights)
    
    l
    }

is.hex.col <- function(string){
    strvec <- strsplit(string, "")[[1]]
    if (strvec[1] != "#") return(FALSE)
    if (length(strvec) != 7 & length(strvec) != 9) return(FALSE)
    for (character in strvec[-1]){
        if (!(character %in% c("0","1","2","3","4","5","6","7","8","9","a","b","c","d","e","f","A","B","C","D","E","F"))) return(FALSE)
        }
    TRUE
    }

hex.transparency <- function(hex, transstring="88"){
    if (is.hex.col(hex)==FALSE){
        print("WARNING: colour not hexadecimal. Cannot modify transparency.")
        return(hex)
        }
    if (nchar(hex) == 7) return(paste0(hex, transstring))
    else {
        substr(hex,8,9) <- transstring
        return(hex)
        }
    }


plot.twisst <- function(twisst_object, show_topos=TRUE, rel_height=3, ncol_topos=NULL, regions=NULL, ncol_weights=1,
                        cols=topo_cols, tree_type="clad", tree_x_lim=c(0,5), xlim=NULL, ylim=NULL, xlab=NULL, xaxt=NULL,
                        mode=2, margins = c(4,4,2,2), concatenate=FALSE, gap=0, include_region_names=FALSE){
    
    #check if there are enough colours
    if (length(twisst_object$topos) > length(cols)){
        print("Not enough colours provided (option 'cols'), using rainbow instead")
        cols = rainbow(length(twisst_object$topos))
        }
    
    if (mode==5) {
        stacked=TRUE
        weights_order=1:length(twisst_object$topos)
        fill_cols = cols
        line_cols = NA
        lwd = 0
        }
    
    if (mode==4) {
        stacked=FALSE
        weights_order=order(apply(twisst_data$weights_mean, 2, mean, na.rm=T)[1:length(twisst_object$topos)], decreasing=T) #so that the largest values are at the back
        fill_cols = cols
        line_cols = NA
        lwd=par("lwd")
        }
    
    if (mode==3) {
        stacked=FALSE
        weights_order=order(apply(twisst_data$weights_mean, 2, mean, na.rm=T)[1:length(twisst_object$topos)], decreasing=T) #so that the largest values are at the back
        fill_cols = sapply(cols, hex.transparency, transstring="80")
        line_cols = cols
        lwd=par("lwd")
        }
    
    if (mode==2) {
        stacked=FALSE
        weights_order=order(apply(twisst_data$weights_mean, 2, mean, na.rm=T)[1:length(twisst_object$topos)], decreasing=T) #so that the largest values are at the back
        fill_cols = sapply(cols, hex.transparency, transstring="80")
        line_cols = NA
        lwd=par("lwd")
        }
    
    if (mode==1) {
        stacked=FALSE
        weights_order=1:length(twisst_object$topos)
        fill_cols = NA
        line_cols = cols
        lwd=par("lwd")
        }
    
    if (is.null(regions)==TRUE) regions <- 1:twisst_object$n_regions
            
    if (concatenate == TRUE) ncol_weights <- 1
    
    if (show_topos==TRUE){
        n_topos <- length(twisst_object$topos)
        
        if (is.null(ncol_topos)) ncol_topos <- n_topos
        
        #if we have too few topologies to fill the spaces in the plot, we can pad in the remainder
        topos_pad <- (n_topos * ncol_weights) %% (ncol_topos*ncol_weights) 
        
        topos_layout_matrix <- matrix(c(rep(1:n_topos, each=ncol_weights), rep(0, topos_pad)),
                                    ncol=ncol_topos*ncol_weights, byrow=T)
        }
    else {
        ncol_topos <- 1
        n_topos <- 0
        topos_layout_matrix <- matrix(NA, nrow= 0, ncol=ncol_topos*ncol_weights)
        }
    
    #if we have too few regions to fill the spaces in the plot, we pad in the remainder
    data_pad <- (length(regions)*ncol_topos) %% (ncol_topos*ncol_weights)
    
    if (concatenate==TRUE) weights_layout_matrix <- matrix(rep(n_topos+1,ncol_topos), nrow=1)
    else {
        weights_layout_matrix <- matrix(c(rep(n_topos+(1:length(regions)), each=ncol_topos),rep(0,data_pad)),
                                        ncol=ncol_topos*ncol_weights, byrow=T)
        }
    
    layout(rbind(topos_layout_matrix, weights_layout_matrix),
           height=c(rep(1, nrow(topos_layout_matrix)), rep(rel_height, nrow(weights_layout_matrix))))
    
    if (show_topos == TRUE){
        if (tree_type=="unrooted"){
            par(mar=c(1,2,3,2), xpd=NA)
            
            for (i in 1:n_topos){
                plot.phylo(twisst_object$topos[[i]], type = tree_type, edge.color=cols[i],
                        edge.width=5, label.offset=0.3, cex=1, rotate.tree = 90, x.lim=tree_x_lim)
                mtext(side=3,text=paste0("topo",i), cex=0.75)
                }
            }
        else{
            par(mar=c(1,1,1,1), xpd=NA)
            
            for (i in 1:n_topos){
                plot.phylo(twisst_object$topos[[i]], type = tree_type, edge.color=cols[i],
                        edge.width=5, label.offset=0.3, cex=1, rotate.tree = 0, x.lim=tree_x_lim)
                mtext(side=3,text=paste0("topo",i,"   "), cex=0.75)
                }
            }
        }
    
    if (is.null(ylim)==TRUE) ylim <- c(0,1)
    
    par(mar=margins, xpd=FALSE)
    
    if (concatenate == TRUE) {
        chrom_offsets = cumsum(twisst_object$lengths + gap) - (twisst_object$lengths + gap)
        chrom_ends <- chrom_offsets + twisst_object$lengths
        
        plot(0, pch = "",xlim = c(chrom_offsets[1],tail(chrom_ends,1)), ylim = ylim,
             ylab = "", yaxt = "n", xlab = xlab, , xaxt = "n", bty = "n", main = "")
        
        for (j in regions) {
            if (is.null(twisst_object$interval_data[[j]])) positions <- twisst_object$pos[[j]] + chrom_offsets[j]
            else positions <- twisst_object$interval_data[[j]][,c("start","end")] + chrom_offsets[j]
            plot.weights(twisst_object$weights[[j]][weights_order], positions, xlim=xlim,
                         fill_cols = fill_cols[weights_order], line_cols=line_cols[weights_order],lwd=lwd,stacked=stacked, add=T)
            }
        }
    else{
        for (j in regions){
            if (is.null(twisst_object$interval_data[[j]])) positions <- twisst_object$pos[[j]]
            else positions <- twisst_object$interval_data[[j]][,c("start","end")]
            plot.weights(twisst_object$weights[[j]][weights_order], positions, xlim=xlim, ylim = ylim,
                         xlab=xlab, fill_cols = fill_cols[weights_order], line_cols=line_cols[weights_order],lwd=lwd,stacked=stacked, xaxt=xaxt)
            if (include_region_names==TRUE) mtext(3, text=names(twisst_object$weights)[j], adj=0, cex=0.75)
            }
        }
    }


#function for plotting tree that uses ape to get node positions
draw.tree <- function(phy, x, y, x_scale=1, y_scale=1, method=1, direction="right",
                      col="black", col.label="black", add_labels=TRUE, add_symbols=FALSE,
                      label_offset = 1, symbol_offset=0, col.symbol="black",symbol_bg="NA",
                      pch=19, cex=NULL, lwd=NULL, label_alias=NULL){

    n_tips = length(phy$tip.label)

    if (direction=="right") {
        node_x = (node.depth(phy, method=method) - 1) * x_scale * -1
        node_y = node.height(phy) * y_scale
        label_x = node_x[1:n_tips] + label_offset
        label_y = node_y[1:n_tips]
        adj_x = 0
        adj_y = .5
        symbol_x = node_x[1:n_tips] + symbol_offset
        symbol_y = node_y[1:n_tips]
        }
    if (direction=="down") {
        node_y = (node.depth(phy, method=method) - 1) * y_scale * 1
        node_x = node.height(phy) * x_scale
        label_x = node_x[1:n_tips]
        label_y = node_y[1:n_tips] - label_offset
        adj_x = .5
        adj_y = 1
        symbol_x = node_x[1:n_tips]
        symbol_y = node_y[1:n_tips] - symbol_offset
        }
    
    #draw edges
    segments(x + node_x[phy$edge[,1]], y + node_y[phy$edge[,1]],
             x + node_x[phy$edge[,2]], y + node_y[phy$edge[,2]], col=col, lwd=lwd)
    
    if (is.null(label_alias) == FALSE) tip_labels <- label_alias[phy$tip.label]
    else tip_labels <- phy$tip.label
    
    if (add_labels=="TRUE") text(x + label_x, y + label_y, col = col.label, labels=tip_labels, adj=c(adj_x,adj_y),cex=cex)
    if (add_symbols=="TRUE") points(x + symbol_x, y + symbol_y, pch = pch, col=col.symbol, bg=symbol_bg)

    }

#code for plotting a summary barplot
plot.twisst.summary <- function(twisst_object, order_by_weights=TRUE, only_best=NULL, cols=topo_cols,
                                x_scale=0.12, y_scale=0.15, direction="right", col="black", col.label="black",
                                label_offset = 0.05, lwd=NULL, cex=NULL){
    
    #check if there are enough colours
    if (length(twisst_object$topos) > length(cols)){
        print("Not enough colours provided (option 'cols'), using rainbow instead")
        cols = rainbow(length(twisst_object$topos))
        }
    
    # Either order 1-15 or order with highest weigted topology first
    
    if (order_by_weights == TRUE) {
        ord <- order(twisst_object$weights_overall_mean[1:length(twisst_object$topos)], decreasing=T)
        if (is.null(only_best) == FALSE) ord=ord[1:only_best]
        }
    else ord <- 1:length(twisst_object$topos)
    
    N=length(ord)
    
    #set the plot layout, with the tree panel one third the height of the barplot panel
    layout(matrix(c(2,1)), heights=c(1,3))
    
    par(mar = c(1,4,.5,1))
    
    #make the barplot
    x=barplot(twisst_object$weights_overall_mean[ord], col = cols[ord],
            xaxt="n", las=1, ylab="Average weighting", space = 0.2, xlim = c(0.2, 1.2*N))

    #draw the trees
    #first make an empty plot for the trees. Ensure left and right marhins are the same
    par(mar=c(0,4,0,1))
    plot(0,cex=0,xlim = c(0.2, 1.2*N), xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(0,1), bty="n")
    
    #now run the draw.tree function for each topology. You can set x_scale and y_scale to alter the tree width and height.
    for (i in 1:length(ord)){
        draw.tree(twisst_object$topos[[ord[i]]], x=x[i]+.2, y=0, x_scale=x_scale, y_scale=y_scale,
                col=cols[ord[i]], label_offset=label_offset, cex=cex, lwd=lwd)
        }
    
    #add labels for each topology
    text(x,.9,names(twisst_object$topos)[ord],col=cols[ord])
    }


#code for plotting a summary boxplot
plot.twisst.summary.boxplot <- function(twisst_object, order_by_weights=TRUE, only_best=NULL, trees_below=FALSE, cols=topo_cols,
                                x_scale=0.12, y_scale=0.15, direction="right", col="black", col.label="black",
                                label_offset = 0.05, lwd=NULL, label_alias=NULL, cex=NULL, outline=FALSE,
                                cex.outline=NULL, lwd.box=NULL, topo_names=NULL){
    
    #check if there are enough colours
    if (length(twisst_object$topos) > length(cols)){
        print("Not enough colours provided (option 'cols'), using rainbow instead")
        cols = rainbow(length(twisst_object$topos))
        }
    
    # Either order 1-15 or order with highest weigted topology first
    
    if (order_by_weights == TRUE) {
        ord <- order(twisst_object$weights_overall_mean[1:length(twisst_object$topos)], decreasing=T)
        if (is.null(only_best) == FALSE) ord=ord[1:only_best]
        }
    else ord <- 1:length(twisst_object$topos)
    
    N=length(ord)
    
    #set the plot layout, with the tree panel one third the height of the barplot panel
    if (trees_below==FALSE) layout(matrix(c(2,1)), heights=c(1,3))
    else layout(matrix(c(1,2)), heights=c(3,1))
    
    
    par(mar = c(1,4,.5,1), bty="n")
    
    #make the barplot
    boxplot(as.data.frame(rbindlist(twisst_object$weights))[,ord], col = cols[ord],
            xaxt="n", las=1, xlim = c(.5, N+.5), ylab="Weighting", outline=outline, cex=cex.outline, lwd=lwd.box)
    
    #draw the trees
    #first make an empty plot for the trees. Ensure left and right marhins are the same
    par(mar=c(0,4,0,1))
    plot(0,cex=0, xlim = c(.5, N+.5), xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(0,1), bty="n")
    
    #now run the draw.tree function for each topology. You can set x_scale and y_scale to alter the tree width and height.
    for (i in 1:N){
        draw.tree(twisst_object$topos[[ord[i]]], i+.1, y=0, x_scale=x_scale, y_scale=y_scale,
                col=cols[ord[i]], label_offset=label_offset, cex=cex, lwd=lwd, label_alias=label_alias)
        }
    
    if (is.null(topo_names)==TRUE) topo_names <- names(twisst_object$topos)
    
    #add labels for each topology
    text(1:N,.9,topo_names[ord],col=cols[ord], cex=cex)
    }

#function for subsetting the twisst object by a set of topologies
subset.twisst.by.topos <- function(twisst_object, topos){
    l <- list()
    regions <- names(twisst_object$weights)
    l$interval_data <- twisst_object$interval_data
    l$n_regions <- twisst_object$n_regions
    l$lengths <- twisst_object$lengths
    l$pos <- twisst_object$pos
    l$topocounts <- sapply(regions, function(region) twisst_data$topocounts[[region]][,topos], simplify=F)
    l$weights <- sapply(regions, function(region) twisst_object$weights[[region]][,topos], simplify=F)
    l$weights_mean <- l$weights_mean[,topos]
    l$weights_overall_mean <- twisst_object$weights_overall_mean[topos]
    l$topos <- twisst_object$topos[topos]
    l
    }

#function for subsetting the twisst object by specific regions
subset.twisst.by.regions <- function(twisst_object, regions){
    l <- list()
    regions <- names(twisst_object$weights[regions])
    l$interval_data <- twisst_object$interval_data[regions]
    l$n_regions <- length(regions)
    l$lengths <- twisst_object$lengths[regions]
    l$pos <- twisst_object$pos[regions]
    l$topocounts <- twisst_object$topocounts[regions]
    l$weights <- twisst_object$weights[regions]
    l$weights_mean <- twisst_object$weights_mean[regions]
    weights_totals <- apply(t(sapply(l$weights, apply, 2, sum, na.rm=T)), 2, sum)
    l$weights_overall_mean <- weights_totals / sum(weights_totals)
    l$topos <- twisst_object$topos
    l
    }

rbindlist <- function(l){
    df <- l[[1]]
    if (length(l) > 1){
        for (i in 2:length(l)){
            df <- rbind(df, l[[i]])
            }
        }
    df
    }


palettes = list(
sashamaps = c( #https://sashamaps.net/docs/resources/20-colors/
'#4363d8', #azure blue
'#3cb44b', #grass green
'#ffe119', #custard yellow
'#e6194b', #pomergranate red
'#f58231', #chalk orange
'#911eb4', #violet
'#46f0f0', #cyan
'#f032e6', #hot pink
'#000075', #navy blue
'#fabebe', #pale blush
'#008080', #teal
'#e6beff', #mauve
'#9a6324', #coconut bown
'#800000', #red leather
'#aaffc3', #light sage
'#808000', #olive
'#ffd8b1', #white skin
'#bcf60c', #green banana
'#808080', #grey
'#fffac8'), #beach sand

alphabet = c( #https://en.wikipedia.org/wiki/Help:Distinguishable_colors
"#0075DC", #Blue
"#FFA405", #Orpiment
"#2BCE48", #Green
"#003380", #Navy
"#F0A3FF", #Amethyst
"#990000", #Wine
"#9DCC00", #Lime
"#8F7C00", #Khaki
"#FFE100", #Yellow
"#C20088", #Mallow
"#FFCC99", #Honeydew
"#5EF1F2", #Sky
"#00998F", #Turquoise
"#FF5005", #Zinnia
"#4C005C", #Damson
"#426600", #Quagmire
"#005C31", #Forest
"#FFA8BB", #Pink
"#740AFF", #Violet
"#E0FF66", #Uranium
"#993F00", #Caramel
"#94FFB5", #Jade
"#FFFF80", #Xanthin
"#FF0010", #Red
"#191919", #Ebony
"#808080"), #Iron


krzywinski = c( #https://mk.bcgsc.ca/colorblind/palettes.mhtml#15-color-palette-for-colorbliness
"#68023F", #104   2  63 imperial purple, tyrian purple, nightclub, pompadour, ribbon, deep cerise, mulberry, mulberry wood, pansy purple, merlot
"#008169", #  0 129 105 deep sea, generic viridian, observatory, deep sea, elf green, deep green cyan turquoise, tropical rain forest, tropical rain forest, blue green, elf green
"#EF0096", #239   0 150 vivid cerise, persian rose, fashion fuchsia, hollywood cerise, neon pink, luminous vivid cerise, shocking pink, deep pink, deep pink, fluorescent pink
"#00DCB5", #  0 220 181 aquamarine, aqua marine, caribbean green, aqua, eucalyptus, vivid opal, brilliant turquoise, shamrock, shamrock, caribbean green
"#FFCFE2", #255 207 226 light pink, azalea, classic rose, pale pink, classic rose, pastel pink, pink lace, pig pink, orchid pink, chantilly
"#003C86", #  0  60 134 flat medium blue, royal blue, bay of many, bondi blue, congress blue, cobalt, elvis, darkish blue, submerge, yale blue
"#9400E6", #148   0 230 vivid purple, violet, purple, veronica, purple, electric purple, purple, vivid mulberry, vivid purple, vivid violet
"#009FFA", #  0 159 250 azure, luminous vivid cornflower blue, vivid cornflower blue, brilliant azure, bleu de france, dark sky blue, brilliant azure, united nations blue, light brilliant cobalt blue, cornflower
"#FF71FD", #255 113 253 blush pink, shocking pink, ultra pink, pink flamingo, fuchsia pink, light brilliant orchid, pink flamingo, candy pink, light brilliant magenta, light magenta
"#7CFFFA", #124 255 250 light brilliant cyan, electric blue, dark slate grey, very light cyan, very light opal, bright cyan, dark slate grey, brilliant cyan, light brilliant opal, bright light blue
"#6A0213", #106   2  19 burnt crimson, rosewood, claret, venetian red, dark tan, persian plum, prune, deep amaranth, deep reddish brown, crown of thorns
"#008607", #  0 134   7 india green, ao, green, office green, web green, green, emerald green, islamic green, forest green, forest green
"#F60239", #246   2  57 american rose, red, carmine, electric crimson, luminous vivid crimson, neon red, carmine red, scarlet, torch red, tractor red
"#00E307", #  0 227   7 vivid sap green, vivid emerald green, vivid green, vibrant green, vivid green, vivid harlequin, green, radioactive green, vivid malachite green, vivid pistachio
"#FFDC3D"), #255 220  61 gargoyle gas, filmpro lemon yellow, banana yellow, golden dream, bright sun, broom, banana split, golden dream, twentyfourseven, wild thing


safe = c( #from cartomap package
"#88CCEE",
"#CC6677",
"#DDCC77",
"#117733",
"#332288",
"#AA4499",
"#44AA99",
"#999933",
"#882255",
"#661100",
"#6699CC",
"#888888"),

cud = c(
"#E69F00",
"#56B4E9",
"#009E73",
"#F0E442",
"#0072B2",
"#D55E00",
"#CC79A7",
"#000000",
"#999999"))

show_palettes <- function(){
    par(mfrow = c(length(palettes), 1), mar = c(2,0,2,0))

    for (name in names(palettes)){
        n = length(palettes[[name]])
        plot(0, cex=0, xlim = c(1,n+1), ylim = c(0,1), xaxt="n", yaxt="n", bty="n", main=name)
        rect(1:n, 0, 2:(n+1), 1, col=palettes[[name]])
        axis(1, at=(1:n)+0.5, labels = palettes[[name]], tick=F)
        }
    }

