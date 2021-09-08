#!/usr/bin/env Rscript

plot_paired_sample_compositions = function(highest_sample,highest_sample_paired=NULL, otu_table, reads_info = NULL, covar_info=NULL){
  #' takes one or two sample names, and an otu table and plots the composition of the pair. Doesn't necessarily need to be paired samples, but will find the second member of the pair if one is not specified. Option to specify information abouut covariates and read abundances, in the specific format they were generated in. 
  
  highest_samp = highest_sample
  if (is.null(highest_sample_paired)) {
    highest_samp_paired = ifelse(endsWith(highest_samp, "A"), paste(substr(highest_samp, 1, nchar(highest_samp)-1), "B", sep = "" ), paste(substr(highest_samp, 1, nchar(highest_samp)-1), "A", sep = "" )) }
  
  #What percentage of the sample had been covariates? How many reads?
  if (!is.null(reads_info, covar_info)) {
    covar_abund = round(1-covar_info[covar_info$Sample==highest_samp,2],3)
    high_samp_reads = reads_info[reads_info$Sample==highest_samp,2]
    
    #If a paired sample exists, what is the relative abundance of covariate contaminants in it? how many reads?
  if (length(highest_samp_paired)==1) {
    covar_abund_paired =round(1-covar_info[covar_info$Sample==highest_samp_paired,2],3)
    high_samp_paired_reads = reads_info[reads_info$Sample==highest_samp_paired,2]}
  }
  
  maxsamp_comp = otu_table[rownames(otu_table) %in% c(highest_samp, highest_samp_paired),]
  
  maxsamp_comp = melt(t(maxsamp_comp))
  colnames(maxsamp_comp) = c("Taxa", "Sample", "Freq")
  maxsamp_comp$alpha = sapply(maxsamp_comp$Taxa, function(x) {ifelse(x==i, TRUE, FALSE)})
  if (length(covar_abund ==1)){
  maxsamp_comp$Covar_Abund = as.factor(sapply(maxsamp_comp$Sample, function(x) {ifelse(x==highest_samp, covar_abund, covar_abund_paired)}))}
  if (length(high_samp_reads) ==1) {
  maxsamp_comp$Reads =  as.factor(sapply(maxsamp_comp$Sample, function(x) {ifelse(x==highest_samp, paste("reads=",high_samp_reads, sep = ""), paste("reads=",high_samp_paired_reads, sep = ""))})) }
  
  
  maxsamp_comp = maxsamp_comp[maxsamp_comp$Freq>0,] #plot only values that are not 0
  
  top5tax = as.character(unique(maxsamp_comp[order(-maxsamp_comp$Freq),]$Taxa)[c(1:5)])
  
  mycols =colorRampPalette(brewer.pal(8, "Dark2")) #make a palette
  ncolors = nrow(maxsamp_comp) #we will use this to expand the palette above
  
  composition_plot = ggplot(maxsamp_comp, aes(x=Sample, y= Freq))+
    geom_bar(aes(fill = Taxa), stat="identity")+
    theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    scale_alpha_discrete(range=c(0.4, 0.90), guide = "none")+
    scale_fill_manual(name = "top 5 taxa" , breaks=top5tax, values= mycols(ncolors)) +
    scale_y_continuous(breaks=seq(0,1,0.2))+
    geom_text(aes(y = 1.1, label = paste("contam abund", Covar_Abund, sep = " ")), vjust = 0)+
    geom_text(aes(y = 1.2, label = Reads), vjust = 0)#+
  # theme(legend.position = "none")
  
  ## get the legend 
  tmp <- ggplot_gtable(ggplot_build(composition_plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) ==  "guide-box")
  legend <- tmp$grobs[[leg]]
  
  
  summary = plot_grid(composition_plot +theme(legend.position = "none"), arrangeGrob(legend), align = "v", nrow = 3, rel_heights = c(3/4, 1/4))
  
  return(summary)
}
  
plot_pcoa_axis_compositions = function(ordination_data, greater_than_threshold=NULL, less_than_threshold=NULL, axis, filtered_otu_table){
  #' takes a phyloseq ordination table, complete with metadata if you want, and makes a composition plot of the taxa above an axis value. useful for figuring out which taxa are driving structure in a pcoa
  
  rownames(filtered_otu_table) = make.names(rownames(filtered_otu_table))
  
  if (axis =="Axis.1"){
    if (!is.null(greater_than_threshold)) {axis_samps = ordination_data[ordination_data$Axis.1 > greater_than_threshold, "SampleName"]}
    if (!is.null(less_than_threshold)){axis_samps = ordination_data[ordination_data$Axis.1 < less_than_threshold, "SampleName"]}
    if (!is.null(less_than_threshold) & !is.null(greater_than_threshold)){
      axis_samps = ordination_data[ordination_data$Axis.1 > greater_than_threshold & ordination_data$Axis.1 < less_than_threshold, "SampleName"]}}
  else if (axis == "Axis.2"){
    if (!is.null(greater_than_threshold)) {axis_samps = ordination_data[ordination_data$Axis.2 > greater_than_threshold, "SampleName"]}
    if (!is.null(less_than_threshold)){axis_samps = ordination_data[ordination_data$Axis.2 < less_than_threshold, "SampleName"]}
    if (!is.null(less_than_threshold) & !is.null(greater_than_threshold)){
      axis_samps = ordination_data[ordination_data$Axis.2 > greater_than_threshold & ordination_data$Axis.2 < less_than_threshold,"SampleName"]}}
    
  if (is.null(greater_than_threshold) & is.null(less_than_threshold)){stop("pls supply threshold")}
  
  # some compatibility hacking, sorry
  
  axis_samps = make.names(as.character(axis_samps))
  
  #Dataframe of OTUs
  axis_comp = filtered_otu_table[rownames(filtered_otu_table) %in% axis_samps,]
  # colnames(axis_comp)[apply(axis_comp,1,which.max)] #What bug makes up the highest % of every sample? 
  
  axis_comp = melt(as.matrix(axis_comp))
  colnames(axis_comp) = c("ID", "Taxonomy", "Freq")
  
  #Plot
  top5tax = as.character(unique(axis_comp[order(-axis_comp$Freq),]$Taxonomy)[c(1:5)])
  mycols =colorRampPalette(brewer.pal(8, "Dark2")) #make a palette
  ncolors = nrow(axis_comp) #we will use this to expand the palette above
  
  axis_plot = ggplot(axis_comp)+
    geom_bar(aes(x=ID, y=Freq, fill=Taxonomy), stat="identity" )+
    # facet_grid(.~ID)+
    scale_fill_manual(name = "top 5 taxa" , breaks=top5tax, values= sample(mycols(ncolors))) +
    theme(axis.text.x = element_blank(), strip.text.x = element_blank(), axis.ticks.x = element_blank())
  
  tmp <- ggplot_gtable(ggplot_build(axis_plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) ==  "guide-box")
  legend <- tmp$grobs[[leg]]
  
  summary = plot_grid(axis_plot +theme(legend.position = "none"), arrangeGrob(legend), align = "v", nrow = 2, rel_heights = c(3/4, 1/4), axis = "b")
  return(summary)
}

coerce_to_phyloseq = function(filtered_otu_table, metadata_file){
  
  #' Takes an otu table and associtaed metadata file to turn into a phyloseq class obj.
  #' for the df, the samples must be rows and cols OTUs
  #' make sure sample names are held in the rownmaes.
  
  otumat = data.frame(t(filtered_otu_table)) # Remember to use the re-normalised dataset!

  #Taxa table
  
  all_tax = rownames(otumat)
  
  #Split OTU names into lists of each taxonomic level
  leg_lbl <- strsplit(all_tax , split='.',fixed=TRUE) # split at ".". Very unspecific for everything beyond Genus (field:6)
  
  #Initialise empty taxa table
  taxtable = data.frame(matrix(ncol=7,nrow = 0))
  colnames(taxtable) = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  num_to_domain <- c('D','P','C','O','F','G','S') # 6: Genus, 7:species
  
  
  for (i in 1:length(leg_lbl)){
    for (level in c(1:7)) {
      if(leg_lbl[[i]][level] != "__"){ # tax level seems to be defined    
        if(level != 7 ){
          # non-species level tax
          taxtable[i,level] <- paste(num_to_domain[level], leg_lbl[[i]][level], sep = "_")
        } else {
          # species level: concatenate all strings from species level (7:)
          taxtable[i,level] <- paste(num_to_domain[level], paste(leg_lbl[[i]][ -(1:6) ],collapse="_"), sep = "_")
        }
      }
    }
  }
  
  rownames(taxtable) = all_tax
  
  taxtbl = as.matrix(taxtable)
  
  # sample data
  
  meta_filtered = metadata_file
  rownames(meta_filtered) = make.names(rownames(meta_filtered))
  meta_filtered$Visit = as.factor(meta_filtered$Visit)
  meta_filtered$Kit = as.factor(meta_filtered$Kit)
  
  #Stitching into a phyloseq obj -- note that R hates variable names that start with a number, so it inserts an X in front
  
  OTU = otu_table(otumat, taxa_are_rows = TRUE)
  TAX = tax_table(taxtbl)
  sampledata = sample_data(meta_filtered)
  
  physeq = phyloseq(OTU, TAX, sampledata)
  
  return(physeq) }
 
  