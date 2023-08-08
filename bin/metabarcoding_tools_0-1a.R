#### Custom functions
#Propagate taxonomie for tax_glom to work
propagate_incomplete_taxonomy <- function(phyloseq){
	taxranks <- colnames(tax_table(phyloseq))
	for (i in 2:length(taxranks)){
		tax_table(phyloseq)[tax_table(phyloseq)[,taxranks[i]]=="",taxranks[i]]<-paste(tax_table(phyloseq)[tax_table(phyloseq)[,taxranks[i]]=="",taxranks[i-1]],"_spc",sep="")
	}
	  return(phyloseq)
}


#Make taxa labels nice for plots
replace_tax_prefixes <- function(phyloseq){
  tmp_tax_table <- apply(tax_table(phyloseq), c(1, 2),function(y) gsub("^\\w:","",y))
  tmp_tax_table <- apply(tmp_tax_table, c(1, 2),function(y) gsub("_spc_.*","_spc",y))
  tmp_tax_table <- apply(tmp_tax_table, c(1, 2),function(y) gsub("_"," ",y))
  tmp_tax_table <- apply(tmp_tax_table, c(1, 2),function(y) gsub(";$","",y))
  tax_table(phyloseq)<- tmp_tax_table
  return(phyloseq)
}

#label low throughput samples
label_low_throughput <- function(phyloseq, threshold){
  sample_names(phyloseq)[sample_sums(phyloseq)<threshold]<-paste(sample_names(phyloseq)[sample_sums(phyloseq)<threshold],"|LT",threshold,sep="")
  return(phyloseq)
}

label_sample_by_host <- function(phyloseq, hostcol, projcol="", idcol=""){
  hostlist <- sample_data(phyloseq)[[hostcol]]
  projlist <- sample_data(phyloseq)[[projcol]]
  print(length(idcol))
  if(length(idcol)!=length(hostlist)){
   	uniqIDs <- sample(1:99999,length(hostlist) ,replace=F)
   	uniqIDs <- sprintf(paste(hostlist, projlist,"%05d", sep="|"), uniqIDs)
  }else{
  	uniqIDs <- sprintf(paste(hostlist, projlist,"%05d", sep="|"), idcol)
  }

  sample_names(phyloseq)<-uniqIDs
  return(phyloseq)
}

id_cont_lh <-function(phyloseq, distri=T, neg, pos){
  neg_samples=c()
  if (length(neg)>0){
    neg_samples= subset_samples(test, sample_names(test) %in% neg)
    neg_samples.rel = transform_sample_counts(neg_samples, function(x) x/sum(x))
    sam_samples= subset_samples(test, !(sample_names(test) %in% neg))
    sam_samples.rel = transform_sample_counts(sam_samples, function(x) x/sum(x))

        tail(sort(rowSums(otu_table(neg_samples))))
  }
  pos_samples=c()
  if (length(pos)>0){
    pos_samples=pos
  }


  ## checks
  ### distribution of continuous distribution over majority of samples
  rowSums(otu_table(test)>0)
  rowSums(otu_table(test))

  d1 <- ((rowSums(otu_table(sam_samples)>0))/length(otu_table(sam_samples)[1,]))
  d2 <- rowSums(otu_table(sam_samples.rel))
  d3 <- rowSums(otu_table(sam_samples.rel))/length(otu_table(sam_samples)[1,])

  d1n <- ((rowSums(otu_table(neg_samples)>0))/length(otu_table(neg_samples)[1,]))
  d2n <- rowSums(otu_table(neg_samples.rel))
  d3n <- rowSums(otu_table(neg_samples.rel))/length(otu_table(neg_samples.rel)[1,])

  #
  # names(d1)[d1-d1n<0]
  #
  # sort(d1-d1n)
  #
  # plot(d1n~d1)
  # plot(d2n~d2)

    plot(d1~d2 ) # i would expect contamination to be in many samples, but low abundances, i.e. high y and low x (upper left corner)
    plot(d3n~d3) # i would expect contamination to be low in d2 and high in d2n, NOT cross-contamination though between samples

    # points(d2n~d1n ,add=T, col="red")



  #### relation to low abundance samples

   abu_h <- subset_samples(sam_samples.rel, colSums(otu_table(sam_samples)) > 2500)
   abu_l <- subset_samples(sam_samples.rel, colSums(otu_table(sam_samples)) < 2500)

    plot(rowMeans(otu_table(abu_h))~rowMeans(otu_table(abu_l)), pch=19, alpha=0.3)
      # I would expect low throughput samples to have higher values for contamination taxa

    tail(sort(round(rowMeans(otu_table(abu_h))-rowMeans(otu_table(abu_l)), digits=4)))
    par(mar=c(4,15,1,1), mfrow=c(1,1))

    barplot(head(sort(round(rowMeans(otu_table(sam_samples.rel))-rowMeans(otu_table(neg_samples.rel)), digits=4)), n=20), horiz=T, las=2)
      # i would expect contamination to occur more in negative controls


        # plot(c(otu_table(its2.species.rel)[10,])~ c(colSums(otu_table(its2.species))), type="n", xlim=c(0,150000), ylim=c(0,1))

        # require(RColorBrewer)
        # col_vector = sample(grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)],length(otu_table(its2.species.rel)[,1]), replace=T)
        # #
    # for (i in 1:length(otu_table(its2.species.rel)[,1])){
    #   points(c(otu_table(its2.species.rel)[i,])~ c(colSums(otu_table(its2.species))),col=col_vector[i], add=T)
    # }
    #

    # check lt flag in name

  ### presence in negative controls
  ### presence in positive controls


}


map_interactions_trait <- function(phyloseq, traittable, traitcols="", speccol){
	otu <- otu_table(phyloseq)
	taxa <- taxa_names(phyloseq)
	samples <- sample_names(phyloseq)

	#traittable <- chemie
	#speccol <-"host"
	#traitcols = c(5:31)

	traits <- traittable[, traitcols]
	#traits$spc <- traittable[, speccol]

	tmp_matrix = matrix(ncol=length(traits[1,]), nrow=length(taxa), data=NA, byrow=F)
	rownames(tmp_matrix)<- taxa;		colnames(tmp_matrix)<- names(traits)

	results = data.frame(matrix(ncol=length(traits[1,]), nrow=0))
	colnames(results)<- names(traits)
	results_uncertainity = data.frame(matrix(ncol=length(traits[1,]), nrow=0))
	colnames(results_uncertainity)<- names(traits)
	results_uncertainity_weighted = data.frame(matrix(ncol=length(traits[1,]), nrow=0))
	colnames(results_uncertainity_weighted)<- names(traits)

	results_nomatch = c()

	for (i in 1:length(taxa)){
		matches <- traittable[, speccol] == taxa[i]
		if (sum(matches) == 1){
			tmp_matrix[taxa[i],] <- as.matrix(traits[matches,])#as.matrix(traits[matches,1:length(traits[1,])-1])
		}else if(sum(matches) > 1){
			# not sure what to do is multiple hits
		}else {
			results_nomatch =c(results_nomatch, taxa[i])
		}
	}
	template_matrix <- tmp_matrix

	for (j in 1:length(samples)){
		conc_matrix = template_matrix

#		if (sum(row.names(otu[,j])!=row.names(conc_matrix))>1){
#			print "Taxa inconsistency"
#		}else{
			prob_matrix <- matrix(ncol=length(tmp_matrix[1,]) , nrow=length(tmp_matrix[,1])  , data=otu[,j], byrow=F)
			mult_matrix <- conc_matrix * prob_matrix
			tmp_results <- as.data.frame(t(colSums	(mult_matrix, na.rm=T)))
			tmp_uncertainity <- as.data.frame(t(colSums	(is.na(mult_matrix)))/length(taxa))
			tmp_uncertainity_weighted <- colSums(as.numeric(!is.na(mult_matrix)) * prob_matrix)
			rownames(tmp_results) <- samples[j]
			rownames(tmp_uncertainity) <- samples[j]
			results = rbind(results, tmp_results)
			results_uncertainity = rbind(results_uncertainity, tmp_uncertainity)
			results_uncertainity_weighted = rbind(results_uncertainity_weighted, tmp_uncertainity_weighted)
#		} # end if taxa inconsistency


	} # end for j
	rownames(results_uncertainity_weighted)<- samples

	result_list=list()
	result_list$mapping <- results
	result_list$uncertainity <- results_uncertainity
	result_list$uncertainity_weighted <- results_uncertainity_weighted
	result_list$nomatch <- results_nomatch

  return(result_list)
}
