# Loading in necessary libraries
library(phyloseq)
library(ggplot2)
library(bipartite)

# Setting working directory (check path)
setwd('/Users/ra39huv/tmp_tut/metabarcoding_workshop/data_ITS2')

# Custom functions inclusion (check path)
source('../bin/metabarcoding_tools_0-1a.R')


# Loading in data
## Taxonomy
data.tax <- tax_table(as.matrix(read.table("taxonomy.vsearch", header=T,row.names=1,fill=T,sep=",")))
# for bacteria only:
#data.tax <- data.tax[,1:5]

## Community table
data.otu <- otu_table(read.table("asv_table.merge.txt"), taxa_are_rows=T)

## Sample metadata (second line optional if sample names include "-"):
data.map <- 	sample_data(read.table("ITS2_samples.csv", header=T, row.names=1,  sep=";", fill=T))

## check metadata vs. samples in sequencing data consistency
sample_names(data.map )[!(sample_names(data.map ) %in% sample_names(data.otu))]
sample_names(data.otu )[!(sample_names(data.otu ) %in% sample_names(data.map))]

## merge the three tables to a single phylseq object
(data.comp <- merge_phyloseq(data.otu,data.tax,data.map ))
tail(tax_table(data.comp))

# preprocessing data pt.1 :
## given hierarchical classification options at the end, we have to propagate the taxonomy over taxonomic levels to not throw out stuff only classified to higher tax levels
data.comp <- propagate_incomplete_taxonomy(data.comp)
tail(tax_table(data.comp))

## filtering irrelevant taxa, 
## 16S data: 
#(data.comp.filter = subset_taxa(data.comp, kingdom=="d:Bacteria"))
#(data.comp.filter = subset_taxa(data.comp.filter, genus!="d:Bacteria_spc_spc_spc_spc"))
#(data.comp.filter = subset_taxa(data.comp.filter, family!="f:Chloroplast"))

## ITS2 data: 
(data.comp.filter = subset_taxa(data.comp, phylum=="p:Streptophyta"))
(data.comp.filter = subset_taxa(data.comp.filter, species!="p:Streptophyta_spc_spc_spc_spc"))
(data.comp.filter = subset_taxa(data.comp.filter, species!="k:Viridiplantae_spc_spc_spc_spc_spc"))
(data.comp.filter = subset_taxa(data.comp.filter, kingdom!=""))

## Make taxa labels nice for plots
data.comp.filter <- replace_tax_prefixes(data.comp.filter)

### Check the names
tail(tax_table(data.comp.filter))

## Multiple ASVs might represent the same species, here they are collated
# "genus" for 16s data, species" for ITS2 data
(data.species <- tax_glom(data.comp.filter,taxrank="genus"))
taxa_names(data.species) <- tax_table(data.species)[,"genus"]

# Transform to relative data
data.species.rel = transform_sample_counts(data.species, function(x) x/sum(x))

# low abundance filtering
otu_table(data.species.rel)[otu_table(data.species.rel)<0.01 ]<-0
otu_table(data.species)[otu_table(data.species.rel)<0.01 ]<-0
data.species.filter		= prune_taxa(taxa_sums(data.species)>0, data.species)
(data.species.rel.filter = prune_taxa(rowSums(otu_table(data.species.rel))>0, data.species.rel))

# First diversity and community metrics and graphs
# Distribution of major taxa, accumulated over all samples
par(mar=c(4,15,1,1), mfrow=c(1,1))
barplot(t(as.data.frame(sort(taxa_sums(data.species.rel.filter), decreasing=T)[1:20])), las=2, horiz=T)

## Distribution over samples
data.melt <- psmelt(data.species.rel.filter)
data.melt$Sample <- as.factor(data.melt$Sample)

dir.create("plots")

pdf("plots/sample_rel_abundance.pdf", width=12, height=15)
ggplot(data.melt, aes(OTU, Abundance, fill= family)) +
  facet_grid(Host ~ .)+
  theme_bw()+
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  geom_bar(stat="identity")+
  theme(legend.position="bottom")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),axis.title.x = element_text(family = "sans", size = 15)) + xlab("Taxa")
dev.off()

## Richness (Observed) / Shannon Diversity (replace x/col by group in metadata)
pdf("plots/sample_diversity.pdf", width=15, height=10)
plot_richness(data.species.filter,x= "Host" , col= "Host", measures=c("Observed","Shannon"))+theme_bw()+ theme(axis.text.x = element_text(angle = 60, hjust = 1) ) + geom_point(size=4) +geom_boxplot()
dev.off()

## Ordination
data.species.filter.nmds <-  ordinate(data.species.rel.filter, method="NMDS", "bray",k=4, trymax=200)

#data.species.rel.filter2 <- subset_samples(data.species.rel.filter, !sample_names(data.species.rel.filter) %in% c("p|L|48468","l|B|99952"))
data.species.filter.nmds <-  ordinate(data.species.rel.filter, method="DCA")

pdf("plots/sample_ordination.pdf", width=15, height=10)
plot_ordination(data.species.rel.filter, data.species.filter.nmds , color="Host")+geom_point(size=6)+theme_bw()#+geom_label(label=sample_names(data.species.rel.filter))
dev.off()


## Networks (replace id by group if you want to have them merged by metadata)
#netmat <- t(otu_table(data.species.rel.filter))
netmat <- otu_table(merge_samples(data.species.rel.filter,group="Host"))

### (optional) keep only major links, otherwise it often becomes overwhelming
otu_table(netmat)[otu_table(netmat)<0.05 ]<-0
netmat		= prune_taxa(taxa_sums(netmat)>0, netmat)

### plotting
pdf("plots/sample_network.pdf", width=15, height=5)
plotweb(data.frame(t(otu_table(netmat))))
dev.off()
