library(vegan)
library(ggplot2)
library(tidyverse)
library(reshape2)

setwd('~/OneDrive - Michigan State University/MSU/Ritam/')

otu_gut <- read.table('OTU_table.txt',sep='\t', header=T, row.names = 1)

map_gut <- read.csv('map.csv', row.names = 1, header=T)

tax_gut <- read.table('taxonomy.txt', sep='\t', header=T)

tax_gut_filt <- tax_gut[!grepl("Mitochondria", tax_gut$taxonomy),]
tax_gut_filt <- tax_gut_filt[!grepl("Chloroplast", tax_gut_filt$taxonomy),]
tax_gut_filt <- tax_gut_filt[!grepl("Unassigned", tax_gut_filt$taxonomy),]
tax_gut_filt <- tax_gut_filt[,-c(3,4)]

tax_gut_filt <- tax_gut_filt %>%
  separate(taxonomy, into=c("Kingdom", "Phylum", "Class", 
                            "Order", "Family", "Genus", "Species"), sep=";", remove=T)
tax_gut_filt[2:8] <- lapply(tax_gut_filt[2:8], function(x) gsub(".*__", "", x))

tax_gut_filt[] = lapply(tax_gut_filt, blank2na, na.strings=c('','NA','na','N/A','n/a','NaN','nan'))
lastValue <- function(x) tail(x[!is.na(x)], 1)
last_taxons<- apply(tax_gut_filt, 1, lastValue)
tax_gut_filt$last_taxon <- last_taxons
head(tax_gut_filt)
tax_gut_filt$final_names <- paste(tax_gut_filt$last_taxon, tax_gut_filt$OTUid, sep=' - ')

otus_gut <- tax_gut_filt$OTUid
otu_gut_filt <- otu_gut[rownames(otu_gut) %in% otus_gut,]
rownames(otu_gut_filt)
otu_gut <- otu_gut_filt
map_gut <- map_gut[map_gut$time!=1,]
otu_gut <- otu_gut[,colnames(otu_gut) %in% map_gut$Sample.Name]

# Order the samples
otu_gut <- otu_gut[,order(colnames(otu_gut))]
# Order the samples of the map the same way
map_gut=map_gut[order(as.character(map_gut$Sample.Name)),]
# Check to make sure they all match with each other
map_gut$Sample.Name==colnames(otu_gut)


#Rarefying data to the sample with lowest reads (31255)
rarecurve(t(otu_gut), step=20, label = FALSE)
otu_gut <- otu_gut[rowSums(otu_gut)>0,]
otu_gut_rare <- t(rrarefy(t(otu_gut), min(colSums(otu_gut)))) 

s <- specnumber(otu_gut_rare,MARGIN=2)
h <- vegan::diversity(t(otu_gut_rare), "shannon")
pielou=h/log(s)

map.div <- map_gut
map.div$Richness <- s
map.div$Shannon <- h
map.div$Pielou <- pielou 

map.alpha <- melt(map.div, id.vars=c("time","rep", "Animal", 'treatment', 'Sample.Name', 'cage'), 
                  measure.vars=c("Richness", "Shannon", "Pielou"))

ggplot(map.alpha[map.alpha$variable=='Richness',], aes(y=value, x=as.factor(time), color=treatment, group=cage))+
  xlab('Time (days)')+
  geom_boxplot()+
  facet_wrap(~Animal)+
  theme_bw()+
  theme(axis.text.x=element_text(angle = 45, hjust = 1))

#tidy(aov(value~cage* treatment, data=map.alpha[map.alpha$variable=="Richness" & map.alpha$Animal=="Ferret",]))
#no stat differences between tretment and control in Ferret
#PCoA
str(map_gut)
gut.dist <- vegdist(t(otu_gut_rare), method="bray")
gut.pcoa <- cmdscale(gut.dist, eig=TRUE)

gut_envfit <- envfit(gut.pcoa, map_gut)

ax1.gut<- gut.pcoa$eig[1]/sum(gut.pcoa$eig)
ax2.gut <- gut.pcoa$eig[2]/sum(gut.pcoa$eig)
ax3.gut <- gut.pcoa$eig[3]/sum(gut.pcoa$eig)

plot_color <- rep('hotpink2', nrow(map_gut))
plot_color[map_gut$treatment=='Infected'] <- "olivedrab3"

map_gut$Animal <- as.character(map_gut$Animal)
plot_color <- rep('black', nrow(map_gut))
plot_color[map_gut$Animal=='Chicken'] <- "red"

plot(gut.pcoa$points[,1], gut.pcoa$points[,2], bg=plot_color, pch=21,cex=2,
     xlab=paste("PCoA1: ",100*round(ax1.gut,3),"% var. explained",sep=""), 
     ylab=paste("PCoA2: ",100*round(ax2.gut,3),"% var. explained",sep=""))

#are replicates different when you only compare the same times?
adonis(gut.dist~map_gut$rep[map_gut$Animal=='Ferret'], strata=map_gut$time[map_gut$Animal=='Ferret']) #no significancy with time
adonis(gut.dist~map_gut$time[map_gut$Animal=='Ferret']) #no significancy with time
adonis(gut.dist~map_gut$Animal)
adonis(gut.dist~map_gut$treatment, strata = map_gut$Animal) #no significancy with time
adonis(gut.dist~map_gut$cage, strata = map_gut$treatment)

gut.rel.abun <- decostand(otu_gut_rare, method="total", MARGIN=2) #calculating relative abundance
gut_plot_taxa <- data.frame(OTUid=as.factor(rownames(gut.rel.abun)),gut.rel.abun) %>%
  gather(Sample.Name, abun, -OTUid) %>%
  left_join(map_gut[,c('Sample.Name','treatment', 'cage', 'Animal', 'time')], by='Sample.Name') %>%
  left_join(tax_gut_filt, by='OTUid') %>%
  group_by(Animal,treatment,Class,Family) %>%
  summarise(n_abun=sum(abun),
            n_samples=length(unique(Sample.Name)),
            rel_abun=n_abun/n_samples)

ggplot(gut_plot_taxa,
       aes(as.factor(treatment), rel_abun, fill=Class)) +
  #geom_point(pch=21, size=3) +
  geom_bar(stat = 'identity')+
  facet_wrap(~Animal) +
  theme_bw()+
  labs(y="Normalized Relative abundace", x= "Time") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.text=element_text(size=6),
        legend.position = 'bottom') +
  guides(fill=guide_legend(ncol=3)) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))

data.frame(OTUid=as.factor(rownames(gut.rel.abun)), gut.rel.abun) %>%
  gather(Sample.Name, abun, -OTUid) %>%
  left_join(map_gut[,c('Sample.Name','treatment', 'cage', 'Animal', 'time')], by='Sample.Name') %>%
  left_join(tax_gut_filt, by='OTUid') %>%
  filter(OTUid=='OTU6') %>%
  ggplot(aes(x=treatment, y=abun, color=treatment)) +
  geom_boxplot()+
  facet_wrap(~Animal)+
  labs(y='Relative abundance')+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.text=element_text(size=6),
        legend.position = 'bottom') +
  guides(fill=guide_legend(ncol=3))
  
#install.packages("gplots")
library("gplots")
library("devtools")
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

#install.packages('pheatmap')
library(pheatmap)
otu_gut_rare<- otu_gut_rare[rowSums(otu_gut_rare)>0,]
gut.rel.abun <- decostand(otu_gut_rare, method="total", MARGIN=2) #calculating relative abundance

heat <- t(scale(t(gut.rel.abun)))
pheatmap(heat)

tax_gut_filt$OTUid <- as.character(tax_gut_filt$OTUid)
tax_gut_subset <- tax_gut_filt[tax_gut_filt$OTUid %in% rownames(gut.rel.abun),]
tax_gut_subset <- tax_gut_subset[order(tax_gut_subset$OTUid),]
gut.rel.abun <- gut.rel.abun[order(rownames(gut.rel.abun)),]

rownames(gut.rel.abun) 
tax_gut_subset$OTUid
rownames(gut.rel.abun) <- tax_gut_subset$final_names

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

gut.rel.abun.high.abun <- gut.rel.abun[rowSums(gut.rel.abun)>.01,]
data_subset_norm <- t(apply(gut.rel.abun.high.abun, 1, cal_z_score))

my_sample_col <- data.frame(sample = rep('grey', length(colnames(gut.rel.abun))))
my_sample_col$sample[colnames(data_subset_norm) %in% map_gut$Sample.Name[map_gut$treatment=='Infected']] <- 'red'
my_sample_col <- within(my_sample_col, sample <- ifelse(is.na(sample), 'infected', 'control'))
row.names(my_sample_col) <- colnames(data_subset_norm)
animal_col <- rep('Chicken', length(colnames(gut.rel.abun)))
animal_col[colnames(data_subset_norm) %in% map_gut$Sample.Name[map_gut$Animal=='Ferret']] <- 'Ferret'
my_sample_col$animal <- animal_col

my_heatmap <- pheatmap(data_subset_norm, annotation_col = my_sample_col, cutree_rows = 2, cellwidth = 6, cellheight = 5, fontsize = 6, gaps_col = 6)
my_heatmap <- pheatmap(gut.rel.abun, annotation_col = my_sample_col, cutree_cols = 3, cutree_rows = 2, cellwidth = 6, cellheight = 5, fontsize = 6)

save_pheatmap_pdf <- function(x, filename, width=5, height=5) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(my_heatmap, "heatmap_z_score_v3.pdf")

#Campylobacter abundance 
tax_gut_v2 <- tax_gut_filt
rownames(tax_gut_v2) <- tax_gut_v2$OTUid
names(tax_gut_v2$OTUid) <- NULL
otu_gut_v2 <- otu_gut
map_gut_v2 <- map_gut
rownames(map_gut_v2) <- map_gut_v2$Sample.Name
names(map_gut_v2$Sample.Name) <- NULL

library(DESeq2)
library(phyloseq)
OTU = otu_table(otu_gut, taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(tax_gut_v2))
SAM = sample_data(map_gut_v2)
physeq <- merge_phyloseq(phyloseq(OTU, TAX), SAM)

#Infected Vs Control
physeq_ferret = subset_samples(physeq, Animal == "Ferret")
diagdds = phyloseq_to_deseq2(physeq_ferret, ~ treatment)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05

sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq_ferret)[rownames(sigtab), ], "matrix"))
dim(sigtab)

ferret_otu <- rownames(sigtab)
gut.rel.abun <- decostand(otu_gut_rare, method="total", MARGIN=2) #calculating relative abundance

data.frame(OTUid=as.factor(rownames(gut.rel.abun)), gut.rel.abun) %>%
  gather(Sample.Name, abun, -OTUid) %>%
  left_join(map_gut[,c('Sample.Name','treatment', 'cage', 'Animal', 'time')], by='Sample.Name') %>%
  filter(OTUid %in% ferret_otu,
         Animal == 'Ferret') %>%
  left_join(tax_gut_filt, by='OTUid') %>%
  ggplot(aes(x=treatment, y=abun, color=final_names)) +
  geom_boxplot(alpha=.5)+
  geom_point(size=2)+
  facet_wrap(~final_names, scales = 'free_y', labeller = label_wrap_gen(width=10))+
  labs(y='Relative abundance', x=NULL, color='OTU')+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.text=element_text(size=6),
        legend.position = 'none',
        strip.text = element_text(size = 6)) +
  guides(fill=guide_legend(ncol=3))
