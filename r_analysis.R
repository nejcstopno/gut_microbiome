library(vegan)
library(ggplot2)
library(tidyr)
library(reshape2)
library(dplyr)
library(rstatix)

otu_gut <- read.table('OTU_table.txt',sep='\t', header=T, row.names = 1)

map_gut <- read.csv('map.csv', row.names = 1, header=T)

tax_gut <- read.table('taxonomy.txt', sep='\t', header=T)

tax_gut_filt <- tax_gut[!grepl("Mitochondria", tax_gut$taxonomy),]
tax_gut_filt <- tax_gut_filt[!grepl("Chloroplast", tax_gut_filt$taxonomy),]
tax_gut_filt <- tax_gut_filt[!grepl("Unassigned", tax_gut_filt$taxonomy),]

tax_gut_filt <- tax_gut_filt %>%
  separate(taxonomy, into=c("Kingdom", "Phylum", "Class", 
                            "Order", "Family", "Genus", "Species"), sep=";", remove=T)
tax_gut_filt[2:8] <- lapply(tax_gut_filt[2:8], function(x) gsub(".*__", "", x))

otus_gut <- tax_gut_filt$OTUid
otu_gut_filt <- otu_gut[rownames(otu_gut) %in% otus_gut,]
otu_gut <- otu_gut_filt
map_ferret <- map_gut[map_gut$Animal == 'Ferret',]

otu_ferret <- otu_gut[,colnames(otu_gut) %in% map_ferret$Sample.Name]

# Order the samples
otu_ferret <- otu_ferret[,order(colnames(otu_ferret))]
map_ferret=map_ferret[order(as.character(map_ferret$Sample.Name)),]
map_ferret$Sample.Name==colnames(otu_ferret)

#Rarefying data to the sample with lowest reads (17000)
rarecurve(t(otu_ferret), step=20, label = FALSE)
otu_ferret <- otu_ferret[rowSums(otu_ferret)>0,]
set.seed(003)
otu_ferret_rare <- t(rrarefy(t(otu_ferret), 17000)) 

s <- specnumber(otu_ferret_rare,MARGIN=2)
h <- vegan::diversity(t(otu_ferret_rare), "shannon")
pielou=h/log(s)

map.div <- map_ferret
map.div$Richness <- s
map.div$Shannon <- h
map.div$Pielou <- pielou 

map.alpha <- melt(map.div, id.vars=c("time","rep", "Animal", 'treatment', 'Sample.Name', 'cage'), 
                  measure.vars=c("Richness", "Shannon", "Pielou"))

alpha_diversity=ggplot(map.alpha,aes(y=value, x=as.factor(time), color=treatment))+
  theme_bw()+
  scale_color_manual(values = c('#00429d', '#93003a'), 
                     breaks = c('Control', 'Infected')) +
  geom_boxplot()+
  facet_wrap(~variable, scale="free_y")+
  xlab('Time (days)')

aov.richness <- aov(value ~ time+treatment+cage+time*treatment, #replace Richness with Shannon
                      data = map.alpha[map.alpha$variable == 'Pielou',]) 
aov.shannon <- aov(value ~ time+treatment+cage+time*treatment, #replace Richness with Shannon
                    data = map.alpha[map.alpha$variable == 'Shannon',])
aov.pielou <- aov(value ~ time+treatment+cage+time*treatment, #replace Richness with Shannon
                    data = map.alpha[map.alpha$variable == 'Pielou',])
summary(aov.richness)
summary(aov.shannon)
summary(aov.pielou)

#' no stat differences between treatment and control when considering richness but
#' strong difference with Shannon (p=0.0001) and Pielou (p=0.0001)

stat.test.rich <- map.alpha %>%
  group_by(time) %>%
  filter(variable == 'Richness') %>%
  t_test(value ~ treatment) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test.rich

stat.test.shannon <- map.alpha %>%
  group_by(time) %>%
  filter(variable == 'Shannon') %>%
  t_test(value ~ treatment) %>%
  adjust_pvalue(method = "BH") 
stat.test.shannon

stat.test.pielou <- map.alpha %>%
  group_by(time) %>%
  filter(variable == 'Pielou') %>%
  t_test(value ~ treatment) %>%
  adjust_pvalue(method = "BH") 
stat.test.pielou

#PCoA full dataset (time 1 and 3)
gut.dist <- vegdist(t(otu_ferret_rare), method="bray")
gut.pcoa <- cmdscale(gut.dist, eig=TRUE)

gut_envfit <- envfit(gut.pcoa, map_ferret)

map.alpha$Axis1.BC <- gut.pcoa$points[,1]
map.alpha$Axis2.BC <- gut.pcoa$points[,2]
ax1.bc.otu <- gut.pcoa$eig[1]/sum(gut.pcoa$eig)
ax2.bc.otu <- gut.pcoa$eig[2]/sum(gut.pcoa$eig)

pcoa.BC <- ggplot(map.alpha, aes(x=Axis1.BC, y=Axis2.BC)) +
  theme_classic() +
  geom_point(aes(col=treatment, alpha=as.factor(time)), size=3)+
  scale_color_manual(values = c('#00429d', '#93003a')) + 
  #scale_alpha_continuous(range = c(0.1, 1)) + 
  theme(legend.position = c(.5,.4),
        legend.background = element_rect(fill="transparent"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10))+
  labs(x=paste('PCoA1 (',100*round(ax1.bc.otu,3),'%)',sep=''),
       y=paste('PCoA2 (',100*round(ax2.bc.otu,3),'%)', sep=''), 
       color='Treatment', alpha='Time (day)')

#are replicates different when you only compare the same times?
adonis(gut.dist~map_ferret$rep, strata=map_ferret$time) #no significancy with time
adonis(gut.dist~map_ferret$time) #no significancy with time
adonis(gut.dist~map_ferret$treatment)
adonis(gut.dist~map_ferret$treatment, strata = map_ferret$time) #no significancy with time
adonis(gut.dist~map_ferret$time, strata = map_ferret$treatment)
adonis(gut.dist~map_ferret$cage)

#PCoA only time 3
gut.time=otu_ferret_rare[,map_ferret$time == 3]
map.time=map_ferret[map_ferret$time == 3,]
gut.time <- gut.time[rowSums(gut.time)>0,]

gut.time.dist <- vegdist(t(gut.time), method="bray")
gut.time.pcoa <- cmdscale(gut.time.dist, eig=TRUE)
gut_envfit <- envfit(gut.time.pcoa, map.time)

map.time$Axis1.BC <- gut.time.pcoa$points[,1]
map.time$Axis2.BC <- gut.time.pcoa$points[,2]
ax1.bc.otu <- gut.time.pcoa$eig[1]/sum(gut.time.pcoa$eig)
ax2.bc.otu <- gut.time.pcoa$eig[2]/sum(gut.time.pcoa$eig)

pcoa.time <- ggplot(map.time, aes(x=Axis1.BC, y=Axis2.BC)) +
  theme_classic() +
  geom_point(aes(col=treatment), size=3)+
  scale_color_manual(values = c('#00429d', '#93003a')) + 
  theme(legend.position = c(.5,.4),
        legend.background = element_rect(fill="transparent"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10))+
  labs(x=paste('PCoA1 (',100*round(ax1.bc.otu,3),'%)',sep=''),
       y=paste('PCoA2 (',100*round(ax2.bc.otu,3),'%)', sep=''), 
       color='Treatment')

adonis(gut.time.dist~map.time$rep) 
adonis(gut.time.dist~map.time$treatment)
adonis(gut.time.dist~map.time$cage)
adonis(gut.time.dist~map.time$treatment, strata = map.time$cage)
adonis(gut.time.dist~map.time$cage, strata = map.time$treatment)

#' Microbiome diversity representation
gut.rel.abun <- decostand(otu_ferret_rare, method="total", MARGIN=2) #calculating relative abundance
gut_plot_taxa <- data.frame(OTUid=as.factor(rownames(gut.rel.abun)),gut.rel.abun) %>%
  gather(Sample.Name, abun, -OTUid) %>%
  left_join(map_ferret, by='Sample.Name') %>%
  left_join(tax_gut_filt, by='OTUid') %>%
  group_by(time,treatment,Phylum,Class,Family) %>%
  summarise(n_abun=sum(abun),
            n_samples=length(unique(Sample.Name)),
            rel_abun=n_abun/n_samples)

diversityPlot=ggplot(gut_plot_taxa,
       aes(x=as.factor(treatment), y=rel_abun, fill=Class)) +
  geom_bar(stat = 'identity')+
  theme_bw()+
  labs(y="Relative abundace", x= NULL) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.text=element_text(size=6),
        legend.position = 'bottom') +
  facet_grid(~time)
  #guides(fill=guide_legend(ncol=3))

data.frame(OTUid=as.factor(rownames(gut.rel.abun)), gut.rel.abun) %>%
  gather(Sample.Name, abun, -OTUid) %>%
  left_join(map_ferret, by='Sample.Name') %>%
  left_join(tax_gut_filt, by='OTUid') %>%
  filter(Genus =='Campylobacter') %>%
  ggplot(aes(x=treatment, y=abun, color=treatment)) +
  geom_boxplot()+
  labs(y='Relative abundance')+
  theme_bw()+
  scale_color_manual(values = c('#00429d', '#93003a')) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.text=element_text(size=6),
        legend.position = 'bottom') +
  guides(fill=guide_legend(ncol=3))+
  facet_grid(Genus~time)
  
library("gplots")
library("devtools")
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
library(pheatmap)

gut.rel.abun <- decostand(otu_ferret_rare, method="total", MARGIN=2) #calculating relative abundance
summary(rownames(gut.rel.abun)%in% tax_gut_filt$OTUid)

tax_gut_subset <- tax_gut_filt[tax_gut_filt$OTUid %in% rownames(gut.rel.abun),]
tax_gut_subset <- tax_gut_subset[order(tax_gut_subset$OTUid),]
gut.rel.abun <- gut.rel.abun[order(rownames(gut.rel.abun)),]

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

gut.rel.abun.high.abun <- gut.rel.abun[rowSums(gut.rel.abun)>.01,]
data_subset_norm <- t(apply(gut.rel.abun.high.abun, 1, cal_z_score))

my_sample_col <- data.frame(sample = rep('control', length(colnames(gut.rel.abun))))
my_sample_col$sample[colnames(data_subset_norm) %in% map_ferret$Sample.Name[map_ferret$treatment=='Infected']] <- 'infected'
row.names(my_sample_col) <- colnames(data_subset_norm)
time_col <- rep('Day 3', length(colnames(gut.rel.abun)))
time_col[colnames(data_subset_norm) %in% map_ferret$Sample.Name[map_ferret$time==1]] <- 'Day 1'
my_sample_col$time <- time_col

my_heatmap <- pheatmap(data_subset_norm, annotation_col = my_sample_col, cutree_rows = 2, cellwidth = 6, cellheight = 5, fontsize = 6, gaps_col = 6)
my_heatmap <- pheatmap(gut.rel.abun.high.abun, annotation_col = my_sample_col, cutree_cols = 3, cutree_rows = 2, cellwidth = 6, cellheight = 5, fontsize = 6)
data_subset_norm[rownames(data_subset_norm)=='OTU6',]
save_pheatmap_pdf <- function(x, filename, width=5, height=5) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(my_heatmap, "heatmap_z_score_v3.pdf")


#' Differentially abundant ZOTUs
library(DESeq2)
library(phyloseq)
tax_ferret=tax_gut_filt[tax_gut_filt$OTUid %in% rownames(otu_ferret),]
rownames(otu_ferret) %in% tax_ferret$OTUid
rownames(tax_ferret)=tax_ferret$OTUid
OTU = otu_table(otu_ferret, taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(tax_ferret))
SAM = sample_data(map_ferret)
physeq_ferret <- merge_phyloseq(phyloseq(OTU, TAX), SAM)

#Infected Vs Control
physeq_day3 = subset_samples(physeq, time == 3)
diagdds = phyloseq_to_deseq2(physeq_day3, ~ treatment)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05

sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq_ferret)[rownames(sigtab), ], "matrix"))
dim(sigtab)

ferret_otu <- rownames(sigtab)

data.frame(OTUid=as.factor(rownames(gut.rel.abun)), gut.rel.abun) %>%
  gather(Sample.Name, abun, -OTUid) %>%
  left_join(map_ferret, by='Sample.Name') %>%
  filter(OTUid %in% ferret_otu) %>%
  left_join(tax_gut_filt, by='OTUid') %>%
  ggplot(aes(x=treatment, y=abun, color=treatment)) +
  geom_boxplot(alpha=.5)+
  geom_point(size=2)+
  facet_wrap(~Genus*OTUid, scales = 'free_y')+
  labs(y='Relative abundance', x=NULL, color='OTU')+
  scale_color_manual(values = c('#00429d', '#93003a')) + 
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.text=element_text(size=6),
        legend.position = 'none',
        strip.text = element_text(size = 6)) +
  guides(fill=guide_legend(ncol=3))
