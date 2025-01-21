
# Author: Alicia Garcia-Roldan

# Date: December 2023

# R v.3.6.3 and v.4.1.2

# SqueezeMeta v.1.6.1post1


# An amazing and cool script that selects all the high and medium quality bins (completeness >= 50%, contamination < 10%) in your sample and
# make a heatmap plot based on the abundance (copy number) of certain genes of each bin.


#########################
## BIOTIN TRANSPORTERS ##
#########################


## LIBRARY ##

library(SQMtools)   # v.1.6.2
library(ggplot2)    # v.3.4.0
library(paletteer)  # v.1.5.0


## DATASET ##

samples <- strsplit(Job1$misc$samples, '_SQM')  # Instead of writing, a cool way of having the samples names...hacker ;)


sample2realbin = list() # List of list of list of list...

for (i in samples){
  sample2realbin[[i]] = rownames(projs[[i]]$bins$table)
}


sample2fakename = list() # The fake name consists on sample_bin_no.fa ending

for (sample in names(sample2realbin)) {
  realbin = sample2realbin[[sample]][0:length(sample2realbin[[sample]])]
  sample2fakename[[sample]] = paste(sample, unlist(strsplit(realbin, ".fa"))[seq(1,length(unlist(strsplit(realbin, ".fa"))),
                                                                   by=2)], sep = "_") # Playing to have the perfect fakename
  
  # fake_name = sample2fakename[[sample]][0:length(sample2fakename[[sample]])]
  # real_name = sample2realbin[[sample]][0:length(sample2realbin[[sample]])]
  
}


sample2coolrealbin = list()
sample2coolfakebin = list()
coolfakebin2sample = list()

for (sample in names(sample2realbin)){
  bins = sample2realbin[[sample]]
  for(b in bins) {
    proj = projs[[sample]]
    comp = proj$bin$table[b,'Completeness']
    cont = proj$bin$table[b,'Contamination']
    if(!is.na(comp) & !is.na(cont) & comp >= 50 & cont < 10) { # We only want high and medium quality bins
      #sample2coolrealbin[[b]] = sample
      sample2coolrealbin[[sample]] = c(sample2coolrealbin[[sample]], b)
      sample2coolfakebin[[sample]] = c(sample2coolfakebin[[sample]], paste(sample,
                                                                           unlist(strsplit(b,".fa"))[seq(1,length(unlist(strsplit(b, ".fa"))),
                                                                                                                by=2)],
                                                                           sep="_")) # Playing to have the perfect fakename
      
      tocho = paste(sample, unlist(strsplit(b,".fa"))[seq(1,length(unlist(strsplit(b, ".fa"))),
                                                                                  by=2)], sep="_")
      
      coolfakebin2sample[[tocho]] = c(coolfakebin2sample[[tocho]], sample) # Making and invert list to ease the work in the future...
      
    }
  }
}


# caca = list()
# 
# for (sample in names(sample2coolrealbin)){
#   coolbins = sample2coolrealbin[[sample]]
#   fakebins = sample2coolfakebin[[sample]]
#   for (b in coolbins) {
#     for (c in fakebins) {
#       caca[[c]] = subsetBins(projs[[sample]], b)
#     }
#   }
# }


# Now we have to use subsetBins to work with the genes that we want, so here's an option...

new_bins = list()

for (sample in names(sample2coolrealbin)){
  coolbins = sample2coolrealbin[[sample]]
  fakebins = sample2coolfakebin[[sample]]
  for (i in 1:length(sample2coolrealbin[[sample]])){
    new_bins[[fakebins[i]]] = subsetBins(projs[[sample]], coolbins[i])
  }
}

# Another easier option would be a list with fakenames and realnames so you can call each name without using "i"



ma = matrix(0, nrow=length(coolfakebin2sample), ncol=6, dimnames= list(names(coolfakebin2sample), c('BioY (K03523)',
                                                                                                    'BioN (K16783)',
                                                                                                    'BioM (K16784)',
                                                                                                    'EcfT (K16785)',
                                                                                                    'EcfA1 (K16786)',
                                                                                                    'EcfA2 (K16787)'))) # Genes we want to work with

for(sample in names(sample2coolfakebin)) {
  fakebins = sample2coolfakebin[[sample]]
  # for (i in 1:length(sample2coolfakebin[[sample]])){
  for(b in fakebins) {
    x = new_bins[[b]]

    if ("K03523" %in% rownames(x$functions$KEGG$copy_number)){ # Sometimes the gene you are looking for is not in each bin so it turns out as an error
                                                               # with this statement, we avoid this uncomfortable situation...
      ma[b,'BioY (K03523)'] = ma[b,'BioY (K03523)'] + x$functions$KEGG$copy_number['K03523',]  # Value of copy number will be added to our matrix
    }
    if ("K16783" %in% rownames(x$functions$KEGG$copy_number)){
      ma[b,'BioN (K16783)'] = ma[b,'BioN (K16783)'] + x$functions$KEGG$copy_number['K16783',]
    }
    if ("K16784" %in% rownames(x$functions$KEGG$copy_number)){
      ma[b,'BioM (K16784)'] = ma[b,'BioM (K16784)'] + x$functions$KEGG$copy_number['K16784',]
    }
    if ("K16785" %in% rownames(x$functions$KEGG$copy_number)){
      ma[b,'EcfT (K16785)'] = ma[b,'EcfT (K16785)'] + x$functions$KEGG$copy_number['K16785',]
    }
    if ("K16786" %in% rownames(x$functions$KEGG$copy_number)){
      ma[b,'EcfA1 (K16786)'] = ma[b,'EcfA1 (K16786)'] + x$functions$KEGG$copy_number['K16786',]
    }
    if ("K16787" %in% rownames(x$functions$KEGG$copy_number)){
      ma[b,'EcfA2 (K16787)'] = ma[b,'EcfA2 (K16787)'] + x$functions$KEGG$copy_number['K16787',]
    }
    
}
}


## SAVE ##

write.table(ma, file = "Biotin_transporters_Job1_matrix.txt",sep = "\t")


# We want a cool plot, so we need to use ggplot2...what do we need to use ggplot? Yay! A beautiful dataframe <3

ma = as.data.frame(as.table((ma))) # Let's change our matrix into a dataframe

colnames(ma)<-c("Bins", "Genes", "Copy_number")


# Maybe some day we need to have a separate column with the sample to which our bins belong. Here we have two options:

# OPTION A: A little bit shabby...

# sample=c(rep("IC19'5", length(sample2coolfakebin[["IC19_5"]])),
#          rep("IC22", length(sample2coolfakebin[["IC22"]])),
#          rep("IC36_1", length(sample2coolfakebin[["IC361"]])),
#          rep("IC36_2", length(sample2coolfakebin[["IC362"]])),
#          rep("IC39_1", length(sample2coolfakebin[["IC391"]])),
#          rep("IC39_2", length(sample2coolfakebin[["IC392"]])))


# OPTION B: B of better ;)

# sample=c()
# 
# for (b in coolfakebin2sample){
#   sample <-c(sample, b)
# }
# 


# sample=rep(sample,6) # We have six genes, so we repeat it six times
# 
# ma <- cbind(ma, sample) # And we add it to our dataframe


## SAVE ##

write.table(ma, file = "Biotin_transporters_Job1_df.txt",sep = "\t")


## PLOTTING TIME!! ##  Manzanilla script is my God


Wessie=c(rep("#B40F20", length(sample2coolfakebin[["IC19_5"]])),
         rep("#E58601", length(sample2coolfakebin[["IC22"]])),
         rep("#46ACC8", length(sample2coolfakebin[["IC361"]])),
         rep("#82C7DA", length(sample2coolfakebin[["IC362"]])),
         rep("#E2D200", length(sample2coolfakebin[["IC391"]])),
         rep("#EDE466", length(sample2coolfakebin[["IC392"]]))) # My Wessie till the end of Job1

NoWessie=c("black",
           rep("#027C1E", 2),
           # "#E4764E",
           # "#00BC7F",
           rep("#C03180", 3))


ggplot(ma, aes(x = Genes, y = Bins, fill = Copy_number)) + geom_tile(color="white") + # Plot heatmap
  coord_fixed( # If we want to have squares instead of rectangles, we should use coord_fixed()
    ratio=1/2) + # Rectangles
  #facet_grid(sample1~Bins) + # Having separated plots...maybe someday is useful
  
  scale_fill_gradientn (colours =
                          #paletteer_d("ggsci::cyan_material"),
                          #paletteer_c("grDevices::Spectral"),
                          #paletteer_d("ggsci::lime_material"),
                          #paletteer_d("RColorBrewer::BuGn"),
                          #paletteer_d("RColorBrewer::BuPu"), 
                          #paletteer_d("RColorBrewer::Oranges"),
                          #paletteer_d("RColorBrewer::GnBu"),
                          #paletteer_d("RColorBrewer::PuRd"),
                          #paletteer_d("RColorBrewer::RdPu"),
                          #paletteer_d("RColorBrewer::YlGnBu"),
                          #paletteer_d("RColorBrewer::PuBu"),
                          #paletteer_d("RColorBrewer::Greys"),
                          #paletteer_d("RColorBrewer::Greens"),
                          #paletteer_dynamic("cartography::sand.pal", 20),
                          #paletteer_dynamic("cartography::brown.pal", 20),
                        paletteer_dynamic("cartography::orange.pal", 20),
                        
                        values = scales::rescale(c(0,0.001,0.002,0.003,0.004,
                                                   0.2,0.4,0.6,0.8, 1,
                                                   1.2,1.4,1.6,1.8,2))) +

                        
                        
                        # values = scales::rescale(c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1,
                        #                            1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2))) +
  scale_y_discrete(limits = rev(levels(as.factor(ma$Bins)))) +  # Reverse axis. We want form A to Z. 
  # 
  # 
  xlab(label = "Biotin transporters") +
  ylab(label = NULL) +
  labs(fill="Abundance (copy number)") +
  # 
  theme(
    legend.position = "right",  # Maybe you want it left?
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.text = element_text(size = 10),
    legend.title = element_text(size=10, face="bold"),
    
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(face = 'bold', size = 15,
                                margin = margin(t = 0.5, # Top margin
                                                r = 0.2, # Right margin
                                                b = 1, # Bottom margin
                                                l = 0, # Left margin
                                                unit = "cm")),
    
    axis.text.x = element_text(face="bold",
                               color = NoWessie,
                               angle=45,
                               hjust=1,
                               size = 10),  # It's angled 45 degress. It can be set to '0'.
    axis.text.y = element_text(face="bold",
                               color = rev(Wessie), # Colouring our axis text depending on the sample to which our bins belong
                               angle=0,
                               hjust=1,
                               size =9),
    
    panel.background = element_rect(fill = "transparent", colour = NA),  # Transparent
    plot.background = element_rect(fill = "transparent", colour = NA),   # Transparent
    
  )

# sum(ma$Copy_number == 0) # To check if everything is well represented
# sum(ma$Copy_number != 0)

# new_bins$IC22_maxbin.010$bins$table # To see which taxa belongs our bin


## SAVE ##

ggsave("Biotin_transporters_Job1.png", width = 30, height = 25, units = "cm", bg = "transparent", dpi = 300)
ggsave("Biotin_transporters_Job1.pdf", width = 30, height = 25, units = "cm", bg = "transparent", dpi = 300)

