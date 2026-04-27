#### Environmental Microbiology: Microbial Analysis of Zoo Feces
#### Author: Elle M Barnes
#### Created: 2025-09-29


##### Installing packages #####
    install.packages("ggplot2") # example of how to install package...
    install.packages("vegan")
    install.packages("dplyr")
    install.packages("phyloseq")
      # unique install for "phyloseq" package (if issue with above)
        if (!requireNamespace("BiocManager", quietly = TRUE))
          install.packages("BiocManager")
        BiocManager::install("phyloseq")
    
    # source('http://bioconductor.org/biocLite.R')
    # biocLite('phyloseq')
    
    install.packages("ape")
    install.packages("tidyr")

# unique install for "microbiome" package
    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("microbiome")


##### Setting your working directory #####

# If you did not save the data analysis folder to your desktop, you will
# need to change the path to match where you stored the downloaded folder

    # Run this command if you have Mac
      setwd("~/Desktop/EnviroMicro_ZooPoo_F25Analysis")
    
    # Run this command if you have Windows
      setwd("C:/Desktop/EnviroMicro_ZooPoo_F25Analysis")
  
# Load packages
    library(ggplot2)
    library(vegan)
    library(dplyr)
    library(phyloseq)
    library(ape)
    library(tidyr)
    library(microbiome)


##### Create physeq object #####

# Load data files: SV table (counts), taxonomy table, and metadata
    sv <- read.csv("featuretable16S.csv")
    row.names(sv) <- sv$ASV_ID
    sv$ASV_ID <- NULL
    
    taxa <- read.csv("taxonomy.csv")
    row.names(taxa) <- taxa$ASV_ID
    taxa$ASV_ID <- NULL
    
    metadata <- read.csv("zoo_metadata.csv")

# Convert dataframes into phyloseq format
    sv = otu_table(sv, taxa_are_rows = TRUE)
    
    taxa = tax_table(as.matrix(taxa))
    
    metadata <- sample_data(metadata)
    rownames(metadata) <- metadata$SampleID
    metadata$SampleID <- NULL

# merge SV & taxa data into phyloseq object
    physeq <- phyloseq(sv, taxa)

# create random phylo tree with ASV names as tips
    random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))

# create new phyloseq object with all info
    physeq16s <- phyloseq(sv, taxa, metadata, random_tree)
    
# Run the following line to view a summary of what your physeq object contains
    physeq16s # What does it tell you?


##### Filtering #####

# Filter out reads mapped to chloroplasts and mitochondria (host & 
    # diet contamination)

    physeq16s <- physeq16s %>%
      subset_taxa(
          Family  != "Mitochondria" &
          Order   != "Chloroplast")

# Create table, number of features for each phyla
    table(tax_table(physeq16s)[, "Phylum"], exclude = NULL)

# Remove features with ambiguous phylum annotation 
    physeq16s <- subset_taxa(physeq16s, !is.na(Phylum) &
                               !Phylum %in% c("", "Unknown", "uncharacterized",
                                              "unidentified")) 
    
    physeq16s
    
    # How many ASVs were removed (hint: compare the output of physeq16s with
    # it's previous output in your console)?

# Save your phyloseq object as an .RDS file so you don't have to run all these 
# commands next time you open this code. You will only have to run the first
# line of code in the next section: readRDS().
    saveRDS(physeq16s, "physeq16s_zoopoo.rds")



##### Loading in & Summarizing the data #####
    
## Important: I have included the code only for 2024, but you should repear these 
## analyses with "physeq16s_zoopoo_23.rds"
    
physeq16s <- readRDS("physeq16s_zoopoo_24.rds") 
    # this is a data object that contains all our data: read counts per sample,
    # taxonomy, and sample metadata

  # A) This produces an output for our sample metadata. What type of information
    # does it contain?
  sample_data(physeq16s)
  
  # B) Let's look at a summary of the number of reads we got back from 
    # sequencing...
  summarize_phyloseq(physeq16s)

# Let's ask R to produce some taxonomic summaries of our ASVs (i.e., unique 
  # 16S reads):

  # A) Number of ASVs for each domain. Which domains were identified by our 
    # V3-V4 sequencing?
  table(tax_table(physeq16s)[, "Domain"], exclude = NULL)

  # B) Number of ASVs for each phylum. Which is the most diverse phylum?
  table(tax_table(physeq16s)[, "Phylum"], exclude = NULL)


##### Alpha Diversity #####

# Let's explore alpha diversity using ASV richness (which R calls 'Observed') 
  # and Shannon Diversity:
  
  # Let's remove our controls for now...
  physeq16s.noctrl <- subset_samples(physeq16s, Replicate !="Control")
  
  adiv <- data.frame(
    "Observed" = phyloseq::estimate_richness(physeq16s.noctrl, measures = "Observed"),
    "Shannon" = phyloseq::estimate_richness(physeq16s.noctrl, measures = "Shannon"),
    "Species" = phyloseq::sample_data(physeq16s.noctrl)$Species,
    "Diet" = phyloseq::sample_data(physeq16s.noctrl)$Diet,
    "Individual" = phyloseq::sample_data(physeq16s.noctrl)$Individual,
    "Sex" = phyloseq::sample_data(physeq16s.noctrl)$Sex,
    "Replicate" = phyloseq::sample_data(physeq16s.noctrl)$Replicate,
    "DigestiveType" = phyloseq::sample_data(physeq16s.noctrl)$DigestiveType)
  
  # Create boxplots showing alpha diversity by Species
  adiv %>%
    gather(key = metric, value = value, c("Observed", "Shannon")) %>%
    mutate(metric = factor(metric, levels = c("Observed", "Shannon"))) %>%
    ggplot(aes(x = factor(Species, levels = c("Giraffe", "Zebra", "Rhino", 
                                       "Otter", "Tiger", "Polar Bear")),  # note that you will have to change these labels for 2023
               y = value, 
               fill= Diet)) +
    geom_boxplot() +
    geom_jitter(color="black", size = 1.5, alpha=0.5, position=position_jitter(0.2)) +
    labs(x = "", y = "") +
    facet_wrap(~ metric, scales = "free") +
    theme_bw()+
    scale_fill_manual(values = c("#466e95", "#b7d8ae"))
  
  # You can see that there is one Giraffe replicate that sequenced poorly.
  # Let's remove sample G_R3:
  physeq16s.noctrl.prune <- subset_samples(physeq16s.noctrl, sample_sums(physeq16s.noctrl) > 200)
  
  # Now, let's create a new data.frame and boxplot...
     adiv <- data.frame(
        "Observed" = phyloseq::estimate_richness(physeq16s.noctrl.prune, measures = "Observed"),
        "Shannon" = phyloseq::estimate_richness(physeq16s.noctrl.prune, measures = "Shannon"),
        "Species" = phyloseq::sample_data(physeq16s.noctrl.prune)$Species,
        "Diet" = phyloseq::sample_data(physeq16s.noctrl.prune)$Diet,
        "Individual" = phyloseq::sample_data(physeq16s.noctrl.prune)$Individual,
        "Sex" = phyloseq::sample_data(physeq16s.noctrl.prune)$Sex,
        "Replicate" = phyloseq::sample_data(physeq16s.noctrl.prune)$Replicate,
        "DigestiveType" = phyloseq::sample_data(physeq16s.noctrl.prune)$DigestiveType)
      
     # and we will name the plot: alpha
      alpha <- adiv %>%
        gather(key = metric, value = value, c("Observed", "Shannon")) %>%
        mutate(metric = factor(metric, levels = c("Observed", "Shannon"))) %>%
        ggplot(aes(x = factor(Species, levels = c("Giraffe", "Zebra", "Rhino",
                                                  "Otter", "Tiger", "Polar Bear")), # note that you will have to change these labels for 2023
                   y = value, 
                   fill= Diet)) +
        geom_boxplot() +
        geom_jitter(color="black", size = 1.5, alpha=0.5, position=position_jitter(0.2)) +
        labs(x = "", y = "") +
        facet_wrap(~ metric, scales = "free") +
        theme_bw()+
        scale_fill_manual(values = c("#466e95", "#b7d8ae")) 
  
      alpha # Better--let's use this one in your report!
      
      # let's save this figure to our folder (next 3 lines)
      png(width=8, height=3.5, units="in", res=600, file = "alphadiversity.png")
      print(alpha)
      dev.off()
  
  # But, are these differences statistically significant?
  
    # A) Richess
    summary(aov(Observed ~ adiv$Diet + adiv$Species, data=adiv))
    
    # B) Shannon diversity
    summary(aov(Shannon ~ adiv$Diet + adiv$Species, data=adiv))


  # ...and how does this compare to absolute abundance (# 16S reads by qPCR)  
    
    # I stored that information in our metadata file
      meta <- as.data.frame(as.matrix(metadata))
    
      summary(aov(Copy16S ~ Diet + DigestiveType, data=meta))
    
    
##### Beta Diversity: Composition #####
    
    # Microbial ecologists can compare composition by calculating the pairwise distance between samples.
    # Samples that are more similar in composition will have shorter pairwise distances (i.e., they will
    # be closer together). We can use many distance metrics to calculate compositional difference. Here,
    # we will use weighted Unifrac ("wunifrac") which takes into account both abundance and phylogenetic 
    # distances between ASVs in each sample.
    
    physeq16s.noctrl.log <- transform_sample_counts(physeq16s.noctrl.prune, function(x) log(1 + x))
    physeq16s.ord <- ordinate(physeq16s.noctrl.log, method = "MDS", distance = "wunifrac")
    
    
    # Visualize differences in composition using ordination
    
    comp <- plot_ordination(physeq16s.noctrl.log, physeq16s.ord, color = "Diet", shape = "Species")+ 
      geom_point(size = 4) +
      theme_bw()+
      scale_color_manual(values = c("#75639a", "#d698ae"))
    
    comp # view the plot
    
        # let's save this figure to our folder (next 3 lines)
        png(width=5, height=3, units="in", res=600, file = "betadiversity.png")
        print(comp)
        dev.off()
    
    # Is this difference statistically significant?
    
    wu.dist <- phyloseq::distance(physeq16s.noctrl.log, method="wunifrac") # calculates a pairwise distance for every point
    vegan::adonis2(wu.dist ~ sample_data(physeq16s.noctrl.prune)$Diet * sample_data(physeq16s.noctrl.prune)$Species,
                   by = "terms")   # this looks at the interaction between diet and species
    
    vegan::adonis2(wu.dist ~ sample_data(physeq16s.noctrl.prune)$Diet + sample_data(physeq16s.noctrl.prune)$Species,
                   by = "terms")   # this looks at the contribution of diet and species individually
    
    
    # This tells us if there are differences, but by collapsing each diverse community into a single
    # point, we lose the ability to tell which ASVs are different. We can look into these taxonomic 
    # differences by calculating and plotting ASV relative abundance.
    
      # A) Relative abundance of top phyla
    
        # Convert our reads into relative abundances:
        physeq16s.rel = transform_sample_counts(physeq16s.noctrl.prune, function(x) x/sum(x)*100)
        
        # Agglomerate taxa at phylum level (makes for nicer plots)
        glom <- tax_glom(physeq16s.rel, taxrank = 'Phylum', NArm = FALSE)
        ps.melt <- psmelt(glom)
        
        # rearrange the data and find the median
        ps.melt$Phylum <- as.character(ps.melt$Phylum)
        
        ps.melt <- ps.melt %>%
          group_by(Phylum) %>%
          mutate(median=median(Abundance))
        
        # Identify most abundant phyla (>5% of community)
        keep <- unique(ps.melt$Phylum[ps.melt$Abundance > 5]) 
        
        # Group all other phyla as "Other"
        ps.melt$Phylum[!(ps.melt$Phylum %in% keep)] <- "Other" 
        
        # put all this info together!
        ps.melt2 <- ps.melt %>%
          group_by(Sample, Diet, DigestiveType, Phylum) %>%
          dplyr::summarise(Abundance=sum(Abundance))
        
        # Create a color palette (we need at least 9 colors based on keep + "Other")
        palette <- c("#75639a", "#beb8c9", "#466e95", "#7396a9", "#85c0a3", "#b7d8ae", "#fae588", "#fcefb4","#Dee0e2") # feel free to change your color palette to any hex codes!
        
        
        # Plot relative abundance for each sample (the following chunk of code
        # is all one function)
        rel.abund <- ggplot(ps.melt2, aes(x = Sample, y = Abundance, 
                                                 fill = Phylum)) + 
                      geom_bar(stat = "identity", aes(fill=factor(Phylum, levels = c("Firmicutes",
                                                                                     "Bacteroidota", 
                                                                                     "Fusobacteriota",
                                                                                     "Proteobacteria", 
                                                                                     "Synergistota", 
                                                                                     "Verrucomicrobiota",
                                                                                     "Fibrobacterota", 
                                                                                     "Spirochaetota", 
                                                                                     "Other")))) + 
                      labs(x="", y="Relative Abundance (%)") +
                      theme_bw() + 
                      theme(legend.position = "right", 
                            legend.key.size = unit(0.2,'cm'),
                            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
                      scale_fill_manual(values=palette)+ 
                      facet_grid(.~factor(DigestiveType, levels = c('Ruminant',
                                                                    'Hindgut fermenter',
                                                                    'Monogastric carnivore')), 
                                 scales = "free", switch = "x", space = "free_x")+
                      theme(strip.placement = "outside")
      
        
        rel.abund # view the plot   
        
        png(width=7, height=3, units="in", res=600, file = "relativeabund.phylum.png")
        print(rel.abund)
        dev.off()
        
        
        
    # It looks like Firmicutes is an important phylum in the gut. Let's look at which orders dominate
    # our samples within this phylum. 
        
      # B) Relative abundance of top orders within Firmicutes
        
        # Subset our dataset to just Firmicutes:
        firmicutes <- subset_taxa(physeq16s.noctrl.prune, Phylum=="Firmicutes")
        firmicutes.noNA <- subset_taxa(firmicutes, !is.na(Order) & !Order %in% c("", "Unknown", "uncharacterized", "unidentified"))
        
        # Convert our reads into relative abundances:
        physeq16s.rel.firm = transform_sample_counts(firmicutes.noNA, function(x) x/sum(x)*100)
        
        # Agglomerate taxa at order level
        glom <- tax_glom(physeq16s.rel.firm, taxrank = 'Order', NArm = FALSE)
        ps.melt <- psmelt(glom)
        ps.melt$Order <- as.character(ps.melt$Order)
        
        ps.melt <- ps.melt %>%
          group_by(Order) %>%
          mutate(median=median(Abundance))
        
        # Identify most abundant orders (>10% of community--there are more
        # orders than phyla so we need to set a higher threshold)
        keep <- unique(ps.melt$Order[ps.melt$Abundance > 10]) 
        ps.melt$Order[!(ps.melt$Order %in% keep)] <- "Other" 
        
        ps.melt2 <- ps.melt %>%
          group_by(Sample, DigestiveType, Order) %>%
          dplyr::summarise(Abundance=sum(Abundance))
        
        # Create a NEW color palette (we need at least 9 colors based on keep + "Other")
        palette2 <- c("#565a94", "#b379b0", "#d698ae", "#e8b0b0", "#f0d4bd", "#f7f4df", "#c1e2eb", "#aac4e6","#Dee0e2") 
          # feel free to change your color palette to any hex codes!
        
        
        # Plot relative abundance for each sample
        firmi <- ggplot(ps.melt2, aes(x = Sample, y = Abundance, 
                             fill = Order)) + 
          geom_bar(stat = "identity", aes(fill=factor(Order, levels = c("Peptostreptococcales-Tissierellales",
                                                                         "Lachnospirales", 
                                                                         "Oscillospirales", 
                                                                         "Clostridiales",
                                                                         "Lactobacillales",
                                                                         "Bacillales", 
                                                                         "Christensenellales", 
                                                                         "Other")))) + 
          labs(x="", y="Relative Abundance (%)") +
          theme_bw() + 
          theme(legend.position = "right", 
                legend.key.size = unit(0.2,'cm'),
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
          scale_fill_manual(values=palette2)+ 
          facet_grid(.~factor(DigestiveType, levels = c('Ruminant',
                                                        'Hindgut fermenter',
                                                        'Monogastric carnivore')), 
                     scales = "free", switch = "x", space = "free_x")+
          theme(strip.placement = "outside")
    
        firmi # view the plot
        
        png(width=8, height=4, units="in", res=600, file = "firmicutes.png")
        print(firmi)
        dev.off()
        
        
        
        
    # These four plots should all be included in your lab report alongside a plot or table of the 16S gene 
    # copy abundances from our library prep lab (+ any associated statistical results). 
        
    # But, this is only the beginning! Can you think of another aspect of the data that you want to explore?
    # Ask Dr. Barnes to help you write the R code to answer that question.
    
