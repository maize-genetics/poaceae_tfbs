# Charlie Hale, 4/2/2023
# Plot Andropogoneae environment on phylogeny

library(dplyr)
library(ape)
library(phytools)
library(geiger)
# Plots data on phylogeny
# This function requires a vector of continuous trait values, a vector of taxa names/IDs corresponding
# to the trait values (matching order), and a phylogenetic tree read in from ape (read.tree).
# Tip labels on the tree will be joined with vector of taxa names/IDs, but taxa that
# are present in only the tree or or only the taxa/trait will be filtered out.
# Default is to output a circular phylogeny, otherwise specify "phylogram".
# Phylogeny is colored by trait values, reconstructing ancestral trait values

# Usage example: plotTraitOnTree(myTrait, myTaxa, myTree, 
#                visible_tip_labels = T, legend.title = "MinTemp",
#                colors = c("blue", "white","red"), shape = "fan")

plotTraitOnTree <- function(trait, taxa, tree, visible_tip_labels = FALSE, 
                            plot.dim = c(20, 50), plot.title = NA, legend.title = "",
			                      colors = c("blue", "white","red"), legend = NULL, shape = "fan") {
  # Drop trait observations for taxa that are not present in tree
  taxaNotInTree <- taxa[!taxa %in% tree$tip.label]
  
  trait.df <- data.frame(cbind(trait, taxa)) %>%
    filter(!taxa %in% taxaNotInTree) %>%
    filter(!is.na(trait)) 
  trait.vec <- as.numeric(as.matrix(trait.df)[,1]) # Convert to vector to enable plotting
  names(trait.vec) <- trait.df$taxa
  # Drop tips that are not present in trait/taxa list
  tipsToDrop <- tree$tip.label[!tree$tip.label %in% trait.df$taxa]
  tree.pruned <- drop.tip(tree, tipsToDrop)
  # Set negative branch lengths to 0
  tree.pruned$edge.length[tree.pruned$edge.length < 0] <- 0
  # Check that names in tree match names in data matrix
  name.check(tree.pruned, trait.vec)
  # Determine whether tip labels are visible based on user arguments
  font.size <- c(0.00001,1)
  if(visible_tip_labels) {
    font.size <- c(0.4,1)
  }
  # Generate plot
  options(repr.plot.width=plot.dim[1], repr.plot.height=plot.dim[2])
  plot <- contMap(tree.pruned,trait.vec, plot = FALSE)
  plot <-setMap(plot, colors)
  plot(plot, fsize = font.size,outline=F, type = shape, title = plot.title,
       leg.txt=legend.title, legend = legend)
  }




