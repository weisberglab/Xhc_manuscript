###########################
#  Xanthomonous Phylowave #
###########################

# Loading necessary functions
source(file = 'data/phylowave/2_1_Index_computation_20240909.R')
source(file = 'data/phylowave/2_2_Lineage_detection_20240909.R')
source(file = 'data/phylowave/2_3_Lineage_fitness_20240909.R')
source(file = 'data/phylowave/2_4_Lineage_defining_mutations.R')

# Load necessary packages
library(data.table)
library(ape, quiet = T)
library(phytools, quiet = T)
library(stringr, quiet = T)
library(MetBrewer, quiet = T)
library(parallel, quiet = T) 
library(mgcv, quiet = T)
library(cowplot, quiet = T)
library(ggplot2, quiet = T)
library(ggtree, quiet = T)
library(cmdstanr, quiet = T) 
library(binom, quiet = T)
library(phytools, quiet = T)
library(vcfR, quiet = T)
library(tidyr)
library(tidyverse)
library(dplyr)
library(purrr)
library(ggplot2)

#############
# Load data #
#############

tree_ani1 = read.nexus('data/phylowave/run2_mcc.beast.tre')

#write.tree(tree_ani1, file = 'data/phylowave/time_tree.mcc.nwk')

#tree_ani1 <- drop.tip(tree_ani1, c("17CO03101", "17CO07001"))

## Make sure the tree is binary, and ladderized
tree_ani1 = collapse.singles(ladderize(multi2di(tree_ani1, random = F), right = F))

## Subset meta and reorder based on strain present iin the tree_ani
meta <- read.table("data/phylowave/Xanthomonas_phylowave.tsv", sep="\t", header=TRUE, quote="")

## Reorder meta based on 
tree_tip_labels <- tree_ani1$tip.label

meta_subset <- meta %>%
  filter(strain %in% tree_tip_labels)

meta_reordered <- as.data.frame(meta_subset) %>%
  dplyr::slice(match(tree_tip_labels, strain))

dim(meta_reordered)

## Names all sequences
names_seqs = meta_reordered$strain
(n_seq = length(names_seqs))

## Collection times of all sequences
times_seqs = meta_reordered$year
min(times_seqs)
max(times_seqs)

## Clade assignment
clades_seqs = sapply(names_seqs, function(x)tail(str_split(x, pattern = '/')[[1]],1))


######################
# Detecting lineages #
######################

# Get cophenetic distance matrix from the tree
dist_mat <- cophenetic(tree_ani1)

# Hierarchical clustering
hc <- hclust(as.dist(dist_mat))

# Choose number of clusters (lineages)
k <- 2
clusters <- cutree(hc, k = k)

# Assign strain to lineage
strain_lineage_map <- setNames(paste0("lineage", clusters), names(clusters))

length(strain_lineage_map)

##########################
# Compute index dynamics #
##########################

## Index parameters
## Length genome 
genome_length = 5000000
## Mutation rate 
mutation_rate = 1.38e-7 # substitution per site per year
## Parameters for the index
timescale = 100

## Window of time on which to search for samples in the population
wind = 5 #years

## Compute pairwise distance matrix
genetic_distance_mat = dist.nodes.with.names(tree_ani1)

# Get the time of each internal root
nroot = length(tree_ani1$tip.label) + 1 ## Root number
distance_to_root = genetic_distance_mat[nroot,]
root_height = times_seqs[which(names_seqs == names(distance_to_root[1]))] - distance_to_root[1]
nodes_height = root_height + distance_to_root[n_seq+(1:(n_seq-1))]


# Preperation of data tips and nodes 
## Meta-data with all nodes 
dataset_with_nodes = data.frame('ID' = c(1:n_seq, n_seq+(1:(n_seq-1))),
                                'name_seq' = c(names_seqs, n_seq+(1:(n_seq-1))),
                                'time' = c(times_seqs, nodes_height),
                                'is.node' = c(rep('no', n_seq), rep('yes', (n_seq-1))), 
                                'clade' = c(strain_lineage_map, rep(NA, n_seq-1)))

# Compute index of every tips and nodes
dataset_with_nodes$index = compute.index(time_distance_mat = genetic_distance_mat, 
                                         timed_tree = tree_ani1, 
                                         time_window = wind,
                                         metadata = dataset_with_nodes, 
                                         mutation_rate = mutation_rate,
                                         timescale = timescale,
                                         genome_length = genome_length)


#Plot tree & index below, with colors from NextStrain clades
## Color key for Nextstrain clades
colors_clade = met.brewer(name="Signac", n=length(levels(as.factor(dataset_with_nodes$clade))), type="continuous")

## Color of each node, based on the key
dataset_with_nodes$clade_color = as.factor(dataset_with_nodes$clade)
clade_labels = levels(dataset_with_nodes$clade_color)
levels(dataset_with_nodes$clade_color) = colors_clade
dataset_with_nodes$clade_color = as.character(dataset_with_nodes$clade_color)

# Plot tree and index
par(mfrow = c(2,1), oma = c(0,0,0,0), mar = c(4,4,0,0))

#min_year = -18000
min_year = 1930
max_year = 2023

## Tree
plot(tree_ani1, show.tip.label = F, 
     edge.color = 'grey0', edge.width = 0.25,
     x.lim = c(min_year, max_year)-root_height)

tiplabels(pch = 16, col = dataset_with_nodes$clade_color, cex = 0.3)
axisPhylo_NL(side = 1, root.time = root_height, backward = F,
             at_axis = seq(min_year, max_year, 0.5)-root_height,
             lab_axis = seq(min_year, max_year, 0.5), lwd = 0.5)



## Index
plot(dataset_with_nodes$time, 
     dataset_with_nodes$index, 
     col = adjustcolor(dataset_with_nodes$clade_color, alpha.f = 1),
     bty = 'n', xlim = c(min_year, max_year), cex = 0.4,
     pch = 16, bty = 'n', ylim = c(0, 1), 
     main = paste0(''), 
     ylab = 'Index', xlab = 'Time (years)', xaxt = 'n', yaxt = 'n')

axis(2, las = 2, lwd = 0.5)
axis(1, lwd = 0.5)

# Color key
legend('topright', 
       legend = clade_labels,
       fill = colors_clade, border = colors_clade,
       cex = 0.5, bty = 'n', ncol = 5)


#####################
# Lineage detection #
#####################

## Setting parameters
time_window_initial = 2030;
time_window_increment = 100;
p_value_smooth = 0.05
weight_by_time = 0.1
k_smooth = -1
plot_screening = F
min_descendants_per_tested_node = 5
min_group_size = 5
weighting_transformation = c('inv_sqrt')

parallelize_code = T
number_cores = 4

max_stepwise_deviance_explained_threshold = 0
max_groups_found = 14
stepwise_AIC_threshold = 0

keep_track = T


## Run the detection function
start_time = Sys.time()

potential_splits = find.groups.by.index.dynamics(timed_tree = tree_ani1,
                                                 metadata = dataset_with_nodes,
                                                 node_support = tree_ani1$edge.length[match((n_seq+1):(2*n_seq-1), tree_ani1$edge[,2])],
                                                 threshold_node_support = 1/(29903*0.00081),
                                                 time_window_initial = time_window_initial,
                                                 time_window_increment = time_window_increment,
                                                 min_descendants_per_tested_node = min_descendants_per_tested_node,
                                                 min_group_size = min_group_size,
                                                 p_value_smooth = p_value_smooth,
                                                 stepwise_deviance_explained_threshold = max_stepwise_deviance_explained_threshold,
                                                 stepwise_AIC_threshold = stepwise_AIC_threshold,
                                                 weight_by_time = weight_by_time,
                                                 weighting_transformation = weighting_transformation,
                                                 k_smooth = k_smooth,
                                                 parallelize_code = parallelize_code,
                                                 number_cores = number_cores, 
                                                 plot_screening = plot_screening,
                                                 max_groups_found = max_groups_found, 
                                                 keep_track = keep_track)
end_time = Sys.time()
print(end_time - start_time)


df_explained_dev = data.frame('N_groups' = 0:length(potential_splits$best_dev_explained),
                              'Non_explained_deviance' = (1-c(potential_splits$first_dev, potential_splits$best_dev_explained)),
                              'Non_explained_deviance_log' = log(1-c(potential_splits$first_dev, potential_splits$best_dev_explained)))
df_explained_dev$Non_explained_deviance_log = df_explained_dev$Non_explained_deviance_log-min(df_explained_dev$Non_explained_deviance_log)

par(mfrow = c(1,2), oma = c(2,2,1,1), mar = c(2,2,2,0.5), mgp = c(0.75,0.25,0), cex.axis=0.5, cex.lab=0.5, cex.main=0.7, cex.sub=0.5)
plot(df_explained_dev$N_groups,
     df_explained_dev$Non_explained_deviance,
     bty = 'n', ylim = c(0, ceiling(10*max(df_explained_dev$Non_explained_deviance))/10),
     xaxt = 'n', yaxt = 'n', pch = 16, main = 'linear scale', cex = 0.5, 
     ylab = 'Non-explained deviance (%)', xlab = 'Number of groups')
axis(1, lwd = 0.5, tck=-0.02)
axis(2, las = 2, at = seq(0,ceiling(10*max(df_explained_dev$Non_explained_deviance))/10,0.1),
     labels = seq(0, ceiling(10*max(df_explained_dev$Non_explained_deviance))/10,0.1)*100, lwd = 0.5, tck=-0.02)

plot(df_explained_dev$N_groups,
     (df_explained_dev$Non_explained_deviance),
     log = 'y',
     ylim = c(0.01, 1),
     bty = 'n',
     xaxt = 'n', yaxt = 'n', pch = 16, main = 'log scale', cex = 0.5, 
     ylab = 'Non-explained deviance (%) - log scale', xlab = 'Number of groups')
axis(1, lwd = 0.5, tck=-0.02)
axis(2, las = 2, at = c(0.01, 0.1, 0.25, 0.5, 1),
     labels = c(0.01, 0.1, 0.25, 0.5, 1)*100, lwd = 0.5, tck=-0.02)


split = merge.groups(timed_tree = tree_ani1, metadata = dataset_with_nodes, 
                     initial_splits = potential_splits$potential_splits, 
                     group_count_threshold = 30, group_freq_threshold = 0.01)


## Label sequences with new groups
dataset_with_nodes$groups = as.factor(split$groups)
## Reorder labels by time of emergence
name_groups = levels(dataset_with_nodes$groups)
time_groups_world = NULL
for(i in 1:length(name_groups)){
  time_groups_world = c(time_groups_world, min(dataset_with_nodes$time[which(dataset_with_nodes$groups == name_groups[i] &
                                                                               dataset_with_nodes$is.node == 'no')]))
}

levels(dataset_with_nodes$groups) = match(name_groups, order(time_groups_world, decreasing = T))
dataset_with_nodes$groups = as.numeric(as.character(dataset_with_nodes$groups))
dataset_with_nodes$groups = as.factor(dataset_with_nodes$groups)

## Update names in split list
split$tip_and_nodes_groups = match(split$tip_and_nodes_groups, order(time_groups_world, decreasing = T))
names(split$tip_and_nodes_groups) = 1:length(split$tip_and_nodes_groups)
split$groups = as.factor(split$groups)
levels(split$groups) = match(name_groups, order(time_groups_world, decreasing = T))
split$groups = as.numeric(as.character(split$groups))

## Choose color palette
n_groups <- length(name_groups)
colors_groups = met.brewer(name="Cross", n=n_groups, type="continuous")



## Color each group
#dataset_with_nodes$group_color = dataset_with_nodes$groups
#levels(dataset_with_nodes$group_color) = colors_groups
#dataset_with_nodes$group_color = as.character(dataset_with_nodes$group_color)



lineage_labels <- c("1"="G", "2"="D", "3"="F",
                    "4"="E", "5"="B", "6"="A", "7"="C")

lineage_colors <- c("G"="#D85952", "D"="#EB7D41",  "F"="#FFBB44", "E"="#79987B", "B"="#206575", "A"="#122451", "C"="#C969A1")
                     

dataset_with_nodes$lineage <- lineage_labels[as.character(dataset_with_nodes$groups)]
dataset_with_nodes$lineage_colors <- lineage_colors[as.character(dataset_with_nodes$lineage)]



## Index colored by group
par(mfrow = c(2,1), oma = c(0,0,0,0), mar = c(4,4,0,0))

min_year = 1600 #1930

## Tree
plot(tree_ani1, show.tip.label = FALSE, 
     edge.color = 'grey', edge.width = 0.25,
     x.lim = c(min_year, max_year)-root_height)

tiplabels(pch = 16, col = dataset_with_nodes$lineage_colors, cex = 0.3)

#tiplabels(text = dataset_with_nodes2$source, col = dataset_with_nodes2$group_color, cex = 0.3, offset=1, frame = "none", bg = "transparent")
axisPhylo_NL(side = 1, root.time = root_height, backward = F,
             at_axis = seq(min_year, max_year, 0.5)-root_height,
             lab_axis = seq(min_year, max_year, 0.5), lwd = 0.5)


plot(dataset_with_nodes$time, 
     dataset_with_nodes$index, 
     col = dataset_with_nodes$lineage_colors,
     #col = adjustcolor(dataset_with_nodes$lineage_colors, alpha.f = 1),
     bty = 'n', xlim = c(min_year, max_year), cex = 0.5,
     pch = 16, bty = 'n', ylim = c(0, 1), 
     main = paste0(''), #log = 'y',
     ylab = 'Index', xlab = 'Time (years)', 
     #ylab = '', xlab = '', 
     yaxt = 'n')

axis(2, las=1, hadj = 2)
#axis(2, las = 2)


sorted_colors <- lineage_colors[sort(names(lineage_colors))]

# Color key
legend('topright', 
       #legend = name_groups,
       legend = names(sorted_colors),
       #fill = colors_groups, border = colors_groups,
       fill = sorted_colors, border = sorted_colors,
       cex = 0.5, bty = 'n', ncol = 5)



## Color sampling efforts
#dataset_with_nodes2 <- merge(dataset_with_nodes, meta[,c('strain', 'source')], by.x='name_seq', by.y='strain', all.x=T)

#dataset_with_nodes2$group_color = dataset_with_nodes2$groups
#levels(dataset_with_nodes2$group_color) = colors_groups
#dataset_with_nodes2$group_color = as.character(dataset_with_nodes2$group_color)

ggtree(tree_ani1, mrsd="2021-12-01", size=0.1)  %<+% dataset_with_nodes[,-1] +
  theme_tree2() +
  geom_tippoint(aes(color=as.character(lineage)), size=2, offset=-30) +
  scale_color_manual(values = lineage_colors) +
  coord_cartesian(xlim = c(1930, 2022), clip = "off") 



## Tree with index-defined groups
groups = matrix(dataset_with_nodes$groups[which(dataset_with_nodes$is.node == 'no')], ncol = 1)
colnames(groups) = 'groups'
rownames(groups) = dataset_with_nodes$name_seq[which(dataset_with_nodes$is.node == 'no')]
cols = as.character(colors_groups)
names(cols) = as.character(1:max(as.numeric(name_groups)))


write.csv(dataset_with_nodes, file = 'data/phylowave/lineages.csv', quote = F, row.names=F)

dataset_with_nodes <- read.csv('data/phylowave/lineages.csv')
traits <- dataset_with_nodes[which(dataset_with_nodes$is.node == 'no'), c('name_seq', 'lineage')]
traits$linE <- 0
traits$linF <- 0
traits$linG <- 0
traits$linFG <- 0
traits$linEFG <- 0
traits$linDEFG <- 0

traits[which(traits$lineage %in% c('E')), ]$linE <- 1
traits[which(traits$lineage %in% c('G')), ]$linG <- 1
traits[which(traits$lineage %in% c('F')), ]$linF <- 1

traits[which(traits$lineage %in% c('F','G')), ]$linFG <- 1
traits[which(traits$lineage %in% c('E','F','G')), ]$linEFG <- 1
traits[which(traits$lineage %in% c('D','E','F','G')), ]$linDEFG <- 1

#write.csv(traits[,c('name_seq', 'linE', 'linF', 'linG', 'linFG', 'linEFG', 'linDEFG')], 'data/phylowave/traits.txt', quote = F, row.names = F)

####################################
# Quantify the fitness of lineages #
####################################

#dataset_with_nodes$groups <- dataset_with_nodes$clade

start_time = Sys.time()
## Load and compile stan code (this can take a few minutes)
model_compiled <- cmdstan_model(stan_file = '2_Functions/Model_multinomial_logistic_birthdeath_lineage_fitness_20231220.stan')

## Run model on BV3
res_fitness = estimate_rel_fitness_groups_with_branches(dataset_with_nodes = dataset_with_nodes,
                                                        tree = tree_ani1,
                                                        min_year = min(times_seqs), 
                                                        window = 5,
                                                        model_compiled = model_compiled,
                                                        iter_warmup = 250, iter_sampling = 500, refresh = 50, seed = 1)
end_time = Sys.time()
print(end_time - start_time)


#dataset_with_nodes$groups <- dataset_with_nodes$clade

order_colors = order(as.numeric(split$tip_and_nodes_groups))
colour_lineage = colors_groups[match(split$tip_and_nodes_groups[order_colors], name_groups)]

plot_fit_data_new(data = res_fitness$data,
                  Chains = res_fitness$chains,
                  colour_lineage = colour_lineage,
                  xmin = 1930, xmax = 2023.5)

plot_estimated_fitness_ref_ancestral(data = res_fitness$data,
                                     Chains = res_fitness$chains,
                                     colour_lineage = colour_lineage, 
                                     gentime = 1)


legend('topright', 
       legend = name_groups,
       fill = colour_lineage, border = lineage_colors, #colour_lineage,
       cex = 0.5, bty = 'n', ncol = 5)

# Get the data behind the plot
Chains = res_fitness$chains
gentime = 1
(betas = apply((Chains$beta*gentime), MARGIN = 2, function(x) mean.and.ci(x)))


# Plot relative fitness data
relative_fitness <- data.frame(lineages=c('G', 'F', 'E', 'A', 'B', 'D', 'C'),
           mean=c(0.1407866, 0.09151681, 0.02769257, -0.001666361, -0.00307039, -0.015791068, NA), 
           lower=c(0.1220803, 0.07983693, 0.01878275, -0.008193635, -0.008363283, -0.030674307, NA), 
           upper=c(0.1577785, 0.10472012, 0.03715668, 0.004632169, 0.001847676, -0.000825646, NA))



ggplot(relative_fitness, aes(x=lineages, y=mean, fill=lineages)) +
  geom_bar(stat='identity') +
  scale_fill_manual(values=lineage_colors, name='Lineage') +
  xlab('Lineage') + ylab('Relative Fitness') +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, position=position_dodge(.9)) +
  theme_classic()

ggsave('results/Main_figures/Lineage_fitness.pdf', width = 5, heigh = 4)
