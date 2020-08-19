######################
##### HISTOGRAMS #####
######################

library(igraph)
library(network)
library(ergm)
library(RColorBrewer)
library(intergraph)

# --------------------
setwd("/home/noe/Università/in_corso/SocialNetworkAnalysis/")

data = read.table("network_data/skill.txt")
edgelist = as.matrix(data[, c(1,2)])
skill_w = graph_from_edgelist(edgelist, directed = T)

weight = data[,3]
sum(weight == 0)     #check for missing values (no one)
bin.weight = as.numeric(weight>3) #binary weights
skill = delete_edges(skill_w, E(skill_w)[bin.weight == 0])

V(skill)$name = seq(1,77)
V(skill)$location = scan("network_data/location.txt")
V(skill)$org_level = scan("network_data/org_level.txt")
V(skill)$tenure = scan("network_data/tenure.txt")

# let us extract the adjacency matrix
Y = get.adjacency(skill, sparse = F)
diag(Y) = NA

# number of nodes
n = nrow(Y)


## observed statistics:
# ------------------------------------------------------------------------------------
# DENSITY
density_obs = graph.density(skill)

# RECIPROPCITY
reciprocity_obs = reciprocity(skill)

# TRANSITIVITY
transitivity_obs = transitivity(as.undirected(skill, mode = 'collapse'))

# ASSORTATIVE MIXING
ass_loc_obs = assortativity(skill, V(skill)$location)
ass_orglev_obs = assortativity(skill, V(skill)$org_level)
ass_ten_obs = assortativity(skill, V(skill)$tenure)

# IN/OUT DEGREE
ci_indeg_obs = centr_degree(skill, loops = F, mode = 'in')$centralization
ci_outdeg_obs = centr_degree(skill, loops = F, mode = 'out')$centralization

# IN/OUT CLOSENESS
ci_inclo_obs = centr_clo(skill, mode = 'in')$centralization
ci_outclo_obs = centr_clo(skill, mode = 'out')$centralization

# BETWEENESS CENTRALITY
ci_betw_obs= centr_betw(skill, directed = T)$centralization

# EIGENVECTOR CENTRALITY
eig_obs = mean(eigen_centrality(skill, directed = T, scale = F)$vector)

mod_loc_obs = assortativity(skill, V(skill)$location)
mod_lev_obs = assortativity(skill, V(skill)$org_level)
mod_ten_obs = assortativity(skill, V(skill)$tenure)

######################################
Y = get.adjacency(skill, sparse = F)
diag(Y) = NA

# number of nodes
n = nrow(Y)
net = network(Y, directed = T)
net

# add attributes
net %v% 'location' = scan("network_data/location.txt")
net %v% 'org_level' = scan("network_data/org_level.txt")
net %v% 'tenure' = scan("network_data/tenure.txt")

##################
##### MODELS #####
##################

### SRG MODEL ###
mod00 = ergm(net ~ edges)
summary(mod00)
set.seed(1)
p.MLE = mean(Y, na.rm = T)
B = 1000

dens.sim = trans.sim = reci.sim = modu.sim = ass_loc.sim = ass_orglev.sim = ass_ten.sim = c()
ci_indeg.sim = ci_outdeg.sim = ci_inclo.sim = ci_outclo.sim = ci_betw.sim = eig.sim = c()

for(b in 1:B){
  Y.sim = matrix(NA, n, n) 
  tmp = rbinom(n^2,1,p.MLE ) 
  diag(Y.sim) = NA
  Y.sim = matrix(tmp, n,n)
  g.sim = graph_from_adjacency_matrix(Y.sim)
  
  dens.sim[b] = graph.density(g.sim)
  reci.sim[b] = reciprocity(g.sim)
  trans.sim[b] = transitivity(as.undirected(g.sim, mode = 'collapse'))
  
  ass_loc.sim[b] = assortativity(g.sim, V(skill)$location)
  ass_orglev.sim[b] = assortativity(g.sim, V(skill)$org_level)
  ass_ten.sim[b] = assortativity(g.sim, V(skill)$tenure)
  
  ci_indeg.sim[b] = centr_degree(g.sim, loops = F, mode = 'in')$centralization
  ci_outdeg.sim[b] = centr_degree(g.sim, loops = F, mode = 'out')$centralization
  
  ci_inclo.sim[b] = centr_clo(g.sim, mode = 'in')$centralization
  ci_outclo.sim[b] = centr_clo(g.sim, mode = 'out')$centralization
  
  ci_betw.sim[b] = centr_betw(g.sim, directed = T)$centralization
  
  eig.sim[b] = mean(eigen_centrality(g.sim, directed = T, scale = F)$vector)
}
stats_mod00 = cbind(dens.sim, reci.sim, trans.sim, ass_loc.sim, ass_orglev.sim, ass_ten.sim, 
                   ci_indeg.sim, ci_outdeg.sim, ci_inclo.sim, ci_outclo.sim, ci_betw.sim, eig.sim)

### NH-SRG MODEL ###

y = c(Y)
rowIdx = row(Y)
colIdx = col(Y)

rowidx = c(rowIdx)
colidx = c(colIdx)
mod1 = glm(y ~ factor(rowidx) + factor(colidx), family = "binomial")
mod1
summary(mod1)
n = nrow(Y)
n
B = 1000
pij = mod1$fitted.values
dens.sim = trans.sim = reci.sim = modu.sim = ass_loc.sim = ass_orglev.sim = ass_ten.sim = c()
ci_indeg.sim = ci_outdeg.sim = ci_inclo.sim = ci_outclo.sim = ci_betw.sim = eig.sim = c()

for(b in 1:B){
  tmp = rbinom(n*(n-1),1,pij)
  Y.sim = matrix(, n,n); Y.sim[row(Y.sim) != col(Y.sim)] = tmp
  g.sim = graph_from_adjacency_matrix(Y.sim)
  
  dens.sim[b] = graph.density(g.sim)
  reci.sim[b] = reciprocity(g.sim)
  trans.sim[b] = transitivity(as.undirected(g.sim, mode = 'collapse'))
  
  ass_loc.sim[b] = assortativity(g.sim, V(skill)$location)
  ass_orglev.sim[b] = assortativity(g.sim, V(skill)$org_level)
  ass_ten.sim[b] = assortativity(g.sim, V(skill)$tenure)
  
  ci_indeg.sim[b] = centr_degree(g.sim, loops = F, mode = 'in')$centralization
  ci_outdeg.sim[b] = centr_degree(g.sim, loops = F, mode = 'out')$centralization
  
  ci_inclo.sim[b] = centr_clo(g.sim, mode = 'in')$centralization
  ci_outclo.sim[b] = centr_clo(g.sim, mode = 'out')$centralization
  
  ci_betw.sim[b] = centr_betw(g.sim, directed = T)$centralization
  
  eig.sim[b] = mean(eigen_centrality(g.sim, directed = T, scale = F)$vector)
  
}
stats_mod1 = cbind(dens.sim, reci.sim, trans.sim, ass_loc.sim, ass_orglev.sim, ass_ten.sim, 
                   ci_indeg.sim, ci_outdeg.sim, ci_inclo.sim, ci_outclo.sim, ci_betw.sim, eig.sim)

### P1 MODEL ###
#mod2 = ergm(net ~ edges+sender+receiver+mutual, control=control.ergm(seed=1))
#summary(mod2)

#sim = simulate(mod2, burnin = 1000, nsim = 1000, verbose = TRUE, seed = 1)
#stats_mod2= as.matrix(t(sapply(sim, function(xx){fnc(xx)})))
load("results/P1_sender_reciver.RData")
summary(mod0)

fnc = function(xx){
  ig = asIgraph(xx)
  dens.sim = graph.density(ig)
  reci.sim = reciprocity(ig)
  trans.sim = transitivity(as.undirected(ig, mode = 'collapse'))
  
  ass_loc.sim = assortativity(ig, V(skill)$location)
  ass_orglev.sim = assortativity(ig, V(skill)$org_level)
  ass_ten.sim = assortativity(ig, V(skill)$tenure)
  
  ci_indeg.sim = centr_degree(ig, loops = F, mode = 'in')$centralization
  ci_outdeg.sim = centr_degree(ig, loops = F, mode = 'out')$centralization
  
  ci_inclo.sim = centr_clo(ig, mode = 'in')$centralization
  ci_outclo.sim = centr_clo(ig, mode = 'out')$centralization
  
  ci_betw.sim = centr_betw(ig, directed = T)$centralization
  eig.sim = mean(eigen_centrality(ig, directed = T, scale = F)$vector)
  
  return(cbind(dens.sim, reci.sim, trans.sim, ass_loc.sim, ass_orglev.sim, ass_ten.sim,
               ci_indeg.sim, ci_outdeg.sim, ci_inclo.sim, ci_outclo.sim,
               ci_betw.sim, eig.sim))
}


sim = simulate(mod0, burnin = 1000, nsim = 1000, verbose = TRUE, seed = 1)
stats_mod2 = as.matrix(t(sapply(sim, function(xx){fnc(xx)})))

dim(stats_mod2)

### COVARIATES MODEL ###

load("results/Markov_agg.RData")

summary(mod4a)
B = 1000

library(intergraph)

sim = simulate(mod4a, burnin = 1000, nsim = 1000, verbose = TRUE, seed = 1)

stats_mod4a = as.matrix(t(sapply(sim, function(xx){fnc(xx)})))

dim(stats_mod4a)

###############

load("results/social_agg.RData")

summary(mod5a)
B = 1000

library(intergraph)

sim = simulate(mod5a, burnin = 1000, nsim = 1000, verbose = TRUE, seed = 1)

stats_mod5a = as.matrix(t(sapply(sim, function(xx){fnc(xx)})))
dim(stats_mod5a)


library(ggplot2)
library(ggpubr)
library(reshape2)
library(gridExtra) #grid.arrange(plot1, plot2, ncol=2)

##################
### HISTOGRAMS ###
##################
#cols = brewer.pal(5, name = 'Pastel2')
a = stats_mod00[,1]; b = stats_mod1[,1]; c = stats_mod2[,1]; d = stats_mod4a[,1]; dd = stats_mod5a[,1]
e = data.frame(mod_0 = a, mod_1 = b, mod_2 = c, mod_3 = d, mod_4 = dd)
stats = melt(e)

g1 = ggplot(stats, aes(x=value, fill=variable)) + 
  geom_density(alpha=0.25) + 
  geom_vline(aes(xintercept = density_obs, color = "red"), linetype = "dashed") +
  scale_colour_discrete(name="Data", breaks= "red", labels= "Observed\nStatistic") +
  labs(title="Density", x="value", y = "density") +
  theme_classic() +
  scale_fill_discrete(name="Model\nSimulations",
                      breaks=c("mod_0", "mod_1", "mod_2", "mod_3","mod_4"),
                      labels=c("BRG", "NHSRG", "P1", "Markov", "Social circuit"))  
g1


a = stats_mod00[,2]; b = stats_mod1[,2]; c = stats_mod2[,2]; d = stats_mod4a[,2]; dd = stats_mod5a[,2]
e = data.frame(mod_0 = a, mod_1 = b, mod_2 = c, mod_3 = d, mod_4 = dd)
stats = melt(e)

g2 = ggplot(stats,aes(x=value, fill=variable)) + 
  geom_density(alpha=0.25) + 
  #scale_color_brewer(palette="Pastel2") +
  geom_vline(aes(xintercept = reciprocity_obs, color = "red"), linetype = "dashed") +
  scale_colour_discrete(name="Data", breaks= "red", labels= "Observed\nReciprocity") +
  labs(title="Reciprocity", x="value", y = "density") +
  theme_classic() +
  scale_fill_discrete(name="Model\nSimulations",
                      breaks=c("mod_0", "mod_1", "mod_2", "mod_3", "mod_4"),
                      labels=c("BRG", "NHSRG", "P1", "Markov", "Social circuit"))   
g2


a = stats_mod00[,3]; b = stats_mod1[,3]; c = stats_mod2[,3]; d = stats_mod4a[,3]; dd = stats_mod5a[,3]
e = data.frame(mod_0 = a, mod_1 = b, mod_2 = c, mod_3 = d, mod_4 = dd)
stats = melt(e)

g3 = ggplot(stats,aes(x=value, fill=variable)) + 
  geom_density(alpha=0.25) + 
  #scale_color_brewer(palette="Pastel2") +
  geom_vline(aes(xintercept = transitivity_obs, color = "red"), linetype = "dashed") +
  scale_colour_discrete(name="Data", breaks= "red", labels= "Observed\nTransitivity") +
  labs(title="Transitivity", x="value", y = "density") +
  theme_classic() +
  scale_fill_discrete(name="Model\nSimulations",
                      breaks=c("mod_0", "mod_1", "mod_2", "mod_3", "mod_4"),
                      labels=c("BRG", "NHSRG", "P1", "Markov", "Social circuit"))  
g3


a = stats_mod00[,4]; b = stats_mod1[,4]; c = stats_mod2[,4]; d = stats_mod4a[,4]; dd = stats_mod5a[,4]
e = data.frame(mod_0 = a, mod_1 = b, mod_2 = c, mod_3 = d, mod_4 = dd)
stats = melt(e)

g4 = ggplot(stats,aes(x=value, fill=variable)) + 
  geom_density(alpha=0.25) + 
  #scale_color_brewer(palette="Pastel2") +
  geom_vline(aes(xintercept = ass_loc_obs, color = "red"), linetype = "dashed") +
  scale_colour_discrete(name="Data", breaks= "red", labels= "Observed\nLocation\nModularity") +
  labs(title="Location Modularity", x="value", y = "density") +
  theme_classic() +
  scale_fill_discrete(name="Model\nSimulations",
                      breaks=c("mod_0", "mod_1", "mod_2", "mod_3", "mod_4"),
                      labels=c("BRG", "NHSRG", "P1", "Markov", "Social circuit"))  
g4



a = stats_mod00[,5]; b = stats_mod1[,5]; c = stats_mod2[,5]; d = stats_mod4a[,5]; dd = stats_mod5a[,5]
e = data.frame(mod_0 = a, mod_1 = b, mod_2 = c, mod_3 = d, mod_4 = dd)
stats = melt(e)

g5 = ggplot(stats,aes(x=value, fill=variable)) + 
  geom_density(alpha=0.25) + 
  #scale_color_brewer(palette="Pastel2") +
  geom_vline(aes(xintercept = ass_orglev_obs, color = "red"), linetype = "dashed") +
  scale_colour_discrete(name="Data", breaks= "red", labels= "Observed\nOrg. Level\n Modularity") +
  labs(title="Organization Level Modularity", x="value", y = "density") +
  theme_classic() +
  scale_fill_discrete(name="Model\nSimulations",
                      breaks=c("mod_0", "mod_1", "mod_2", "mod_3", "mod_4"),
                      labels=c("BRG", "NHSRG", "P1", "Markov", "Social circuit"))
g5


a = stats_mod00[,6]; b = stats_mod1[,6]; c = stats_mod2[,6]; d = stats_mod4a[,6]; dd = stats_mod5a[,6]
e = data.frame(mod_0 = a, mod_1 = b, mod_2 = c, mod_3 = d, mod_4 = dd)
stats = melt(e)

g6 = ggplot(stats,aes(x=value, fill=variable)) + 
  geom_density(alpha=0.25) + 
  #scale_color_brewer(palette="Pastel2") +
  geom_vline(aes(xintercept = ass_ten_obs, color = "red"), linetype = "dashed") +
  scale_colour_discrete(name="Data", breaks= "red", labels= "Observed\nTenure\n Modularity") +
  labs(title="Tenure Modularity", x="value", y = "density") +
  theme_classic() +
  scale_fill_discrete(name="Model\nSimulations",
                      breaks=c("mod_0", "mod_1", "mod_2", "mod_3","mod_4"),
                      labels=c("BRG", "NHSRG", "P1", "Markov", "Social circuit"))
g6



a = stats_mod00[,7]; b = stats_mod1[,7]; c = stats_mod2[,7]; d = stats_mod4a[,7]; dd = stats_mod5a[,7]
e = data.frame(mod_0 = a, mod_1 = b, mod_2 = c, mod_3 = d, mod_4 = dd)
stats = melt(e)

g7 = ggplot(stats,aes(x=value, fill=variable)) + 
  geom_density(alpha=0.25) + 
  #scale_color_brewer(palette="Pastel2") +
  geom_vline(aes(xintercept = ci_indeg_obs, color = "red"), linetype = "dashed") +
  scale_colour_discrete(name="Data", breaks= "red", labels= "Observed\nIn Degree") +
  labs(title="In Degree", x="value", y = "density") +
  theme_classic() +
  scale_fill_discrete(name="Model\nSimulations",
                      breaks=c("mod_0", "mod_1", "mod_2", "mod_3", "mod_4"),
                      labels=c("BRG", "NHSRG", "P1", "Markov", "Social circuit"))
g7

a = stats_mod0[,8]; b = stats_mod1[,8]; c = stats_mod2[,8]; d = stats_mod4a[,8]; dd = stats_mod5a[,8]
e = data.frame(mod_0 = a, mod_1 = b, mod_2 = c, mod_3 = d, mod_4 = dd)
stats = melt(e)

g8 = ggplot(stats,aes(x=value, fill=variable)) + 
  geom_density(alpha=0.25) + 
  #scale_color_brewer(palette="Pastel2") +
  geom_vline(aes(xintercept = ci_outdeg_obs, color = "red"), linetype = "dashed") +
  scale_colour_discrete(name="Data", breaks= "red", labels= "Observed\nOut Degree") +
  labs(title="Out Degree", x="value", y = "density") +
  theme_classic() +
  scale_fill_discrete(name="Model\nSimulations",
                      breaks=c("mod_0", "mod_1", "mod_2", "mod_3", "mod_4"),
                      labels=c("BRG", "NHSRG", "P1", "Markov", "Social circuit"))
g8




a = stats_mod0[,9]; b = stats_mod1[,9]; c = stats_mod2[,9]; d = stats_mod4a[,9]; d = stats_mod5a[,9]
e = data.frame(mod_0 = a, mod_1 = b, mod_2 = c, mod_3 = d, mod_4 = dd)
stats = melt(e)

g9 = ggplot(stats,aes(x=value, fill=variable)) + 
  geom_density(alpha=0.25) + 
  #scale_color_brewer(palette="Pastel2") +
  geom_vline(aes(xintercept = ci_inclo_obs, color = "red"), linetype = "dashed") +
  scale_colour_discrete(name="Data", breaks= "red", labels= "Observed\nIn \nCloseness") +
  labs(title="In Closeness", x="value", y = "density") +
  theme_classic() +
  scale_fill_discrete(name="Model\nSimulations",
                      breaks=c("mod_0", "mod_1", "mod_2", "mod_3", "mod_4"),
                      labels=c("BRG", "NHSRG", "P1", "Markov", "Social circuit"))
g9


a = stats_mod0[,10]; b = stats_mod1[,10]; c = stats_mod2[,10]; d = stats_mod4a[,10]; dd = stats_mod5a[,10]
e = data.frame(mod_0 = a, mod_1 = b, mod_2 = c, mod_3 = d, mod_4 = dd)
stats = melt(e)

g10 = ggplot(stats,aes(x=value, fill=variable)) + 
  geom_density(alpha=0.25) + 
  #scale_color_brewer(palette="Pastel2") +
  geom_vline(aes(xintercept = ci_outclo_obs, color = "red"), linetype = "dashed") +
  scale_colour_discrete(name="Data", breaks= "red", labels= "Observed\nOut \nCloseness") +
  labs(title="Out Closeness", x="value", y = "density") +
  theme_classic() +
  scale_fill_discrete(name="Model\nSimulations",
                      breaks=c("mod_0", "mod_1", "mod_2", "mod_3", "mod_4"),
                      labels=c("BRG", "NHSRG", "P1", "Markov", "Social circuit"))
g10

a = stats_mod0[,11]; b = stats_mod1[,11]; c = stats_mod2[,11]; d = stats_mod4a[,11]; dd = stats_mod5a[,11]
e = data.frame(mod_0 = a, mod_1 = b, mod_2 = c, mod_3 = d, mod_4 = dd)
stats = melt(e)

g11 = ggplot(stats,aes(x=value, fill=variable)) + 
  geom_density(alpha=0.25) + 
  #scale_color_brewer(palette="Pastel2") +
  geom_vline(aes(xintercept = ci_betw_obs, color = "red"), linetype = "dashed") +
  scale_colour_discrete(name="Data", breaks= "red", labels= "Observed\nBetweenness") +
  labs(title="Betweenness", x="value", y = "density") +
  theme_classic() +
  scale_fill_discrete(name="Model\nSimulations",
                      breaks=c("mod_0", "mod_1", "mod_2", "mod_3", "mod_4"),
                      labels=c("BRG", "NHSRG", "P1", "Markov", "Social circuit"))
g11


a = stats_mod0[,12]; b = stats_mod1[,12]; c = stats_mod2[,12]; d = stats_mod4a[,12];dd = stats_mod5a[,12]
e = data.frame(mod_0 = a, mod_1 = b, mod_2 = c, mod_3 = d, mod_4 = dd)
stats = melt(e)

g12 = ggplot(stats,aes(x=value, fill=variable)) + 
  geom_density(alpha=0.25) + 
  #scale_color_brewer(palette="Pastel2") +
  geom_vline(aes(xintercept = eig_obs, color = "red"), linetype = "dashed") +
  scale_colour_discrete(name="Data", breaks= "red", labels= "Observed\nEigenvector") +
  labs(title="Eigenvector", x="value", y = "density") +
  theme_classic() +
  scale_fill_discrete(name="Model\nSimulations",
                      breaks=c("mod_0", "mod_1", "mod_2", "mod_3", "mod_4"),
                      labels=c("BRG", "NHSRG", "P1", "Markov", "Social circuit"))
g12

pdf(file = "model_simulations.pdf", height = 5)
ggarrange(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12, nrow = 4, ncol = 3, common.legend = T, legend="right")
dev.off()
