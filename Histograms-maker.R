######################
##### HISTOGRAMS #####
######################

library(igraph)
library(network)
library(ergm)
library(RColorBrewer)

setwd("/home/noe/Università/in_corso/SocialNetworkAnalysis/")

data = read.table("network_data/skill.txt")
edgelist = as.matrix(data[, c(1,2)])
skill_w = graph_from_edgelist(edgelist, directed = T)

weight = data[,3]
sum(weight == 0)     #check for missing values (no one)
bin.weight = as.numeric(weight>3) #binary weights
skill = delete_edges(skill_w, E(skill_w)[bin.weight == 0])

# let us extract the adjacency matrix
Y = get.adjacency(skill, sparse = F)
diag(Y) = NA

# number of nodes
n = nrow(Y)


## observed statistics:
# DENSITY
density_obs = graph.density(skill)

# RECIPROPCITY
reciprocity_obs = reciprocity(skill)

# TRANSITIVITY
transitivity_obs = transitivity(as.undirected(skill, mode = 'collapse'))

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
mod0 = ergm(net ~ edges)
summary(mod0)
set.seed(1)

B = 5000
m =  sum(Y, na.rm = TRUE)

dens.sim = trans.sim = reci.sim = modu.sim = ass_loc.sim = ass_orglev.sim = ass_ten.sim = c()
ci_indeg.sim = ci_outdeg.sim = ci_inclo.sim = ci_outclo.sim = ci_betw.sim = eig.sim = c()
mod_loc.sim = mod_lev.sim = mod_ten.sim = c()

for(b in 1:B){
  Y.sim = matrix(NA, n, n) 
  ones = rep(1, m)
  zeros = rep(0, n*(n-1) - m)
  all = c(ones, zeros)
  Y.sim[col(Y.sim) != row(Y.sim)] = sample(all, n*(n-1))
  g.sim = graph_from_adjacency_matrix(Y.sim)
  
  dens.sim[b] = graph.density(g.sim)
  reci.sim[b] = reciprocity(g.sim)
  trans.sim[b] = transitivity(as.undirected(g.sim, mode = 'collapse'))
  
  ci_indeg.sim[b] = centr_degree(g.sim, loops = F, mode = 'in')$centralization
  ci_outdeg.sim[b] = centr_degree(g.sim, loops = F, mode = 'out')$centralization
  
  ci_inclo.sim[b] = centr_clo(g.sim, mode = 'in')$centralization
  ci_outclo.sim[b] = centr_clo(g.sim, mode = 'out')$centralization
  
  ci_betw.sim[b] = centr_betw(g.sim, directed = T)$centralization
  
  eig.sim[b] = mean(eigen_centrality(g.sim, directed = T, scale = F)$vector)
  mod_loc.sim[b] = assortativity(ig, V(skill)$location)
  mod_lev.sim[b] = assortativity(ig, V(skill)$org_level)
  mod_ten.sim[b] = assortativity(ig, V(skill)$tenure)
}
stats_mod0 = cbind(dens.sim, reci.sim, trans.sim, ci_indeg.sim, ci_outdeg.sim, 
                   ci_inclo.sim, ci_outclo.sim, ci_betw.sim, eig.sim)


### NH-SRG MODEL ###

mod1 = ergm(net~edges + sender + receiver)
summary(mod1)

m = ecount(skill)
B = 1000

dens.sim = trans.sim = reci.sim = modu.sim = ass_loc.sim = ass_orglev.sim = ass_ten.sim = c()
ci_indeg.sim = ci_outdeg.sim = ci_inclo.sim = ci_outclo.sim = ci_betw.sim = eig.sim = c()

for(b in 1:B){
  ones = rep(1, m)
  zeros = rep(0, n*(n-1) - m)
  all = c(ones, zeros)
  tmp = sample(all, n*(n-1))
  Y.sim = matrix(tmp, n,n)
  diag(Y.sim) = NA
  g.sim = graph_from_adjacency_matrix(Y.sim)
  
  dens.sim[b] = graph.density(g.sim)
  reci.sim[b] = reciprocity(g.sim)
  trans.sim[b] = transitivity(as.undirected(g.sim, mode = 'collapse'))
  
  ci_indeg.sim[b] = centr_degree(g.sim, loops = F, mode = 'in')$centralization
  ci_outdeg.sim[b] = centr_degree(g.sim, loops = F, mode = 'out')$centralization
  
  ci_inclo.sim[b] = centr_clo(g.sim, mode = 'in')$centralization
  ci_outclo.sim[b] = centr_clo(g.sim, mode = 'out')$centralization
  
  ci_betw.sim[b] = centr_betw(g.sim, directed = T)$centralization
  
  eig.sim[b] = mean(eigen_centrality(g.sim, directed = T, scale = F)$vector)
  
}
stats_mod1 = cbind(dens.sim, reci.sim, trans.sim, ci_indeg.sim, ci_outdeg.sim, 
                   ci_inclo.sim, ci_outclo.sim, ci_betw.sim, eig.sim)



### P1 MODEL ###

.f = function(){
mod2 = ergm(net ~ edges+sender+receiver+mutual, control=control.ergm(seed=1))
load("results/mod0.RData")
summary(mod0)
save(mod0, )

sim = simulate(mod2, burnin = 1000, nsim = 1000, verbose = TRUE, seed = 1)

fnc = function(xx){
  ig = asIgraph(xx)
  dens.sim = graph.density(ig)
  reci.sim = reciprocity(ig)
  trans.sim = transitivity(as.undirected(ig, mode = 'collapse'))
  
  ci_indeg.sim = centr_degree(ig, loops = F, mode = 'in')$centralization
  ci_outdeg.sim = centr_degree(ig, loops = F, mode = 'out')$centralization
  
  ci_inclo.sim = centr_clo(ig, mode = 'in')$centralization
  ci_outclo.sim = centr_clo(ig, mode = 'out')$centralization
  
  ci_betw.sim = centr_betw(ig, directed = T)$centralization
  eig.sim = mean(eigen_centrality(ig, directed = T, scale = F)$vector)
  
  return(cbind(dens.sim, reci.sim, trans.sim, ci_indeg.sim, ci_outdeg.sim, ci_inclo.sim, ci_outclo.sim,
               ci_betw.sim, eig.sim))
}


stats_mod2= as.matrix(t(sapply(sim, function(xx){fnc(xx)})))
stats_mod2

write.table(stats_mod2, file = "stats_mod2.txt", sep = "\t", )
}

stats_mod2 = read.table("stats_mod2.txt")

##################
### HISTOGRAMS ###
##################

cols = brewer.pal(5, name = 'Pastel2')

par(mfrow = c(3,4))

low = pmin(min(stats_mod0[,1]),min(stats_mod1[,1]), min(stats_mod2[,1]), density_obs) - 0.05
up = pmax(max(stats_mod0[,1]),max(stats_mod1[,1]), max(stats_mod2[,1]), density_obs) + 0.05
hist(stats_mod2[,1], col = cols[3], xlim = c(low, up), xlab = NA, main = "Density")
hist(stats_mod1[,1], col = cols[2], add = T)
#hist(stats_mod2[,1], col = cols[2], add = T)
segments(x0=density_obs, x1=density_obs, y0=-2, y1=322, col='red', lwd = 2)
box()

low = pmin(min(stats_mod0[,2]),min(stats_mod1[,2]), min(stats_mod2[,2]), density_obs) - 0.05
up = pmax(max(stats_mod0[,2]),max(stats_mod1[,2]), max(stats_mod2[,2]), density_obs) + 0.05
hist(stats_mod0[,2], col = cols[1], xlim = c(low, up), xlab = NA, main = "Reciprocity")
hist(stats_mod1[,2], col = cols[2], add = T)
hist(stats_mod2[,2], col = cols[3], add = T)
segments(x0=reciprocity_obs, x1=reciprocity_obs, y0=-2, y1=1502, col='red', lwd = 2)
box()

low = pmin(min(stats_mod0[,3]),min(stats_mod1[,3]), min(stats_mod2[,3]), reciprocity_obs) - 0.05
up = pmax(max(stats_mod0[,3]),max(stats_mod1[,3]), max(stats_mod2[,3]), reciprocity_obs) + 0.05
hist(stats_mod0[,3], col = cols[1], xlim = c(low, up), xlab = NA, main = "Transitivity")
hist(stats_mod1[,3], col = cols[2], add = T)
hist(stats_mod2[,3], col = cols[3], add = T)
segments(x0=reciprocity_obs, x1=reciprocity_obs, y0=-2, y1=1502, col='red', lwd = 2)
box()

low = pmin(min(stats_mod0[,4]),min(stats_mod1[,4]), min(stats_mod2[,4]), ci_indeg_obs) - 0.05
up = pmax(max(stats_mod0[,4]),max(stats_mod1[,4]), max(stats_mod2[,4]), ci_indeg_obs) + 0.05
hist(stats_mod0[,4], col = cols[1], xlim = c(low, up), xlab = NA, main = "In-Degree")
hist(stats_mod1[,4], col = cols[2], add = T)
hist(stats_mod2[,4], col = cols[3], add = T)
segments(x0=ci_indeg_obs, x1=ci_indeg_obs, y0=-2, y1=1502, col='red', lwd = 2)
box()

low = pmin(min(stats_mod0[,5]),min(stats_mod1[,5]), min(stats_mod2[,5]), ci_outdeg_obs) - 0.05
up = pmax(max(stats_mod0[,5]),max(stats_mod1[,5]), max(stats_mod2[,5]), ci_outdeg_obs) + 0.05
hist(stats_mod0[,5], col = cols[1], xlim = c(low, up), xlab = NA, main = "Out-Degree")
hist(stats_mod1[,5], col = cols[2], add = T)
hist(stats_mod2[,5], col = cols[3], add = T)
segments(x0=ci_outdeg_obs, x1=ci_outdeg_obs, y0=-2, y1=2502, col='red', lwd = 2)
box()

low = pmin(min(stats_mod0[,6]),min(stats_mod1[,6]), min(stats_mod2[,6]), ci_inclo_obs) - 0.05
up = pmax(max(stats_mod0[,6]),max(stats_mod1[,6]), max(stats_mod2[,6]), ci_inclo_obs) + 0.05
hist(stats_mod0[,6], col = cols[1], xlim = c(low, up), xlab = NA, main = "In-Closeness")
hist(stats_mod1[,6], col = cols[2], add = T)
hist(stats_mod2[,6], col = cols[3], add = T)
segments(x0=ci_inclo_obs, x1=ci_inclo_obs, y0=-2, y1=2502, col='red', lwd = 2)
box()

low = pmin(min(stats_mod0[,7]),min(stats_mod1[,7]), min(stats_mod2[,7]), ci_outclo_obs) - 0.05
up = pmax(max(stats_mod0[,7]),max(stats_mod1[,7]), max(stats_mod2[,7]), ci_outclo_obs) + 0.05
hist(stats_mod0[,7], col = cols[1], xlim = c(low, up), xlab = NA, main = "Out-Closeness")
hist(stats_mod1[,7], col = cols[2], add = T)
hist(stats_mod2[,7], col = cols[3], add = T)
segments(x0=ci_outclo_obs, x1=ci_outclo_obs, y0=-2, y1=2502, col='red', lwd = 2)
box()

low = pmin(min(stats_mod0[,8]),min(stats_mod1[,8]), min(stats_mod2[,8]), ci_betw_obs) - 0.005
up = pmax(max(stats_mod0[,8]),max(stats_mod1[,8]), max(stats_mod2[,8]), ci_betw_obs) + 0.005
hist(stats_mod0[,8], col = cols[1], xlim = c(low, up), xlab = NA, main = "Betweeness")
hist(stats_mod1[,8], col = cols[2], add = T)
hist(stats_mod2[,8], col = cols[3], add = T)
segments(x0=ci_betw_obs, x1=ci_betw_obs, y0=-2, y1=2502, col='red', lwd = 2)
box()

low = pmin(min(stats_mod0[,9]),min(stats_mod1[,9]), min(stats_mod2[,9]), eig_obs) - 0.001
up = pmax(max(stats_mod0[,9]),max(stats_mod1[,9]), max(stats_mod2[,9]), eig_obs) + 0.001
hist(stats_mod0[,9], col = cols[1], xlim = c(low, up), xlab = NA, main = "Eigenvector")
hist(stats_mod1[,9], col = cols[2], add = T)
hist(stats_mod2[,9], col = cols[3], add = T)
segments(x0=eig_obs, x1=eig_obs, y0=-2, y1=2502, col='red', lwd = 2)
box()

