library(igraph)
library(network)
library(ergm)
library(RColorBrewer)

# --------------------
# Load the data
# --------------------
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
# ------------------------------------------------------------------------------------
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

####################
# MODEL ESTIMATION #
####################
mod0 = ergm(net ~ edges)
mod0
summary(mod0)

mod2 = ergm(net ~ edges+sender+receiver+mutual, control=control.ergm(seed=1))
summary(mod2)

mod3 = ergm(net ~ edges+sender+mutual, control=control.ergm(seed=1))
summary(mod3)

mod4 = ergm(net ~ edges+receiver+mutual, control=control.ergm(seed=1))
summary(mod4)

mod5 = ergm(net ~ edges+mutual, control=control.ergm(seed=1))
summary(mod5)

c(BIC(mod0), BIC(mod2), BIC(mod3), BIC(mod4), BIC(mod5))

# check for convergence
dev.new()
mcmc.diagnostics(mod4)

graph.density(skill)
# or
mean(Y, na.rm = T)
# or, via the expit function
exp(mod5$coef)
exp(mod5$coef)/(1+exp(mod5$coef))

#############################
### ASSESS GOD OF FITNESS ###
#############################
library(intergraph)


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


statistics = as.matrix(t(sapply(sim, function(xx){fnc(xx)})))
statistics

par(mfrow = c(3,3))

dens.sim = statistics[,1]
low = pmin(min(dens.sim), density_obs) - 0.05
up = pmax(max(dens.sim), density_obs) + 0.05
hist(dens.sim, col = "lightgray", xlim = c(low, up), xlab = NA, 
     main = c("Density", paste('pvalue',mean(dens.sim >= density_obs))))
segments(x0=density_obs, x1=density_obs, y0=-2, y1=322, col='red', lwd = 2)
box()

reci.sim = statistics[,2]
low = pmin(min(reci.sim), reciprocity_obs) - 0.05
up = pmax(max(reci.sim), reciprocity_obs) + 0.05
hist(reci.sim, col = "lightgray", xlim = c(low, up), xlab = NA,
     main = c("Reciprocity", paste('pvalue', mean(reci.sim >= reciprocity_obs))))
segments(x0=reciprocity_obs, x1=reciprocity_obs, y0=-2, y1=1002, col='red', lwd = 2)
box()

trans.sim = statistics[,3]
low = pmin(min(trans.sim), transitivity_obs) - 0.05
up = pmax(max(trans.sim), transitivity_obs) + 0.05
hist(trans.sim, col = "lightgray", xlim = c(low, up), xlab = NA,
     main = c('Transitivity', paste('pvalue', mean(trans.sim >= transitivity_obs))))
segments(x0=transitivity_obs, x1=transitivity_obs, y0=-2, y1=902, col='red', lwd = 2)
box()

ci_indeg.sim = statistics[,4]
low = pmin(min(ci_indeg.sim), ci_indeg_obs) - 0.05
up = pmax(max(ci_indeg.sim), ci_indeg_obs) + 0.05
hist(ci_indeg.sim, col = "lightgray",xlim = c(low, up),xlab = NA,
     main = c('In-Degree', paste('pvalue',mean(ci_indeg.sim >= ci_indeg_obs))))
segments(x0=ci_indeg_obs, x1=ci_indeg_obs, y0=-2, y1=2352, col='red', lwd = 2)
box()

ci_outdeg.sim = statistics[,5]
low = pmin(min(ci_outdeg.sim), ci_outdeg_obs) - 0.05
up = pmax(max(ci_outdeg.sim), ci_outdeg_obs) + 0.05
hist(ci_outdeg.sim, col = "lightgray", xlim = c(low, up),xlab = NA,
     main = c("Out-Degree", paste('pvalue', mean(ci_outdeg.sim >= ci_outdeg_obs))))
segments(x0=ci_outdeg_obs, x1=ci_outdeg_obs, y0=-2, y1=2332, col='red', lwd = 2)
box()

ci_inclo.sim = statistics[,6]
low = pmin(min(ci_inclo.sim), ci_inclo_obs) - 0.05
up = pmax(max(ci_inclo.sim), ci_inclo_obs) + 0.05
hist(ci_inclo.sim, col = "lightgray", xlim = c(low, up),xlab = NA,
     main = c("In-Closeness", paste('pvalue', mean(ci_inclo.sim >= ci_inclo_obs))))
segments(x0=ci_inclo_obs, x1=ci_inclo_obs, y0=-2, y1=1322, col='red', lwd = 2)
box()

ci_outclo.sim = statistics[,7]
low = pmin(min(ci_outclo.sim), ci_outclo_obs) - 0.05
up = pmax(max(ci_outclo.sim), ci_outclo_obs) + 0.05
hist(ci_outclo.sim, col = "lightgray", xlim = c(low, up), xlab = NA,
     main = c("Out-Closeness", paste('pvalue', mean(ci_outclo.sim >= ci_outclo_obs))))
segments(x0=ci_outclo_obs, x1=ci_outclo_obs, y0=-2, y1=1222, col='red', lwd = 2)
box()

ci_betw.sim = statistics[,8]
low = pmin(min(ci_betw.sim), ci_betw_obs) - 0.005
up = pmax(max(ci_betw.sim), ci_betw_obs) + 0.005
hist(ci_betw.sim, col = "lightgray", xlim = c(low, up), xlab = NA,
     main = c("Betweenness", paste('pvalue', mean(ci_betw.sim >= ci_betw_obs))))
segments(x0=ci_betw_obs, x1=ci_betw_obs, y0=-2, y1=1502, col='red', lwd = 2)
box()

eig.sim = statistics[,9]
low = pmin(min(eig.sim), eig_obs) - 0.005
up = pmax(max(eig.sim), eig_obs) + 0.005
hist(eig.sim, col = "lightgray",xlim = c(low, up), xlab = NA,
     main = c("Eigenvector Centrality", paste('pvalue', mean(eig.sim >= eig_obs))))
segments(x0=eig_obs, x1=eig_obs, y0=-2, y1=802, col='red', lwd = 2)
box()

# compute an approximate p-value
p_den = mean(dens.sim >= density_obs)
p_tr = mean(trans.sim >= transitivity_obs)
p_rec = mean(reci.sim >= reciprocity_obs)
p_ideg = mean(ci_indeg.sim >= ci_indeg_obs)
p_odeg = mean(ci_outdeg.sim >= ci_outdeg_obs)
p_iclo = mean(ci_inclo.sim >= ci_inclo_obs)
p_oclo = mean(ci_outclo.sim >= ci_outclo_obs)
p_bet = mean(ci_betw.sim >= ci_betw_obs)
p_eig = mean(eig.sim >= eig_obs)

table = as.matrix(round(rbind(density= c(density_obs, mean(dens.sim), p_den), 
                              transitivity= c(transitivity_obs, mean(trans.sim), p_tr),
                              reciprocity = c(reciprocity_obs, mean(reci.sim), p_rec),
                              in_degree = c(ci_indeg_obs, mean(ci_indeg.sim), p_ideg),
                              out_degree = c(ci_outdeg_obs, mean(ci_outdeg.sim), p_odeg),
                              in_closeness = c(ci_inclo_obs, mean(ci_inclo.sim), p_iclo),
                              out_closeness = c(ci_outclo_obs, mean(ci_outclo.sim), p_oclo),
                              betweeness = c(ci_betw_obs, mean(ci_betw.sim), p_bet),
                              eigenvector = c(eig_obs, mean(eig.sim), p_eig)), 4))
library(xtable)
to_latex = xtable(table)
to_latex


