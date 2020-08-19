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

mod0 = ergm(net ~ edges+sender+receiver+mutual, control=control.ergm(seed=1))
#save(mod0, file = './results/P1_sender_reciver.RData')

load('./results/P1_sender_reciver.RData') #è salvato come mod0
summary(mod0)
mod0$mle.lik

exp(mod0$loglikelihood)

mod = ergm(net ~ edges , control=control.ergm(seed=1))
summary(mod)
mod$mle.lik
mod1 = ergm(net ~ edges+mutual+ nodefactor('location', base= c(1,3,4))+nodefactor("org_level", base=c(3,4))+
              nodefactor("tenure", base=c(2,3,4))+nodematch('location'), control=control.ergm(seed=1))
#save(mod1, file = './results/P1_attr_no_send_rec.RData') 
load("results/P1_attr_no_send_rec.RData")  # è salvato come mod1
summary(mod1)
mod1$mle.lik

dev.new()
pdf(file = "./diagnostics_modp1_attr.pdf")
mcmc.diagnostics(mod1)
dev.off()


  loc.new = rep(0, n)
  loc.new[V(skill)$location == 2] = 1
  loc.new[V(skill)$location == 4] = 2
  net %v% 'loc.new' = loc.new

mod2 = ergm(net ~ edges+mutual+ nodefactor('loc.new')+nodefactor("org_level")+
              nodefactor("tenure")+nodematch('loc.new'), control=control.ergm(seed=1))


#save(mod2, file = './results/P1_newattr_no_send_rec.RData')
load("results/P1_newattr_no_send_rec.RData") # è salvato come mod2
summary(mod2)

mod4 = ergm(net ~ edges + sender + receiver + mutual + 
              nodefactor("loc.new") + nodefactor("org_level") + nodefactor("tenure") +
              nodematch("loc.new"),
            control = control.ergm(seed = 1))
#save(mod4, file = './results/P1_newattr_send_rec.RData')
load("results/P1_newattr_no_send_rec.RData") # è salvato come mod4
summary(mod4)


mod5 = ergm(net ~ edges + mutual + istar(2) + istar(3) + ostar(2) + ostar(3) + triangle + 
              nodefactor("loc.new") + nodefactor("org_level") + nodefactor("tenure") +
              nodematch("loc.new")+ gwesp(0.25, fixed = T),
            control = control.ergm(seed = 1, MCMC.samplesize = 4096))
summary(mod5)


library(intergraph)
sim = simulate(mod6, burnin = 1000, nsim = 1000, verbose = TRUE, seed = 1)

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
  mod_loc.sim = assortativity(ig, V(skill)$location)
  mod_lev.sim = assortativity(ig, V(skill)$org_level)
  mod_ten.sim = assortativity(ig, V(skill)$tenure)
  
  return(cbind(dens.sim, reci.sim, trans.sim, ci_indeg.sim, ci_outdeg.sim, ci_inclo.sim, ci_outclo.sim,
               ci_betw.sim, eig.sim,mod_loc.sim, mod_lev.sim, mod_ten.sim))
}


statistics = as.matrix(t(sapply(sim, function(xx){fnc(xx)})))
statistics

par(mfrow = c(3,4))

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
mod_loc.sim = statistics[,10]
low = pmin(min(mod_loc.sim), mod_loc_obs) - 0.005
up = pmax(max(mod_loc.sim), mod_loc_obs) + 0.005
hist(mod_loc.sim, col = "lightgray",xlim = c(low, up), xlab = NA,
     main = c("Modularity - location", paste('pvalue', mean(mod_loc.sim >= mod_loc_obs))))
segments(x0=mod_loc_obs, x1=mod_loc_obs, y0=-2, y1=802, col='red', lwd = 2)
box()

mod_lev.sim = statistics[,11]
low = pmin(min(mod_lev.sim), mod_lev_obs) - 0.005
up = pmax(max(mod_lev.sim), mod_lev_obs) + 0.005
hist(mod_lev.sim, col = "lightgray",xlim = c(low, up), xlab = NA,
     main = c("Modularity - organization level", paste('pvalue', mean(mod_lev.sim >= mod_lev_obs))))
segments(x0=mod_lev_obs, x1=mod_lev_obs, y0=-2, y1=802, col='red', lwd = 2)
box()

mod_ten.sim = statistics[,12]
low = pmin(min(mod_ten.sim), mod_ten_obs) - 0.005
up = pmax(max(mod_ten.sim), mod_ten_obs) + 0.005
hist(mod_ten.sim, col = "lightgray",xlim = c(low, up), xlab = NA,
     main = c("Modularity - tenure", paste('pvalue', mean(mod_ten.sim >= mod_ten_obs))))
segments(x0=mod_ten_obs, x1=mod_ten_obs, y0=-2, y1=802, col='red', lwd = 2)
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
