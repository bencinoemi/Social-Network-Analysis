# ------------------------------------
# ---  SRG Assessing significance  ---
# ------------------------------------

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

E = E(skill)
E.count = length(E)
V = V(skill)
n = length(V)

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

# MODULARITY
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

# Is the observed network coherent with the family of Binomial random graph models G(n, p)?
# a naive approach based on the best case scenario
# ------------------------------------------------------------------------------------

# maximum likelihood estimate of p
set.seed(1)
p.MLE = mean(Y, na.rm = T)

B = 1000
dens.sim = trans.sim = reci.sim = modu.sim = ass_loc.sim = ass_orglev.sim = ass_ten.sim = c()
ci_indeg.sim = ci_outdeg.sim = ci_inclo.sim = ci_outclo.sim = ci_betw.sim = eig.sim =  c()
for(b in 1:B){
  tmp = rbinom(n^2,1,p.MLE)  
  Y.sim = matrix(tmp, n,n)
  diag(Y.sim) = NA
  
  g.sim = graph_from_adjacency_matrix(Y.sim)
  print(is.connected(g.sim))
  
  dens.sim[b] = graph.density(g.sim)
  reci.sim[b] = reciprocity(g.sim)
  trans.sim[b] = transitivity(as.undirected(g.sim, mode = 'collapse'))
  #modu.sim[b] = modularity(g.sim)
  
  # ASSORTATIVE MIXING
  #ass_loc.sim[b] = assortativity(g.sim, V(g.sim)$location)
  #ass_orglev.sim[b] = assortativity(g.sim, V(g.sim)$org_level)
  #ass_ten.sim[b] = assortativity(g.sim, V(g.sim)$tenure)
  
  # IN/OUT DEGREE
  ci_indeg.sim[b] = centr_degree(g.sim, loops = F, mode = 'in')$centralization
  ci_outdeg.sim[b] = centr_degree(g.sim, loops = F, mode = 'out')$centralization
  
  # IN/OUT CLOSENESS
  ci_inclo.sim[b] = centr_clo(g.sim, mode = 'in')$centralization
  ci_outclo.sim[b] = centr_clo(g.sim, mode = 'out')$centralization
  
  # BETWEENESS CENTRALITY
  ci_betw.sim[b] = centr_betw(g.sim, directed = T)$centralization
  eig.sim[b] = mean(eigen_centrality(g.sim, directed = T, scale = F)$vector)
                  
}

# Graphical comparison
par(mfrow = c(1,1))

low = pmin(min(dens.sim), density_obs) - 0.05
up = pmax(max(dens.sim), density_obs) + 0.05
hist(dens.sim, col = "lightgray", main = paste('pvalue',mean(dens.sim >= density_obs)),
     xlab = "Density", xlim = c(low, up))
segments(x0=density_obs, x1=density_obs, y0=-2, y1=322, col='red', lwd = 2)

low = pmin(min(reci.sim), reciprocity_obs) - 0.05
up = pmax(max(reci.sim), reciprocity_obs) + 0.05
hist(reci.sim, col = "lightgray", main = paste('pvalue', mean(reci.sim >= reciprocity_obs)),
     xlab = "Reciprocity", xlim = c(low, up))
segments(x0=reciprocity_obs, x1=reciprocity_obs, y0=-2, y1=302, col='red', lwd = 2)


low = pmin(min(trans.sim), transitivity_obs) - 0.05
up = pmax(max(trans.sim), transitivity_obs) + 0.05
hist(trans.sim, col = "lightgray", main = paste('pvalue', mean(trans.sim >= transitivity_obs)),
     xlab = 'Transitivity', xlim = c(low, up))
segments(x0=transitivity_obs, x1=transitivity_obs, y0=-2, y1=202, col='red', lwd = 2)


low = pmin(min(ci_indeg.sim), ci_indeg_obs) - 0.05
up = pmax(max(ci_indeg.sim), ci_indeg_obs) + 0.05
hist(ci_indeg.sim, col = "lightgray", main = paste('pvalue',mean(ci_indeg.sim >= ci_indeg_obs)),
     xlab = 'In-Degree', xlim = c(low, up))
segments(x0=ci_indeg_obs, x1=ci_indeg_obs, y0=-2, y1=352, col='red', lwd = 2)

low = pmin(min(ci_outdeg.sim), ci_outdeg_obs) - 0.05
up = pmax(max(ci_outdeg.sim), ci_outdeg_obs) + 0.05
hist(ci_outdeg.sim, col = "lightgray", main = paste('pvalue', mean(ci_outdeg.sim >= ci_outdeg_obs)),
     xlab = "Out-Degree", xlim = c(low, up))
segments(x0=ci_outdeg_obs, x1=ci_outdeg_obs, y0=-2, y1=332, col='red', lwd = 2)


low = pmin(min(ci_inclo.sim), ci_inclo_obs) - 0.05
up = pmax(max(ci_inclo.sim), ci_inclo_obs) + 0.05
hist(ci_inclo.sim, col = "lightgray", main = paste('pvalue', mean(ci_inclo.sim >= ci_inclo_obs)),
     xlab = "In-Closeness", xlim = c(low, up))
segments(x0=ci_inclo_obs, x1=ci_inclo_obs, y0=-2, y1=422, col='red', lwd = 2)


low = pmin(min(ci_outclo.sim), ci_outclo_obs) - 0.05
up = pmax(max(ci_outclo.sim), ci_outclo_obs) + 0.05
hist(ci_outclo.sim, col = "lightgray", main = paste('pvalue', mean(ci_outclo.sim >= ci_outclo_obs)),
     xlab = "Out-Closeness", xlim = c(low, up))
segments(x0=ci_outclo_obs, x1=ci_outclo_obs, y0=-2, y1=412, col='red', lwd = 2)


low = pmin(min(ci_betw.sim), ci_betw_obs) - 0.05
up = pmax(max(ci_betw.sim), ci_betw_obs) + 0.05
hist(ci_betw.sim, col = "lightgray", main = paste('pvalue', mean(ci_betw.sim >= ci_betw_obs)),
     xlab = "Betweenness", xlim = c(low, up))
segments(x0=ci_betw_obs, x1=ci_betw_obs, y0=-2, y1=302, col='red', lwd = 2)

low = pmin(min(eig.sim), eig_obs) - 0.005
up = pmax(max(eig.sim), eig_obs) + 0.005
hist(eig.sim, col = "lightgray", main = paste('pvalue', mean(eig.sim >= eig_obs)),
     xlab = "Eigenvector Centrality", xlim = c(low, up))
segments(x0=eig_obs, x1=eig_obs, y0=-2, y1=302, col='red', lwd = 2)

# compute an approximate p-value
mean(abs(dens.sim >= density_obs))
mean(abs(trans.sim >= transitivity_obs))
mean(abs(reci.sim >= reciprocity_obs))
mean(abs(ci_indeg.sim >= ci_indeg_obs))
mean(abs(ci_outdeg.sim >= ci_outdeg_obs))
mean(abs(ci_inclo.sim >= ci_inclo_obs))
mean(abs(ci_outclo.sim >= ci_outclo_obs))
mean(abs(ci_betw.sim >= ci_betw_obs))



# 3. Is the observed network coherent with a Binomial random graph model G(n, p)?
# a more formal approach based on the conditional uniform distribution
# ------------------------------------------------------------------------------------
set.seed(1)

B = 5000
m =  sum(Y, na.rm = TRUE)

dens.sim = trans.sim = reci.sim = modu.sim = ass_loc.sim = ass_orglev.sim = ass_ten.sim = c()
ci_indeg.sim = ci_outdeg.sim = ci_inclo.sim = ci_outclo.sim = ci_betw.sim = eig.sim = c()

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
}

# Graphical comparison
par(mfrow = c(3,3))

##NON HA SENSO FARE L'HIST DELLA DENS
low = pmin(min(dens.sim), density_obs) - 0.05
up = pmax(max(dens.sim), density_obs) + 0.05
hist(dens.sim, col = "lightgray", xlim = c(low, up), xlab = NA, 
     main = c("Density", paste('pvalue',mean(dens.sim >= density_obs))))
segments(x0=density_obs, x1=density_obs, y0=-2, y1=322, col='red', lwd = 2)
box()

low = pmin(min(reci.sim), reciprocity_obs) - 0.05
up = pmax(max(reci.sim), reciprocity_obs) + 0.05
hist(reci.sim, col = "lightgray", xlim = c(low, up), xlab = NA,
     main = c("Reciprocity", paste('pvalue', mean(reci.sim >= reciprocity_obs))))
segments(x0=reciprocity_obs, x1=reciprocity_obs, y0=-2, y1=1002, col='red', lwd = 2)
box()

low = pmin(min(trans.sim), transitivity_obs) - 0.05
up = pmax(max(trans.sim), transitivity_obs) + 0.05
hist(trans.sim, col = "lightgray", xlim = c(low, up), xlab = NA,
     main = c('Transitivity', paste('pvalue', mean(trans.sim >= transitivity_obs))))
segments(x0=transitivity_obs, x1=transitivity_obs, y0=-2, y1=902, col='red', lwd = 2)
box()

low = pmin(min(ci_indeg.sim), ci_indeg_obs) - 0.05
up = pmax(max(ci_indeg.sim), ci_indeg_obs) + 0.05
hist(ci_indeg.sim, col = "lightgray",xlim = c(low, up),xlab = NA,
     main = c('In-Degree', paste('pvalue',mean(ci_indeg.sim >= ci_indeg_obs))))
segments(x0=ci_indeg_obs, x1=ci_indeg_obs, y0=-2, y1=2352, col='red', lwd = 2)
box()

low = pmin(min(ci_outdeg.sim), ci_outdeg_obs) - 0.05
up = pmax(max(ci_outdeg.sim), ci_outdeg_obs) + 0.05
hist(ci_outdeg.sim, col = "lightgray", xlim = c(low, up),xlab = NA,
     main = c("Out-Degree", paste('pvalue', mean(ci_outdeg.sim >= ci_outdeg_obs))))
segments(x0=ci_outdeg_obs, x1=ci_outdeg_obs, y0=-2, y1=2332, col='red', lwd = 2)
box()

low = pmin(min(ci_inclo.sim), ci_inclo_obs) - 0.05
up = pmax(max(ci_inclo.sim), ci_inclo_obs) + 0.05
hist(ci_inclo.sim, col = "lightgray", xlim = c(low, up),xlab = NA,
     main = c("In-Closeness", paste('pvalue', mean(ci_inclo.sim >= ci_inclo_obs))))
segments(x0=ci_inclo_obs, x1=ci_inclo_obs, y0=-2, y1=1322, col='red', lwd = 2)
box()

low = pmin(min(ci_outclo.sim), ci_outclo_obs) - 0.05
up = pmax(max(ci_outclo.sim), ci_outclo_obs) + 0.05
hist(ci_outclo.sim, col = "lightgray", xlim = c(low, up), xlab = NA,
     main = c("Out-Closeness", paste('pvalue', mean(ci_outclo.sim >= ci_outclo_obs))))
segments(x0=ci_outclo_obs, x1=ci_outclo_obs, y0=-2, y1=1222, col='red', lwd = 2)
box()

low = pmin(min(ci_betw.sim), ci_betw_obs) - 0.005
up = pmax(max(ci_betw.sim), ci_betw_obs) + 0.005
hist(ci_betw.sim, col = "lightgray", xlim = c(low, up), xlab = NA,
     main = c("Betweenness", paste('pvalue', mean(ci_betw.sim >= ci_betw_obs))))
segments(x0=ci_betw_obs, x1=ci_betw_obs, y0=-2, y1=1502, col='red', lwd = 2)
box()

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
##################################
###### MODEL SRG WITH ERGM #######
##################################
library(network)
library(ergm)

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
# look in more depth
summary(mod0)

mod1 = ergm(net ~ edges + sender + receiver)
summary(mod1)

load("results/P1_sender_reciver.RData")
summary(mod0)

mod2 = ergm(net ~ edges + mutual + 
               nodefactor("location", base = c(1,3,4)) + nodefactor("org_level", base = c(3,4)) + nodefactor("tenure", base = c(2,3,4)) +
               nodematch("location"), control = control.ergm(seed = 1))
summary(mod2)
# 1) is it significant? Yes -> there is a significant difference 
# between Pr(Y_ij = 1) and Pr(Y_ij = 0)  
# That is, these quantities are significantly different from 0.5
# 2) sign? Negative ->  Pr(Y_ij = 1)< 0.5

# higly significant and negative parameter for edges -> there is a lower tendency in observing
# relation in the network than we would expect by chance.
exp(mod0$coef)
# the odds of observing a relation between two randomly 
# selected nodes is about 55% lower than that of not observing it

graph.density(skill)
# or
mean(Y, na.rm = T)
# or, via the expit function
exp(mod0$coef)/(1+exp(mod0$coef))

# goodness of fit
# ----------------------------------

sim = simulate(mod2, burnin = 1000, nsim = 100, verbose = TRUE, seed = 1)

library(intergraph)

fnc = function(xx){
  ig = asIgraph(xx)
  
  den = graph.density(ig)
  rec = reciprocity(ig)
  tr = transitivity(as.undirected(ig, mode = 'collapse'))
  #modu.sim[b] = modularity(g.sim)
  
  #ass_loc.sim[b] = assortativity(g.sim, V(g.sim)$location)
  #ass_orglev.sim[b] = assortativity(g.sim, V(g.sim)$org_level)
  #ass_ten.sim[b] = assortativity(g.sim, V(g.sim)$tenure)
  
  # IN/OUT DEGREE
  ideg = mean(degree(ig, mode = 'in'))
  odeg = mean(degree(ig, mode = "out"))
  
  # IN/OUT CLOSENESS
  iclo = mean(closeness(ig, mode = 'in'))
  oclo = mean(closeness(ig, mode = 'out'))
  
  # BETWEENESS CENTRALITY
  bet = mean(betweenness(ig, directed = T))
  
  return(cbind(den = den,rec =  rec, tr = tr, ideg = ideg, odeg = odeg, iclo = iclo, oclo = oclo,bet = bet))
}
statistics = as.matrix(t(sapply(sim, function(xx){fnc(xx)})))

p_den = mean(statistics[,1]>= density_obs)
p_tr = mean(statistics[,3] >= transitivity_obs)
p_rec = mean(statistics[,2] >= reciprocity_obs)
p_mod_loc = mean(statistics[,4] >= mod_loc_obs)
p_mod_org = mean(statistics[,5] >= mod_lev_obs)
p_mod_ten = mean(statistics[,6] >= mod_ten_obs)
p_ideg = mean(statistics[,7] >= ci_indeg_obs)
p_odeg = mean(statistics[,6] >= ci_outdeg_obs)
p_iclo = mean(statistics[,9] >= ci_inclo_obs)
p_oclo = mean(statistics[,10] >= ci_outclo_obs)
p_bet = mean(statistics[,11] >= ci_betw_obs)
p_eig = mean(statistics[,12] >= eig_obs)

table = as.matrix(round(rbind(density= c(density_obs, mean(statistics[,1]), p_den), 
                              transitivity= c(transitivity_obs, mean(statistics[,2]), p_tr),
                              reciprocity = c(reciprocity_obs, mean(statistics[,3]), p_rec),
                              mod_location = c(mod_loc_obs, mean(statistics[,4]), p_mod_loc),
                              mod_organzation = c(mod_lev_obs, mean(statistics[,5]), p_mod_org),
                              mod_tenure = c(mod_ten_obs, mean(statistics[,6]), p_mod_ten),
                              in_degree = c(ci_indeg_obs, mean(statistics[,7]), p_ideg),
                              out_degree = c(ci_outdeg_obs, mean(statistics[,8]), p_odeg),
                              in_closeness = c(ci_inclo_obs, mean(statistics[,9]), p_iclo),
                              out_closeness = c(ci_outclo_obs, mean(statistics[,10]), p_oclo),
                              betweeness = c(ci_betw_obs, mean(statistics[,11]), p_bet),
                              eigenvector = c(eig_obs, mean(statistics[,12]), p_eig)), 4))


par(mfrow = c(1,3))
low = pmin(min(statistics[,1]), density_obs) - 0.005
up = pmax(max(statistics[,1]), density_obs) + 0.005
hist(statistics[,1], xlab = "Density", main = paste('pvalue', p_den), xlim = c(low, up))
segments(x0=density_obs, x1=density_obs, y0=-2, y1=1082, col='red', lwd = 2)

low = pmin(min(statistics[,3]), transitivity_obs) - 0.05
up = pmax(max(statistics[,3]), transitivity_obs) + 0.05
hist(statistics[,3], xlab = "Transitivity", main =  paste('pvalue', p_tr), xlim = c(low, up))
segments(x0=transitivity_obs, x1=transitivity_obs, y0=-2, y1=1082, col='red', lwd = 2)

low = pmin(min(statistics[,2]), reciprocity_obs) - 0.05
up = pmax(max(statistics[,2]), reciprocity_obs) + 0.05
hist(statistics[,2], xlab = "Reciprocity", main =  paste('pvalue', p_rec), xlim = c(low, up))
segments(x0=reciprocity_obs, x1=reciprocity_obs, y0=-2, y1=1082, col='red', lwd = 2)

low = pmin(min(statistics[,4]), ci_indeg_obs) - 0.05
up = pmax(max(statistics[,4]), ci_indeg_obs) + 0.05
hist(statistics[,4], xlab = "In-Degree", main =  paste('pvalue', p_ideg), xlim = c(low, up))
segments(x0=ci_indeg_obs, x1=ci_indeg_obs, y0=-2, y1=1082, col='red', lwd = 2)

low = pmin(min(statistics[,5]), ci_outdeg_obs) - 0.05
up = pmax(max(statistics[,5]), ci_outdeg_obs) + 0.05
hist(statistics[,5], xlab = "In-Degree", main =  paste('pvalue', p_odeg), xlim = c(low, up))
segments(x0=ci_outdeg_obs, x1=ci_outdeg_obs, y0=-2, y1=1082, col='red', lwd = 2)

hist(unlist(prova[,2]), xlab = "in-degree"); abline(v = mean(degree(skill, mode = "in")), col = "red")
hist(unlist(prova[,3]), xlab = "out-degree"); abline(v = mean(degree(skill, mode = "out")), col = "red")


