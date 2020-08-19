library(igraph)
library(network)
library(ergm)
library(RColorBrewer)
# --------------------
# Load the data
# --------------------
setwd("G:/Social Network Analysis/Progetto/Dati")
data = read.table("skill.txt")

edgelist = as.matrix(data[, c(1,2)])
skill_w = graph_from_edgelist(edgelist, directed = T)

weight = data[,3]
sum(weight == 0)     #check for missing values (no one)
bin.weight = as.numeric(weight>3) #binary weights
skill = delete_edges(skill_w, E(skill_w)[bin.weight == 0])

V(skill)$name = seq(1,77)
V(skill)$location = scan("/location.txt")
V(skill)$org_level = scan("/org_level.txt")
V(skill)$tenure = scan("/tenure.txt")

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



######################################
Y = get.adjacency(skill, sparse = F)
diag(Y) = NA

# number of nodes
n = nrow(Y)
net = network(Y, directed = T)
net

# add attributes
net %v% 'location' = scan("location.txt")
net %v% 'org_level' = scan("org_level.txt")
net %v% 'tenure' = scan("tenure.txt")

####################
# MODEL ESTIMATION #
####################


#Model1: Best Markov Model
mod1 = ergm(net ~ edges + mutual + ostar(2) + ostar(3) + triangle + 
              nodefactor("location") + nodematch("location"),
            estimate = "MPLE",
            control = control.ergm(seed = 1))

summary(mod1)


#First looking for Best BIC model, so introducing ostars and triangles

#Model2: Social Model: add gwesp and gwdsp

mod2 = ergm(net ~ edges + mutual + ostar(2) + ostar(3) + triangle +
                 nodematch("location") + nodefactor("location") +
                 gwesp(fixed = T) + gwdsp(fixed = T),
               estimate = "MPLE",
               control = control.ergm(seed=1))


summary(mod2)

#Model3: Remove  gwesp(fixed = T) from Model2
mod3 = ergm(net ~ edges + mutual + ostar(2) + ostar(3) + triangle +
                 nodematch("location") + nodefactor("location") +
                 gwdsp(fixed = T),
               estimate = "MPLE",
               control = control.ergm(seed=1))


summary(mod3)
BIC(mod1, mod2, mod3)


#Mod3 Best BIC model; now try to estimate models with MCMC algorithm; in this way
#we can get simulations later, but we have to remove stars and triangles parameters.

#------------------------------------------------------------------------------#

#Model4: 

mod4 = ergm(net ~ edges + mutual +
               nodefactor("location") + nodematch("location") + 
               #gwidegree(decay = 1, fixed = TRUE) + gwodegree(decay = log(2), fixed = TRUE) +
               gwesp(fixed = T) + gwdsp(fixed = T),
            estimate = "MLE",
            control = control.ergm(seed=1))


#No result



#Model5: Remove all except typical configuration of the model
mod5 = ergm(net ~ edges +
              gwesp(fixed = T) + gwdsp(fixed = T),
            control = control.ergm(seed=1, main.method = "Robbins-Monro"))


summary(mod5)





sim = simulate(mod3, burnin = 1000, nsim = 1000, verbose = TRUE, seed = 1)

# let us assume we want to verify whether the model is appropriate to represent the degree and the 
# transitivity in the network

# install.packages("intergraph")
library(intergraph)
?asIgraph

fnc = function(xx){
  ig = asIgraph(xx)
  tr = transitivity(ig)
  ideg = mean(degree(ig, mode = "in"))
  odeg = mean(degree(ig, mode = "out"))
  return(list(tr, ideg, odeg))
}
prova = as.matrix(t(sapply(sim, function(xx){fnc(xx)})))
dev.new()
par(mfrow = c(3,1))
hist(unlist(prova[,1]), xlab = "transitivity"); abline(v = transitivity(friend.net), col = "red")
hist(unlist(prova[,2]), xlab = "in-degree"); abline(v = mean(degree(friend.net, mode = "in")), col = "red")
hist(unlist(prova[,3]), xlab = "out-degree"); abline(v = mean(degree(friend.net, mode = "out")), col = "red")
