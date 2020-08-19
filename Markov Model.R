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



library(ergm)


# -------------------------------------------
# let us move towards the Markov graph model
# -------------------------------------------
# Let us add to the model the triangle term, the in- and the out-stars of order 2
# (indicating the tendency to form clusters in the network)


#Model 1: principal configuration included, estimate MPLE
mod1 = ergm(net ~ edges + istar(2), control = control.ergm(seed = 1))
summary(mod1)

mod2 = ergm(net ~ edges + ostar(2), control = control.ergm(seed = 1))
summary(mod2)

mod3 = ergm(net ~ edges + ostar(2) +
              nodefactor("location") + nodefactor("org_level") + nodefactor("tenure") +
              nodematch("location"), control = control.ergm(seed = 1))
summary(mod3)

mod4 = ergm(net ~ edges + istar(2) +
              nodefactor("location") + nodefactor("org_level") + nodefactor("tenure") +
              nodematch("location"), control = control.ergm(seed = 1))
summary(mod4)

mod4 = ergm(net ~ edges + ostar(2) + ostar(3) +
              nodefactor("location") + nodefactor("org_level") + nodefactor("tenure") +
              nodematch("location"), control = control.ergm(seed = 1))
summary(mod4)
#save(mod4, file="results/markov_no_agg.RData")

mod4a = ergm(net ~ edges + ostar(2) + ostar(3) + 
              nodefactor("location", base = c(1,3,4)) + nodefactor("org_level", base = c(3,4)) + nodefactor("tenure", base = c(2,3,4)) +
              nodematch("location"), control = control.ergm(seed = 1))
summary(mod4a)
#save(mod4a, file="results/Markov_agg.RData")
dev.new()
pdf(file = "./diagnostics_mod4a.pdf")
mcmc.diagnostics(mod4a)
dev.off()

mod5 = ergm(net ~ edges + ostar(2) + ostar(3) + 
              nodefactor("location") + nodefactor("org_level") + nodefactor("tenure") +
              nodematch("location")+ gwesp(fixed = T), control = control.ergm(seed = 1))
summary(mod5)

#save(mod5, file="results/social_no_agg.RData")

mod5a = ergm(net ~ edges + ostar(2) + ostar(3) + 
              nodefactor("location", base=c(1,3,4)) + nodefactor("org_level", base=c(3,4)) + 
              nodefactor("tenure", base = c(2,3,4)) + nodematch("location") +
              gwesp(fixed = T), control = control.ergm(seed = 1))
summary(mod5a)
mod5a$mle.lik

#save(mod5a, file="results/social_agg.RData")


