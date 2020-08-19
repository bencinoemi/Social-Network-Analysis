# link del dataset https://toreopsahl.com/datasets/#Cross_Parker
library(igraph)
library(RColorBrewer)


setwd('/home/noe/Università/in_corso/SocialNetworkAnalysis/')

# ---------------------------
# LOAD NETWORK AND ATTRIBUTES
# ---------------------------

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
 
cols = c(brewer.pal(4, 'Pastel2'))
V(skill)$color = cols[V$location]
V(skill)$names = V
shapes = c("star", "circle", "square", "triangle")
V(skill)$shape = shapes[V$org_level]

# --------------------
# PLOTTING THE NETWORK
# --------------------

mytriangle <- function(coords, v=NULL, params) {
   vertex.color <- params("vertex", "color")
   if (length(vertex.color) != 1 && !is.null(v)) {
      vertex.color <- vertex.color[v]
   }
   vertex.size <- 1/150 * params("vertex", "size")
   if (length(vertex.size) != 1 && !is.null(v)) {
      vertex.size <- vertex.size[v]
   }
   
   symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
           stars=cbind(vertex.size, vertex.size, vertex.size),
           add=TRUE, inches=FALSE)
}
add_shape("triangle", clip=shapes("circle")$clip,
          plot=mytriangle)

mystar <- function(coords, v=NULL, params) {
   vertex.color <- params("vertex", "color")
   if (length(vertex.color) != 1 && !is.null(v)) {
      vertex.color <- vertex.color[v]
   }
   vertex.size  <- 1/150 * params("vertex", "size")
   if (length(vertex.size) != 1 && !is.null(v)) {
      vertex.size <- vertex.size[v]
   }
   norays <- params("vertex", "norays")
   if (length(norays) != 1 && !is.null(v)) {
      norays <- norays[v]
   }
   
   mapply(coords[,1], coords[,2], vertex.color, vertex.size, norays,
          FUN=function(x, y, bg, size, nor) {
             symbols(x=x, y=y, bg=bg,
                     stars=matrix(c(size,size/2), nrow=1, ncol=nor*2),
                     add=TRUE, inches=FALSE)
          })
}
add_shape("star", clip=shape_noclip,
          plot=mystar, parameters=list(vertex.norays=5))

l = layout_with_lgl(skill)
c = norm_coords(l, ymin=-1, ymax=1, xmin=-2, xmax=2)
set.seed(1)
plot(skill, vertex.size = V(skill)$tenure*5, vertex.shape = shapes, rescale = F,
     layout = c, edge.arrow.size = 0.15)

# ------------------
# CENTRAL STATISTICS
# ------------------
skill_adj = get.adjacency(skill, sparse = F)
diag(skill_adj) = NA
skill_adj

# -----------------
# DEGREE CENTRALITY
# -----------------
# in-degree
# ---------
ziind = colSums(skill_adj, na.rm = TRUE)
summary(ziind)

# in_degree standardized
ziind_st = ziind/(n-1)
summary(ziind_st)

# most popular nodes
most_pop = V(skill)$names[order(ziind_st, decreasing = T)[1:3]]

# out_degree
# ----------
zioutd = rowSums(skill_adj, na.rm = TRUE)
summary(zioutd)

# out_degree standardized
zioutd_st = zioutd/(n-1)
summary(zioutd_st)

# most expansive nodes
most_exp = V(skill)$names[order(zioutd_st, decreasing = T)[1:3]]

# PLOT BY DEGREE
par(mfrow = c(1,1))
central = union(most_exp, most_pop)
V(skill)$shape = 'circle'
V(skill)$color[order(ziind_st, decreasing = T)] = c(heat.colors(n))

V(skill)$labelfont = rep(1, n)
V(skill)$labelfont[central] = 2
set.seed(1)
plot(skill, vertex.size = zioutd_st*40, vertex.label.font = V(skill)$labelfont,
     vertex.label.color = 'black', edge.arrow.size = 0.15, layout = c, 
     vertex.color = c(heat.colors(1000, rev = T)[1000*ziind_st]))

# --------------------
# CLOSENESS CENTRALITY
# --------------------

dist_in = distances(skill, mode = 'in')
diag(dist_in) = NA

in_clo = 1/rowSums(dist_in, na.rm = T)

# in-closeness
ziinc = closeness(skill, mode = 'in')
summary(ziinc)

# in-closeness standardized
in_clo_st = in_clo*(n-1)
ziinc_st = closeness(skill, mode = 'in', normalized = T)
summary(ziinc_st)

# most in_close nodes
most_cloin = V(skill)$names[order(ziinc_st, decreasing = T)[1:3]]

# out-closeness
dist_c = distances(skill, mode = 'out')
clo_o = 1/rowSums(dist_c)

zioutc = closeness(skill, mode = 'out')
summary(zioutc)

# out-closeness standardized
zioutc_st = closeness(skill, mode = 'out', normalized = T)
summary(zioutc_st)

# most out_close nodes
most_cloout = V(skill)$names[order(zioutc_st, decreasing = T)[1:3]]

# PLOT BY CLOSENESS
uni = union(most_cloin, most_cloout)

V(skill)$labelfont = rep(1, n)
V(skill)$labelfont[uni] = 2

set.seed(1)
plot(skill, vertex.size = zioutc_st*40, edge.arrow.size = 0.15, layout = c, 
     vertex.label.font = V(skill)$labelfont, vertex.label.color = 'black',
     vertex.color = c(heat.colors(100, rev = T)[round(100*ziinc_st)]))

# ---------------------
# BETWEENESS CENTRALITY
# ---------------------

zib = betweenness(skill, directed = T)
summary(zib)

# standardized betweenness
zib_st = betweenness(skill, directed = T, normalized = T)
summary(zib_st)

# most betweenness central nodes
most_betw = V(skill)$names[order(zib_st, decreasing = T)[1:3]]

V(skill)$color = 'lightblue'

V(skill)$labelfont = 1
V(skill)$labelfont[most_betw] = 2

set.seed(1)
par(mfrow = c(1,1), mar = c(0,0,2,0))
plot(skill, vertex.size = zib_st*40, edge.arrow.size = 0.15,layout = c, 
     vertex.label.font = V$labelfont, vertex.label.color = 'black', vertex.label.cex = 1.5)

# ---------------------
# EIGENVECTOR CENTRALITY
# ---------------------

zie = eigen_centrality(skill, directed = T, scale = F)$vector
summary(zie)

# standardized eigenvector centrality
zie_st = eigen_centrality(skill, directed = T, scale = T)$vector
summary(zie_st)

# most central nodes
most_eig = order(zie_st, decreasing = T)[1:3]

V(skill)$labelfont = 1
V(skill)$labelfont[most_eig] = 2
set.seed(1)
par(mfrow = c(1,1), mar = c(0,0,2,0))
plot(skill, vertex.size = zie_st*40, edge.arrow.size = 0.15,layout = c,
     vertex.label.font = V$labelfont, vertex.label.color = 'black')

# PLOT BY EIGENVECTOR AND BETWEENESS
uni = union(most_betw, most_eig)

V(skill)$labelfont = rep(1, n)
V(skill)$labelfont[uni] = 2

set.seed(1)
plot(skill, vertex.size = zie_st*40, edge.arrow.size = 0.15, layout = c, 
     vertex.label.font = V(skill)$labelfont, vertex.label.color = 'black',
     vertex.color = c(heat.colors(100, rev = T)[round(100*zib_st+1)]))

cbind(V(skill)$names, zib_st, zib_st)

# ALL SUMMARY TOGETHER
# --------------------
cbind(in_deg = summary(ziind_st), out_deg = summary(zioutd_st), clo_in = summary(ziinc_st),
      clo_out = summary(zioutc_st), betw = summary(zib_st), eig = summary(zie_st))

# ----------------------
# CENTRALIZATION INDECES
# ----------------------
max_indeg = max(ziind)
ind_ci = sum(max_indeg - ziind) / (n-1)^2
ci_indeg = centr_degree(skill, loops = F, mode = 'in')$centralization
c(ci_indeg, ind_ci)


max_oudeg = max(zioutd)
out_ci = sum(max_oudeg - zioutd) / (n-1)^2
ci_outdeg = centr_degree(skill, loops = F, mode = 'out')$centralization
c(ci_outdeg, out_ci)

max_inclo = max(ziinc)
inc_ci = sum(max_inclo - ziinc) 
ci_inclo = centr_clo(skill, mode = 'in')$centralization
c(ci_inclo, inc_ci)


max_ouclo = max(zioutc)
ouc_ci = sum(max_ouclo - zioutc)
ci_outclo = centr_clo(skill, mode = 'out')$centralization
c(ci_outclo, ouc_ci)

bet_ci = sum(max(zib) - zib) / ((n-1)^2*(n-2))
ci_betw = centr_betw(skill, directed = T)$centralization
c(ci_betw, bet_ci)

cbind(in_deg = ci_indeg, out_deg = ci_outdeg, clo_in = ci_inclo,
      clo_out = ci_outclo, betw = ci_betw)

# ----------------------------
# IDENTIKIT MOST CENTRAL NODES
# ----------------------------
cbind(in_deg = most_pop, out_deg = most_exp, clo_in = most_cloin,
      clo_out = most_cloout, betw = most_betw, eig = most_eig)

loc = c('Paris', 'Frankfurt', 'Warsaw', 'Geneva')
ten = c('1-12 months', '13-36 months', '37-60 months', '61+ months') 
org = c('Global Dept Manager', 'Local Dept Manager', 'Project Leader', 'Researcher')

cbind(most_pop_nodes = most_pop,
      location = loc[V(skill)$location[most_pop]],
      tenure = ten[V(skill)$tenure[most_pop]], 
      org_level = org[V(skill)$org_level[most_pop]] )

cbind(most_exp_nodes = most_exp,
      location = loc[V(skill)$location[most_exp] ],
      tenure = ten[V(skill)$tenure[most_exp]], 
      org_level = org[V(skill)$org_level[most_exp]]) 

cbind(clo_in_nodes = most_cloin,
      location = loc[V(skill)$location[most_cloin]],
      tenure = ten[V(skill)$tenure[most_cloin] ],
      org_level = org[V(skill)$org_level[most_cloin]]) 

cbind(clo_out_nodes = most_cloout,
      location = loc[V(skill)$location[most_cloout]],
      tenure = ten[V(skill)$tenure[most_cloout]],  
      org_level = org[V(skill)$org_level[most_cloout]])

cbind(most_betw_nodes = most_betw,
      location = loc[V(skill)$location[most_betw]],
      tenure = ten[V(skill)$tenure[most_betw]],
      org_level = org[V(skill)$org_level[most_betw]])

cbind(most_eig_nodes = most_eig,
      location = loc[V(skill)$location[most_eig]],
      tenure = ten[V(skill)$tenure[most_eig]], 
      org_level = org[V(skill)$org_level[most_eig]])



