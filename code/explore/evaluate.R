volume_ssd <- function(v,act,G){
  # Return the volume of the symmetric set difference between
  # a (series of) clusters v[1], v[1:2], v[1:3], ...
  # and the actual cluster act.
  
  # v: ordering (i.e. ranking) of vertices 1:n
  # act: subset of 1:n
  # G: graph
  
  # vectorized for speed
  loss_denom <- sum(degree(G)[act])
  
  belongs_to_act <- ifelse(v %in% act,1,-1)
  deg_v <- degree(G)[v]
  loss_num <- sum(degree(G)[act]) + cumsum(-deg_v*belongs_to_act)
  loss_num/loss_denom
}

cardinality_ssd <- function(v,act,G){
  # Return the cardinality of the symmetric set difference between
  # a (series of) clusters v[1], v[1:2], v[1:3], ...
  # and the actual cluster act.
  
  # v: ordering (i.e. ranking) of vertices 1:n
  # act: subset of 1:n
  # G: graph
  loss_denom <- length(act)
  
  belongs_to_act <- ifelse(p_order %in% act,1,-1)
  loss_num <- length(act) + cumsum(-belongs_to_act)
  loss_num/loss_denom
}