
#### distance matrix for location i and its neighbors ####
i_dist <- function(i, neighbor_index, s)	dist(s[c(i, neighbor_index[[i - 1]]), ])

get_neardistM <- function (ind, ind_distM_d) {
  if (ind < M ){l = ind } else {l = M}; 
  M_i <- rep(0, M * (M - 1) / 2);
  if (l == 1) {}
  else{
    M_i[1: (l * (l - 1) / 2)] <- 
      c(ind_distM_d[[ind]])[(l + 1): (l * (l + 1) / 2)]
  }
  return(M_i)
}

get_neardist <- function (ind, ind_distM_d) {
  if (ind < M ){l = ind } else {l = M}; 
  D_i <- rep(0, M);
  D_i[1:l]<-c(ind_distM_d[[ind]])[1:l]
  return( D_i)
}

get_nearind <- function (ind, ind_distM_i) {
  if (ind < M ){l = ind } else {l = M}; 
  D_i <- rep(0, M);
  D_i[1:l]<-c(ind_distM_i[[ind]])[1:l]
  return( D_i)
}


#### wrap up in function NNMatrix ####

NNMatrix <- function(N, coords.ord, n.indx){
    
    nearind <- t(sapply(1: (N - 1), get_nearind, n.indx))
    neighbor_dist <- sapply(2:N, i_dist, n.indx, coords.ord)
    neardistM <- t(sapply(1: (N - 1), get_neardistM, neighbor_dist))
    neardist <- t(sapply(1: (N - 1), get_neardist, neighbor_dist))
    
    return(list(nearind = nearind, neardistM = neardistM, neardist = neardist))
}



