
correct_integer<- function(data_){
  if(!is.integer(data_)){
    return( apply(data_, 2, function(x) round(x,0)))
  }else{
    return (data_)
  }
}

correct_continuous <- function(data_){
  data_[data_ ==0] =  0.001
  return(data_)
}
correct_domain <- function(data_, c_){
  if(c_ == 'c'){return(correct_continuous(data_))
  }else{
    return(correct_integer(data_))
  }
}



gene_trueLV <- function(q_star, n,  p, power_){
  #eigen_list = seq((p + q_star) * (n^(power_)), (p+1)*(n^(power_)), by = -n^(power_))
  eigen_list = seq((q_star) * (n^(power_)), (1)*(n^(power_)), by = -n^(power_))
  d_ = sqrt(eigen_list)
  sim_U = rustiefel(n, q_star)
  D_matrix = diag(d_)
  sim_L = tcrossprod(sim_U, D_matrix)
  sim_V = rustiefel(p, q_star)

  eta_star = tcrossprod(sim_L, sim_V)
  param_star = dmf_identify(eta_star, q_star)
  return (param_star)
}
