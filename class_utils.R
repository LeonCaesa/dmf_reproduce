
TrainTest_Flag<-function(plot_df, seed_, train_ratio =0.7){
  train_size = floor(train_ratio * nrow(plot_df))
  set.seed(seed_)
  split1<- sample(c(rep(0, train_size), rep(1, nrow(plot_df)- train_size)))
  list(train = split1==0, test = split1==1)}

TrainTest_Split<-function(plot_df, splitflag){
  p = dim(plot_df)[2]; colnames(plot_df)[p] = 'y'
  return (list(train = plot_df[splitflag$train, ], test = plot_df[!splitflag$test, ]))}

Tree_tuned<-function(data){
  tree_fit = rpart(y ~ ., data = data, parms = list(split ="information"), 
                   control = c(cp = 0, xval =10), method= 'class') 
  
  best_cp = tree_fit$cptable[which.min(tree_fit$cptable[, "xerror"]),"CP"]
  tree_fit = prune(tree_fit, best_cp)
  return(tree_fit)}

center_identify = function(fit_dmf){
  dmf_centered = dmf_center(fit_dmf)
  eta_hat = tcrossprod(dmf_centered$L, dmf_centered$V)
  dmf_identified = dmf_identify(eta_hat, dim(dmf_centered$L)[2])
  return (list(L = dmf_identified$L, V = dmf_identified$V, family = fit_dmf$family))
}
