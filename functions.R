## Rearrange the labels on the responses and permute to get statistical significance test. 
permute_response = function(data,sel) {

  data$response[sel]  = sample( data$response[sel]  )
  data$response[!sel] = sample( data$response[!sel] )
  data
}

## Run the linear mixed model with stacked responses; epigenetic aging and chidlhood trauma (or other lme_outomeomes e.g. depression, anxiety etc.)
run_lme = function(data,sel) {
  
  # data = temp
  
  results = rep(NA,6)
  
  mod = try( lme( response~-1+epi+lme_outome, #should I add pubertal stage here?
             random=~-1+epi+lme_outome|subject,
             weights=varIdent(form=~1| epi),
             #corr=corCompSymm(form=~1 | subject/wave),
             data=data, control = list(opt = "optim"),  method = 'ML', na.action = na.exclude),silent=TRUE) #REML
  
  if ( class(mod) != "try-error" ) {
    
    rand_var = var(cbind(data$response[sel], data$response[!sel]) ) #total variance/covariance matrix;  

    subj_var     = getVarCov(mod)[1:2,1:2] #random effect variance/covariance matrix at the subject level; this is also printed in the random effects section of summary(mod).
    results[1:2] = diag(subj_var/rand_var) #percentage of variance attributed to subjects; change within person.
    wave_var     = (rand_var-subj_var)
    results[3:4] = diag( wave_var / rand_var) #percentage of variance attributed to wave; change over time. 
  
    diag(subj_var) = diag(rand_var)
    results[5]     = cov2cor(subj_var)[2,1]
    diag(wave_var) = diag(rand_var)
    results[6]     = cov2cor(wave_var)[2,1]
  }
  
  results
}

##Calculate p=values. 
calc_pval = function( obs_results, per_results, stat ) {
  
  pvals = rep(NA,3)
  names(pvals) = paste0(gsub('_cor','',stat),c("_P_r>0","_P_r<0","_P_not0"))
 per_results = na.omit(per_results)
 n_permutation = sum( !is.na( per_results[,1] ))
  pvals[1] = sum( per_results[,stat] >= obs_results[stat] )           / n_permutation
  pvals[2] = sum( per_results[,stat] <= obs_results[stat] )           / n_permutation
  pvals[3] = sum( abs(per_results[, stat]) >= abs(obs_results[stat])) / n_permutation
  
  pvals
}
