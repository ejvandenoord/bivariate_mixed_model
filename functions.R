prep_permute_response_stratify = function(data,max_records) {
  
  # for each subject/wave both epi and correlate needs to be present 
  data = data[ order(data[,"subject"],data[,"wave"],data[,"epi"]) ,]
  
  temp            = rle( as.vector(data[data$epi==1,"subject"]) )
  gsmsid_records  = temp$lengths
  ind  = match(data[,"subject"],temp$values)
  data = cbind(data,gsmsid_records[ind])
  colnames(data)[ncol(data)] = 'nrecords'
  data = data[data$nrecords<=max_records,]
  
  temp        = rle( as.vector(data[,"subject"]) )
  data$rownr  = unlist(lapply(temp$lengths, function(x) seq(1,x)))
  
  #    write.csv(data,"temp.csv")
  data
  
}

permute_response_stratify = function(data,sample_records,records) {
  
  permute = function(temp) {
    subject      = sample_records[sample_records$nrecords == records[i],"subject"]
    new_subject  = subject[sample( 1:length(subject) )]
    ind          = match(new_subject,temp$subject)
    ind          = ind[ !is.na(ind) ]
    if (i>1) ind = unlist( lapply(ind, function(x) {
      x= seq(x,x+(i-1)) 
      x[order(runif(i))]
      x} ))
    ind
  }
  
  for (i in seq_along(records)) { # i = 2
    
    sel          = data$epi == 1 & data$nrecords == records[i]
    temp         = data[sel,]
    ind          = permute(temp)
    data[ sel ,"response"] = temp[ind,"response"]   
    
    sel          = data$epi == 0 & data$nrecords == records[i]
    temp         = data[sel,]
    ind          = permute(temp)
    data[ sel ,"response"] = temp[ind,"response"]   
  }     
  
  data
  
}         


## Rearrange the labels on the responses and permute to get statistical significance test. 
permute_response = function(data,sel) {
  data$response[sel]  = sample( data$response[sel]  )
  data$response[!sel] = sample( data$response[!sel] )
  data
}

## Run the linear mixed model with stacked responses; epigenetic aging and childhood trauma (or other correlates e.g. depression, anxiety etc.)
run_lme = function(data,correlate,sel) {
  
  # data = perm_data
  
  results = rep(NA,6)
  
  mod = try( lme(response~-1+epi+correlate,
             random=~-1+epi+correlate|subject,
             weights=varIdent(form=~1| epi),
             #corr=corCompSymm(form=~1 | subject/wave),
             #data=data, control = list(opt = "optim"),  method = 'ML', na.action = na.exclude),silent=TRUE) #REML
             data=data, control = lmeControl(MaxIter = 500, msVerbose = F),  method = 'ML', na.action = na.exclude),silent=TRUE) #REML
  
  if ( class(mod) != "try-error" ) {
    
    rand_var     = var(cbind(data$response[sel], data$response[!sel]) ) #total variance/covariance matrix;  
    subj_var     = getVarCov(mod)[1:2,1:2] #random effect variance/covariance matrix at the subject level; this is also printed in the random effects section of summary(mod).
    results[1:2] = diag(subj_var/rand_var) #proportion of variance attributed to subjects;
    wave_var     = (rand_var-subj_var)
    results[3:4] = diag( wave_var / rand_var) #proportion of variance attributed to wave;
  
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
  names(pvals) = c("P_r>0","P_r<0","P_not0")
  per_results <- na.omit(per_results)
  n_permutation = sum( !is.na( per_results[,1] ))
  pvals[1] = sum( per_results[,stat] >= obs_results[stat] )           / n_permutation
  pvals[2] = sum( per_results[,stat] <= obs_results[stat] )           / n_permutation
  pvals[3] = sum( abs(per_results[, stat]) >= abs(obs_results[stat])) / n_permutation
  
  pvals
}