library(nlme)

work_dir = 'G:\\My Drive\\SOP\\methylation\\laura\\mixmod_dnam\\'

outcome        = "depr"  # outcome variables
covs          = c("Sex","obs_age","obs_age2","race_aa","CD03","CD14","CD15","min_thres","max_thres","peak4","peak5","PeakSQRT","AUC","skewness")
n_permutation = 20
source(paste0(work_dir,"functions.R"))

data = read.csv(paste0(work_dir,"data_long_trauma_20190730.csv"),header=T)

data[,predVAR]
# starts
per_results           = matrix(NA,n_permutation,6)
colnames(per_results) = c("subj_epiAge","subj_outcome","wave_epiAge","wave_outcome","subj_cor","wave_cor") 

sel = data$epi == 1
model               = lm(data$response[ sel] ~ as.matrix(data[sel,covs]))
data$response[sel]  = model$residuals
model               = lm(data$response[!sel] ~ as.matrix(data[!sel,covs]))
data$response[!sel] = model$residuals
colnames(data)[ colnames(data) %in% outcome] = "lme_outome"

# observed results
obs_results         = run_lme(data,sel)
names(obs_results)  = colnames(per_results)
obs_results  

# p values
for (i in 1:n_permutation) { # i=2
  if(i %% 10==0) cat(paste0("Permutation: ", i, "\n"))
  temp            = permute_response(data,sel)
  per_results[i,] = run_lme(temp,sel)
  
}

per_results = na.omit(per_results)
summary(per_results)
subj_p = calc_pval(obs_results, per_results, "subj_cor" )
wave_p = calc_pval(obs_results, per_results, "wave_cor" )
# permutation P values
c( subj_p, wave_p )

