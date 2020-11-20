
rm(list=ls()) 

work_dir = "G:\\My Drive\\SOP\\methylation\\laura\\mixmod_dnam\\code\\"

data     = read.csv(paste0(work_dir,"long_trauma.csv"))

### Data should have long format with following structure
# Subject wave type	      response	epi	correlate
# A00027	1    correlate    2	      0	    1
# A00027	1    DNAm age     22.02	  1	    0
# A00027	2    correlate	  1	      0	    1
# A00027	2    DNAm age	    26.84   1	    0
### thus, response contains the values of the correlate and DNAm age and epi/correlate are elements of the design matrix. 

covs  = c("tanner_stage","Sex","obs_age2.1","obs_age2.2", "race","CD03","CD14","CD15","min_thres","max_thres","peak4","peak5","PeakSQRT","AUC","skewness")

n_permutation = 10

# begin
library(nlme)
set.seed(1)
source(paste0(work_dir,"functions.R"))
colnames_results = c("subj_epi","subj_correlate","wave_epi","wave_correlate","subj_contribution","wave_contribution") 


# prepare data for permutations
data           = prep_permute_response_stratify(data)
sample_records = unique( data[data$rownr==1,c("subject","nrecords")] )
records        = unique(data$nrecords)

# residualize response to speed up permutations (has little impact on results) 
sel = data$epi == 1
model               = lm(data$response[ sel] ~ as.matrix(data[sel,covs]))
data$response[sel]  = model$residuals
model               = lm(data$response[!sel] ~ as.matrix(data[!sel,covs]))
data$response[!sel] = model$residuals

### observed results
obs_results         = run_lme(data,correlate,sel)
names(obs_results)  = colnames_results
obs_results  

# Perform permutations that preserve subject-level variance 
per_results           = matrix(NA,n_permutation,6)
colnames(per_results) = colnames_results 
for (i in 1:n_permutation) { # i=1
  
  perm_data       = permute_response_stratify(data,sample_records,records) 
  per_results[i,] = run_lme(perm_data,correlate,sel)
  
}

# Calculate p values. 
per_results = na.omit(per_results)
subj_p      = calc_pval(obs_results, per_results, "subj_contribution" )
wave_p      = calc_pval(obs_results, per_results, "wave_contribution" )
pcalc       = rbind(subj_p,wave_p)
pcalc
