load('C:/Users/ManvendraSingh/Downloads/TRACERx_NEJM_2017.rda')
tx <- TRACERx_NEJM_2017
pats <- sort(unique(tx$patientID))
cat('unique patients =', length(pats), '\n')
print(pats)
