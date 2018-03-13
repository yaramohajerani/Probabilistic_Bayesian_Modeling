data{
	int p; //extract sample size parameter from data
	int N; //total number of observations
	int ni[p+1]; //extract linear indices for station data
	vector[N] y; //extract y variable
	vector[N] x; //extract x variable
	vector[p] v; //extract env variable variable
}
parameters{
	real beta1[p]; //vector of p slopes
	real beta1mean; //parameter for mean of slope distribution
	real<lower=0> beta1_sd; //bounded variable for slope sd
	real beta0[p];	//vector of p slopes
	real betaV0; //intercept between intercept and env variable
	real betaV1; //slope between intercept and env variable
	real<lower=0> beta0_sd; //standard deviation of intercept
	real<lower=0> sigma; //measurement noise sd
}
model{
	for(i in 1:p){ //loop over locations
		beta0[i] ~ normal(betaV0 + betaV1*v[i],beta0_sd); //location-specific intercept is a function of env variable at the location
		beta1[i] ~ normal(beta0mean,beta0_sd); //location-specific slope
		y[(ni[i]+1):ni[i+1]] ~ normal(beta0[i] + beta1[i]*x[(ni[i]+1):ni[i+1]], sigma); //likelihood of data, station by station
	}
}
//generated quantities{
//	matrix[N,p] y_pred; //allocate space for prediction of y
//	for(i in 1:p){ //loop over locations
//		y_pred[,i] = beta0[i] + beta1*x; //location level prediction
//	}
//}
