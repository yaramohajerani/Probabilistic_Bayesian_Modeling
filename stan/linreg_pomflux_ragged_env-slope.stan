data{
	int p; //extract sample size parameter from data
	int N; //total number of observations
	int ni[p+1]; //extract linear indices for station data
	vector[N] y; //extract y variable
	vector[N] x; //extract x variable
	vector[p] v; //extract env variable variable
}
parameters{
	real beta0[p]; //vector of p intercepts
	real beta0mean; //parameter for mean of intercept distribution
	real<lower=0> beta0_sd; //bounded variable for intercept sd
	real beta1[p];	//vector of p slopes
	real betaV0; //intercept between slope and env variable
	real betaV1; //slope between slope and env variable
	real<lower=0> beta1_sd; //standard deviation of slopes
	real<lower=0> sigma; //measurement noise sd
}
model{
	for(i in 1:p){ //loop over locations
		beta1[i] ~ normal(betaV0 + betaV1*v[i],beta1_sd); //location-specific slope is a function of env variable at the location
		beta0[i] ~ normal(beta0mean,beta0_sd); //location-specific intercept
		y[(ni[i]+1):ni[i+1]] ~ normal(beta0[i] + beta1[i]*x[(ni[i]+1):ni[i+1]], sigma); //likelihood of data, station by station
	}
}
//generated quantities{
//	matrix[N,p] y_pred; //allocate space for prediction of y
//	for(i in 1:p){ //loop over locations
//		y_pred[,i] = beta0[i] + beta1*x; //location level prediction
//	}
//}
