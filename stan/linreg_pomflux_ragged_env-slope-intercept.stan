data{
	int p; //extract sample size parameter from data
	int N; //total number of observations
	int ni[p+1]; //extract linear indices for station data
	vector[N] y; //extract y variable
	vector[N] x; //extract x variable
	vector[p] v0; //extract env variable 0 (for intercept)
	vector[p] v1; //extract env variable 1 (for slope)
}
parameters{
	real beta0[p]; //vector of p intercepts
	real beta1[p];	//vector of p slopes
	real<lower=0> beta0_sd; //bounded variable for intercept sd
	real<lower=0> beta1_sd; //standard deviation of slopes
	real beta0V0; //intercept between intercept and env variable
	real beta0V1; //slope between intercept and env variable
	real beta1V0; //intercept between slope and env variable
	real beta1V1; //slope between slope and env variable
	real<lower=0> sigma; //measurement noise sd
}
model{
	for(i in 1:p){ //loop over locations
		beta1[i] ~ normal(beta1V0 + beta1V1*v1[i],beta1_sd); //location-specific slope is a function of env variable at the location
		beta0[i] ~ normal(beta0V0 + beta0V1*v0[i],beta0_sd); //location-specific intercept is a function of env variable at the location
		y[(ni[i]+1):ni[i+1]] ~ normal(beta0[i] + beta1[i]*x[(ni[i]+1):ni[i+1]], sigma); //likelihood of data, station by station
	}
}
//generated quantities{
//	matrix[N,p] y_pred; //allocate space for prediction of y
//	for(i in 1:p){ //loop over locations
//		y_pred[,i] = beta0[i] + beta1*x; //location level prediction
//	}
//}
