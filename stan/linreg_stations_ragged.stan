data{
	int p; //extract sample size parameter from data
	int N; 
	int n[p]; //extract number of stations parameter
	int ni[p+1]; //extract linear indices for station data
	vector[N] y; //extract y variable
	vector[N] x; //extract x variable
}
parameters{
	real beta0[p]; //declare vector pdf for intercepts
	real beta0mean; //declare parameter for mean of intercept distribution
	real<lower=0> beta0_sd; //declare bounded variable for intercept sd
	real beta1;	//declare single slope parameter
	real<lower=0> sigma; //declare measurement noise sd
}
model{
	for(i in 1:p){ //loop over locations
		beta0[i] ~ normal(beta0mean,beta0_sd); //location-specific intercept
		y[(ni[i]+1):ni[i+1]] ~ normal(beta0[i] + beta1*x[(ni[i]+1):ni[i+1]], sigma); //likelihood of data, station by station
	}
}
//generated quantities{
//	matrix[N,p] y_pred; //allocate space for prediction of y
//	for(i in 1:p){ //loop over locations
//		y_pred[,i] = beta0[i] + beta1*x; //location level prediction
//	}
//}

