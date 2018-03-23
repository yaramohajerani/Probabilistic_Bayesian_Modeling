data{
	int p; //extract sample size parameter from data
	int N; //total number of observations
	int NI; //number of environmental variables for the intercept
	int NS; //number of environmental variables for the slope
	int ni[p+1]; //linear indices for station data
	vector[N] y; //y variable
	vector[N] x; //x variable
	matrix[p,NI+1] M0; //environmental variables for logJ0 (NI+1 for column of 1s)
	matrix[p,NS+1] M1; //environmental variables for b (NS+1 for column of 1s)
}
parameters{
	real<lower=-10,upper=10> logJ0[p]; //vector of p intercepts
	real<lower=0,upper=3> b[p];	//vector of p slopes
	real<lower=1e-15,upper=10> logJ0_sd; //bounded variable for intercept sd
	real<lower=1e-15,upper=3> b_sd; //standard deviation of slopes
	vector[NI+1] beta0; //coefficients for logJ0 environmental dependence
	vector[NS+1] beta1; //coefficients for b environmental dependence
	real<lower=1e-15> sigma; //measurement noise sd
}
model{
	for(i in 1:p){ //loop over locations
		logJ0[i] ~ normal(dot_product(M0[i,],beta0),logJ0_sd); //location-specific intercept is a function of env variable at the location
		b[i] ~ normal(dot_product(M1[i,],beta1),b_sd); //location-specific intercept is a function of env variable at the location
		y[(ni[i]+1):ni[i+1]] ~ normal(logJ0[i] - b[i]*x[(ni[i]+1):ni[i+1]], sigma); //likelihood of data, station by station
	}
}
