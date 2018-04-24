data{
	int p; //extract sample size parameter from data
	int N; //total number of observations
	int Ninf; //number of environmental variables for  bInf
	int N1; //number of environmental variables for b1
	int N2; //number of environmental variables for b2
	int ni[p+1]; //linear indices for station data
	vector[N] y; //y variable
	vector[N] x; //x variable
	matrix[p,Ninf+1] MInf; //environmental variables for bInf
	matrix[p,N1+1] M1; //environmental variables for b1
	matrix[p,N2+1] M2; //environmental variables for b2
}
parameters{
	real bInf[p]; //vector of p intercepts (at infinity)
	real b1[p];	//vector of len p for coefficient of exponential
	real b2[p];	//vector of len p for denominator in exponential
	real<lower=1e-15,upper=10> bInf_sd; //bounded variable for intercept bInf
	real<lower=1e-15,upper=10> b1_sd; //standard deviation of b1
	real<lower=1e-15,upper=10> b2_sd; //standard deviation of b2
	vector[Ninf+1] betaInf; //coefficients for bInf environmental dependence
	vector[N1+1] beta1; //coefficients for b1 environmental dependence
	vector[N2+1] beta2; //coefficients for b2 environmental dependence
	real<lower=1e-15> sigma; //measurement noise sd
}
model{
	for(i in 1:p){ //loop over locations
		bInf[i] ~ normal(dot_product(MInf[i,],betaInf),bInf_sd); //location-specific bInf
		b1[i] ~ normal(dot_product(M1[i,],beta1),b1_sd); //location-specific b1
		b2[i] ~ normal(dot_product(M2[i,],beta2),b2_sd); //location-specific b2
		y[(ni[i]+1):ni[i+1]] ~ normal(log(bInf[i] + (b1[i]-bInf[i])*exp(x[(ni[i]+1):ni[i+1]]/b2[i])), sigma); //likelihood of data, station by station
	}
}
