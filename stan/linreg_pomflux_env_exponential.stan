data{
	int p; //extract sample size parameter from data
	int N; //total number of observations
	int NJ0; //number of environmental variables for  bInf
	int Nf; //number of environmental variables for b1
	int Nl; //number of environmental variables for b2
	int ni[p+1]; //linear indices for station data
	vector[N] y; //y variable
	vector[N] x; //x variable
	matrix[p,NJ0+1] MJ0; //environmental variables for bInf
	matrix[p,Nf+1] Mf; //environmental variables for b1
	matrix[p,Nl+1] Ml; //environmental variables for b2
}
parameters{
	real<lower=1e-15> J0[p]; //vector of p intercepts (at infinity)
	real<lower=0,upper=1> f[p];
	real l[p];	//vector of len p for denominator in exponential
	real<lower=1e-15,upper=10> J0_sd; //bounded variable for intercept bInf
	real<lower=1e-15,upper=10> f_sd; //standard deviation of b1
	real<lower=1e-15,upper=10> l_sd; //standard deviation of b2
	vector[NJ0+1] betaJ0; //coefficients for bInf environmental dependence
	vector[Nf+1] betaf; //coefficients for b1 environmental dependence
	vector[Nl+1] betal; //coefficients for b2 environmental dependence
	real<lower=1e-15> sigma; //measurement noise sd
	real<lower=1e-15> alpha;
}
model{
	for(i in 1:p){ //loop over locations
		J0[i] ~ lognormal(dot_product(MJ0[i,],betaJ0),J0_sd); //location-specific bInf
		f[i] ~ normal(dot_product(Mf[i,],betaf),f_sd); //location-specific bInf
		l[i] ~ normal(dot_product(Ml[i,],betal),l_sd); //location-specific b2
		y[(ni[i]+1):ni[i+1]] ~ normal(log(J0[i]*(f[i]+(1-f[i])*exp(x[(ni[i]+1):ni[i+1]]/l[i]))), sigma); //likelihood of data, station by station
	}
}
