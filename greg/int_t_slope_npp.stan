data{
	int p; //extract sample size parameter from data
	int N; //total number of observations
	int ni[p+1]; //extract linear indices for station data
	vector[N] y; //extract log F(z)
	vector[N] x; //extract log(z/100)
	vector[p] T; //extract temperature variable
	vector[p] npp; //extract npp variable
}
parameters{
	////////////////////////////////////////////////////////
	//-INTECEPTS-///////////////////////////////////////////
	////////////////////////////////////////////////////////
	real<lower=-10,upper=10> beta0[p]; //vector of p intercepts
	real<lower=-1,upper=5>   beta0slp_T; //slope of intercept relation with npp
	real< beta0int; //intercept of intercept relation
	real<lower=1E-15> beta0_sd; //noise in intercept relationship
	/////////////////////////////////////////////////////////
	//-SLOPES-///////////////////////////////////////////////
	/////////////////////////////////////////////////////////
	real<lower=-3,upper=0> beta1[p];	//vector of p slopes
	real beta1slp_npp; //slope between slope and temperature
	real beta1int; //intercept of slope relation
	real<lower=1E-15> beta1_sd; //standard deviation of slopes
	////////////////////////////////////////////////////////
	//-OBSERVATIONS-////////////////////////////////////////
	////////////////////////////////////////////////////////
	real<lower=1E-15> sigma; //measurement noise sd
}
model{
	for(i in 1:p){ //loop over locations
		beta0[i] ~ normal(beta0int + beta0slp_T*T[i],     beta0_sd); //location-specific intercept
		beta1[i] ~ normal(beta1int + beta1slp_npp*npp[i], beta1_sd); //location-specific slope is a function of temperature at the location
		target += normal_lpdf(y[(ni[i]+1):ni[i+1]] | beta0[i] + beta1[i]*x[(ni[i]+1):ni[i+1]], sigma);
	}
}
generated quantities{
//	matrix[N,p] y_pred; //allocate space for prediction of y
	vector[N] log_lik;

//	for(i in 1:p){ //loop over locations
//		y_pred[,i] = beta0[i] + beta1*x; //location level prediction
//	}
	for(i in 1:p){
		for(j in (ni[i]+1):ni[i+1]){
			log_lik[j] = normal_lpdf(y[j] | beta0[i] + beta1[i]*x[j], sigma);
		}
	}
}

