data {
	int N; //extract parameter N from data
	//int Npred; //number of data points to give a prediction
	vector[N] x; //extract x variable from data
	vector[N] y; //extract y variable from data
	//vector[Npred] zpred; //places where you want a prediction
}
parameters {
	real beta0; //declare beta0 as continuous unbounded probability distribution
	real beta1; //declare beta1 as continuous unbounded probability distribution
	real<lower=1E-15> sigma; //declare sigma continuous probability bounded at zero, note 
}
model {
	//beta0 ~ normal(10,100);  //prior on beta0
	//beta1 ~ normal(0,100);  //prior on beta1
	y     ~ normal(beta0 + beta1*x, sigma); //likelihood
}
//generated quantities{ //should comment this section out if you have a big dataset as it generates a ton of output
//	vector[Npred] y_pred; //allocate y_pred vector
//	y_pred = exp(beta0)*(zpred/100)^beta1;  //posterior for 'true y' given x
//}
//generated quantities{ //should comment this section out if you have a big dataset as it generates a ton of output
//	vector[N] y_pred; //allocate y_pred vector
//	y_pred = beta0 + beta1*x;  //posterior for 'true y' given x
//}

