data {
	int T;
	real y0;
	real phi;
	real sigma;
}
parameters{}
model {}
generated quantities {
	vector[T] y_hat;
	y_hat[1] = y0;
	for (t in 2:T)
		y_hat[t] =  normal_rng(phi*y_hat[t-1], sigma);
}

