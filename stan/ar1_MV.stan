data {
	int T;
	int p;
	vector[p] y0;
	matrix[p,p] PHI;
	matrix[p,p] SIGMA;
}
parameters{}
model{}
generated quantities {
	matrix[p,T] y_hat;
	y_hat[,1] = y0;
	for (t in 2:T)
		y_hat[,t] =  multi_normal_rng(PHI*y_hat[,t-1], SIGMA);
}

