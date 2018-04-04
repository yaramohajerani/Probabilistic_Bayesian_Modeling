
T <- 1000
k <- sin(seq(0,2*pi,length.out=T)) - 1
dz <- 0.01

F <- numeric(T)
F[1] <- 100

w <- 0.2

for(i in 2:(T-1)){
	F[i] <- F[i-1] + (w/(2*dz))*(F[i-1] - F[i+1]) 
}



for(i in 2:(T-1)){
	F[i] <- F[i-1] + ((w/dz)*(F[i-1] - F[i]) + k[i]*F[i-1])*dz
}

plot(F)



z <- matrix(NA,10,100)

z[1,1] <- 1
z[2:10,1] <- 0

for(i in 1:100){
	
	




