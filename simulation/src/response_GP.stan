/* Full Gaussian Process without w */

functions{
    matrix fullGP_L(real sigmasq, real tausq, real phi, vector dist, int N){

    int h;
    matrix[N, N] temp_dist;
    matrix[N, N] L;

    h = 0;
    for (j in 1:(N - 1)){
        for (k in (j + 1):N){
            h = h + 1;
            temp_dist[j, k] = sigmasq * exp(- phi * dist[h]);
            temp_dist[k, j] = temp_dist[j, k];
        }
    }

    for(i in 1:N){
        temp_dist[i,i] = sigmasq + tausq;
    }

    L = cholesky_decompose(temp_dist);

    return L;
    }
}

data {
    int<lower=1> N;
    int<lower=1> P;
    vector[N] Y;
    matrix[N, P] X;
    vector[N * (N - 1) / 2] dist;
    real as;
    real bs;
    real at;
    real bt;
    real ap;
    real bp;
}

parameters {
    vector[2] beta;
    real sigmasq;
    real tausq;
    real<lower = ap, upper = bp> phi;
}


model{
    sigmasq ~ inv_gamma(as, bs);
    tausq ~ inv_gamma(at, bt);
    Y ~ multi_normal_cholesky(X * beta, fullGP_L(sigmasq, tausq, phi, dist, N));
}















