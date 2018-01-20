/* Latent GP */

functions{
matrix get_w_L(real sigmasq, real phi, vector dist, int N){

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
        temp_dist[i, i] = sigmasq;
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
    vector[P] beta;
    real sigmasq;
    real tausq;
    real<lower = ap, upper = bp> phi;
    vector[N] w;
}


model{
    vector[N] mu_w;
    sigmasq ~ inv_gamma(as, bs);
    tausq ~ inv_gamma(at, bt);

    for(i in 1:N){
        mu_w[i] = beta[1];
    }

    w ~ multi_normal_cholesky(mu_w, get_w_L(sigmasq, phi, dist, N));
    Y ~ normal(X[, 2] * beta[2] + w, sqrt(tausq));
}











