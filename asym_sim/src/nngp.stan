/* kappa_p_1 is (tausq / sigmasq) + 1 */

functions{
real nngp_lpdf(vector Y, matrix X, vector beta, real sigmasq, real tausq,
                real phi, matrix neardist, matrix neardistM, int[,] nearind,
                int N, int M){

        vector[N] V;
        vector[N] Uw;
        vector[N] temp_w;
        int dim;
        int h;
        real kappa_p_1;
        real out;
        kappa_p_1 = tausq / sigmasq + 1;
        Uw = Y - X * beta;
        temp_w = Uw;

        for (i in 2:N) {

            matrix[ i < (M + 1) ? (i - 1) : M, i < (M + 1) ? (i - 1): M]
            temp_neardistM;
            matrix[ i < (M + 1) ? (i - 1) : M, i < (M + 1) ? (i - 1): M] L;
            vector[ i < (M + 1) ? (i - 1) : M] u;
            vector[ i < (M + 1) ? (i - 1) : M] v;
            row_vector[i < (M + 1) ? (i - 1) : M] v2;

            dim = (i < (M + 1))? (i - 1) : M;

            // get exp(-phi * neardistM) + (tausq / sigmasq) * I

            if(dim == 1){temp_neardistM[1, 1] = kappa_p_1;}
            else{
                h = 0;
                for (j in 1:(dim - 1)){
                    for (k in (j + 1):dim){
                        h = h + 1;
                        temp_neardistM[j, k]
                            = exp(- phi * neardistM[(i - 1), h]);
                        temp_neardistM[k, j] = temp_neardistM[j, k];
                    }
                }
                for(j in 1:dim){
                    temp_neardistM[j, j] = kappa_p_1;
                }
            }

            L = cholesky_decompose(temp_neardistM);

            for (j in 1: dim){
                u[j] = exp(- phi * neardist[(i - 1), j]);
            }
            
            // v = L^-1 * u
            v = mdivide_left_tri_low(L, u);

            V[i] = kappa_p_1 - (v' * v);

            // v2 = L^-T L^-1 * u
            v2 = mdivide_right_tri_low(v', L);

            for (j in 1:dim){
                Uw[i] = Uw[i] - v2[j] * temp_w[nearind[(i - 1), j]];
            }
        }

        V[1] = (kappa_p_1);
        out = 0.0;
        for (i in 1:N){
            out = out - 0.5 * log(V[i]) - 0.5 / sigmasq * (Uw[i] * Uw[i] / V[i]);
        }
        out = out - 0.5 * N * log(sigmasq);
        return out;
    }
}

data {
    int<lower=1> N;
    int<lower=1> M;
    int<lower=1> P;
    vector[N] Y;
    matrix[N, P] X;
    int nearind[N - 1, M];
    matrix[N - 1, M] neardist;
    matrix[N - 1, (M * (M - 1) / 2)] neardistM;
    real as;
    real bs;
    real at;
    real bt;
    real ap;
    real bp;
}


parameters{
    vector[P] beta;
    real sigmasq;
    real tausq;
    real<lower = ap, upper = bp> phi;
}

model{
    sigmasq ~ inv_gamma(as, bs);
    tausq ~ inv_gamma(at, bt);
    Y ~ nngp(X, beta, sigmasq, tausq, phi, neardist, neardistM, nearind, N, M);
}

/*
generated quantities {
    vector[N] log_lik;
    for (i in 1:N){
        log_lik[i] <- normal_lpdf(Y[i] |beta[1] + X[i, 2] * beta[2], sqrt(tau^2 + sigma^2));
    }
}
*/













