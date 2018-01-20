/* Latent NNGP */

functions{
real nngp_w_lpdf(vector w_r, real sigmasq, real phi, matrix neardist, matrix neardistM,
                int[,] nearind, int N, int M, real beta1){

        vector[N] V;
        vector[N] Uw;
        vector[N] w;
        int dim;
        int h;
        real out;

        for (i in 1:N){
            w[i] = w_r[i] - beta1;
        }
        Uw = w;

        for (i in 2:N) {

            matrix[ i < (M + 1)? (i - 1) : M, i < (M + 1)? (i - 1): M] temp_neardistM;
            matrix[ i < (M + 1)? (i - 1) : M, i < (M + 1)? (i - 1): M] L;
            vector[ i < (M + 1)? (i - 1) : M] u;
            vector[ i < (M + 1)? (i - 1) : M] v;
            row_vector[i < (M + 1)? (i - 1) : M] v2;

            dim = (i < (M + 1))? (i-1) : M;

            if(dim == 1){temp_neardistM[1, 1] = 1;}
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
                    temp_neardistM[j, j] = 1;
                }
            }

            L = cholesky_decompose(temp_neardistM);

            for (j in 1: dim){
                u[j] = exp(- phi * neardist[(i - 1), j]);
            }

            v = mdivide_left_tri_low(L, u);

            V[i] = 1 - (v' * v);

            v2 = mdivide_right_tri_low(v', L);

            for (j in 1:dim){
                Uw[i] = Uw[i] - v2[j] * w[nearind[(i - 1), j]];
            }
        }
        V[1] = 1;
        out = 0.0;
        for (i in 1:N){
            out = out - 0.5 * log(V[i]) - 0.5 / sigmasq * (Uw[i] * Uw[i] / V[i]);
        }
        out = out - 0.5 * N * log(sigmasq);
        return out;
    }
}


data {
    int N;
    int M;
    int P;
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
    vector[N] w;
}

model{
    sigmasq ~ inv_gamma(as, bs);
    tausq ~ inv_gamma(at, bt);
    w ~ nngp_w(sigmasq, phi, neardist, neardistM, nearind, N, M, beta[1]);
    Y ~ normal(X[, 2] * beta[2] + w, sqrt(tausq));
}















