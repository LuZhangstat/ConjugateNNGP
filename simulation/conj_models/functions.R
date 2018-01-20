library('Matrix')
library("fields")
library(rstan)

## stan code for A ##
model_code <-
  "
    functions {
        matrix getADstan(matrix neardist, matrix neardistM,
        matrix nearind, int N, int M, real phi) {
        int dim;
        int h;
        matrix[N, M + 1] AD;
        for (i in 2:N) {
            matrix[ i < (M + 1)? (i - 1) : M, i < (M + 1)? (i - 1): M] temp_neardistM;
            matrix[ i < (M + 1)? (i - 1) : M, i < (M + 1)? (i - 1): M] L;
            vector[ i < (M + 1)? (i - 1) : M] u;
            vector[ i < (M + 1)? (i - 1) : M] v;
            row_vector[i < (M + 1)? (i - 1) : M] v2;

            dim = (i < (M + 1))? (i-1) : M;

            // get exp(-phi * neardistM)
            if(dim == 1){temp_neardistM[1, 1] = 1;}
            else{
                h = 0;
                for (j in 1:(dim - 1)){
                    for (k in (j + 1):dim){
                        h = h + 1;
                        temp_neardistM[j, k] = exp(- phi * neardistM[(i - 1), h]);
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

            //vector[dim] v;
            v = mdivide_left_tri_low(L, u);

            AD[i, (M+1)] = (1.0 - (v' * v));

            v2 = mdivide_right_tri_low(v', L);

            for(j in 1:dim){
                AD[i, j] = v2[j];
            }
        }
        AD[1, (M+1)] = 1;
        return AD;
        }
    }

    model {}
"
expose_stan_functions(stanc(model_code = model_code))

### conjugate gradient by RcppEigen ###
library(Rcpp)
library(RcppEigen)
library(inline)

incl<- "
using Eigen::Map;
using Eigen::VectorXd;
using Eigen::MappedSparseMatrix;
using Eigen::SparseMatrix;
using Rcpp::as;
using Eigen::ColMajor;
typedef MappedSparseMatrix<double, ColMajor, int> MapSpMatd;
typedef Map<VectorXd> MapVecd;
"


cgSparsecode <- "
#include <RcppEigen.h>
#include <Rcpp.h>

MapSpMatd    A(as<MapSpMatd >(AA));
const MapVecd      b(as<MapVecd>(bb));

Eigen::ConjugateGradient< SparseMatrix<double>, Eigen::Lower > solver;
return wrap(solver.compute(A).solve(b));

"

cgsparse <- cxxfunction(signature(AA = "dgCMatrix", bb = "numeric"),
                        cgSparsecode, "RcppEigen", incl)

