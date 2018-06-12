library('Matrix')
library("fields")
library(rstan)

## stan code for A ##
model_code <-
  "
    functions {
        matrix getADstan(matrix neardist, matrix neardistM, int N, int M, real phi) {
            int dim;
            int h;
            matrix[M + 1, N] AD;
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
                
                AD[(M + 1), i] = (1.0 - (v' * v));
                
                v2 = mdivide_right_tri_low(v', L);
                
                for(j in 1:dim){
                    AD[j, i] = v2[j];
                }
            }
            AD[(M + 1), 1] = 1;
            return AD;
        }
        
        matrix getADstan_ho(matrix neardist, matrix neardistM, int N, int M, real phi,
        real deltasq) {
            int h;
            matrix[M + 1, N] AD;
            for (i in 1:N) {
                matrix[M, M] temp_neardistM;
                matrix[M, M] L;
                vector[M] u; vector[M] v;
                row_vector[M] v2;
                
                // get exp(-phi * neardistM)
                h = 0;
                for (j in 1:(M - 1)){
                    for (k in (j + 1):M){
                        h = h + 1;
                        temp_neardistM[j, k] = exp(- phi * neardistM[h, i]);
                        temp_neardistM[k, j] = temp_neardistM[j, k];
                    }
                }
                for(j in 1:M){
                    temp_neardistM[j, j] = 1;
                }
                
                L = cholesky_decompose(temp_neardistM);
                
                for (j in 1: M){
                    u[j] = exp(- phi * neardist[i, j]);
                }
                
                //vector[dim] v;
                v = mdivide_left_tri_low(L, u);
                
                AD[(M + 1), i] = (1.0 + deltasq - (v' * v));
                
                v2 = mdivide_right_tri_low(v', L);
                
                for(j in 1: M){
                    AD[j, i] = v2[j];
                }
            }
            return AD;
        }
        
        matrix getADstan2(matrix neardist, matrix neardistM, int N, int M, real phi,
                          real deltasq) {
            int dim;
            int h;
            matrix[M + 1, N] AD;
            for (i in 2:N) {
                matrix[ i < (M + 1)? (i - 1) : M, i < (M + 1)? (i - 1): M] temp_neardistM;
                matrix[ i < (M + 1)? (i - 1) : M, i < (M + 1)? (i - 1): M] L;
                vector[ i < (M + 1)? (i - 1) : M] u;
                vector[ i < (M + 1)? (i - 1) : M] v;
                row_vector[i < (M + 1)? (i - 1) : M] v2;
                
                dim = (i < (M + 1))? (i - 1) : M;
                
                // get exp(-phi * neardistM)
                if(dim == 1){temp_neardistM[1, 1] = 1.0 + deltasq;}
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
                        temp_neardistM[j, j] = 1.0 + deltasq;
                    }
                }
                
                L = cholesky_decompose(temp_neardistM);
                
                for (j in 1: dim){
                    u[j] = exp(- phi * neardist[(i - 1), j]);
                }
                
                //vector[dim] v;
                v = mdivide_left_tri_low(L, u);
                
                AD[(M + 1), i] = (1.0 + deltasq - (v' * v));
                
                v2 = mdivide_right_tri_low(v', L);
                
                for(j in 1:dim){
                    AD[j, i] = v2[j];
                }
            }
            AD[(M + 1), 1] = 1 + deltasq;
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

