### the Euclidean norm

EuclideanNorm <- function(x) {sqrt(colSums(x^2))}


### huberizing a vector to length b

Huberize <- function(x, b, norm = EuclideanNorm, ...)
{  nx <- norm(x, ...)
   Ind0 <- (nx < b  )
   Ind1 <- (nx < b/2)
   x*(Ind0 + (1-Ind0)* b / (nx + Ind1) )
}    

limitS <- function(S, F, Z, Q, V, tol = 10^-4, itmax = 1000)#
## determines lim_{t->infty} S_{t|t-1}
     {SO0 <- S + 1
      S0  <- S
      i   <- 0
      while( (sum( (SO0 - S0)^2 ) > tol^2) && (i < itmax) )
        {i   <- i + 1
         S1  <- .getpredCov(S0, F, Q)
         SO0 <- S0
         K   <- .getKG(S1, Z, V)
         S0  <- .getcorrCov(S1, K, Z)
        }
     S1
     }

rootMatrix <- function (X)
{
###########################################
##
##  R-function: rootMatrix - computes the unique square root 'A' 
##                           of matrix 'X', i.e., A%*%A = X
##              former R-function 'root.matrix' of package 'strucchange' 
##  author: Bernhard Spangl, based on work of Achim Zeileis
##  version: 0.2 (2008-02-24)
##
###########################################

##  Paramters:
##  X ... symmetric and positive semidefinite matrix 

    if ((ncol(X) == 1) && (nrow(X) == 1)) 
        return(list(X.det=X, 
                    X.sqrt=matrix(sqrt(X)), X.sqrt.inv=matrix(1/sqrt(X))))
    else {
        X.eigen <- eigen(X, symmetric = TRUE)
        if (any(X.eigen$values < 0)) 
            stop("matrix is not positive semidefinite")
        sqomega <- sqrt(diag(X.eigen$values))
        sqomega.inv <- diag(1/sqrt(X.eigen$values))
        V <- X.eigen$vectors
        X.sqrt <- V %*% sqomega %*% t(V)
        X.sqrt.inv <- V %*% sqomega.inv %*% t(V)
        return(list(X.det=prod(X.eigen$values), 
                    X.sqrt=X.sqrt, X.sqrt.inv=X.sqrt.inv))
    }
}

