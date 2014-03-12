library(spatialSEIR)

compMatDim = c(20,50)
xDim = c(20*50, 5)

S0 = -1
E0 = -1
I0 = -1
R0 = -1
Sstar0 = -1
Estar0 = -1
Istar0 = -1
Rstar0 = -1
Istar = as.numeric(matrix(rbinom(20*50, 500,0.1), ncol = 50))
X = as.numeric(matrix(rnorm(20*50*5), nrow =20*50 , ncol = 5))


# Doesn't do anything interesting, just testing data transfer
res = spatialSEIRInit(compMatDim,
                      xDim,
                      S0,
                      E0,
                      I0,
                      R0,
                      Sstar0,
                      Estar0,
                      Istar0,
                      Rstar0,
                      Istar,
                      X)

