library(spatialSEIR)

compMatDim = c(200,50)
xDim = c(200*50, 5)

S0 = -1
E0 = -1
I0 = -1
R0 = -1
Sstar0 = -1
Estar0 = -1
Istar0 = -1
Rstar0 = -1
Istar = -1
X = matrix(rnorm(200*50), ncol = 5)

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

