using Stan, StanSample, Distributions
bmds_model = """data {
    int<lower=1> N; // number of nodes
    int<lower=0, upper=(N*(N-1))/2> E; // number of edges
    int<lower=1, upper=N-1> D; // target dimension
    int<lower=1, upper=N> edges[E,2]; // a list of edges (i1, i2)
    vector<lower=0>[E] distances;
    int<lower=0, upper=3> censoring[E]; // censoring of distances
    /* 0: uncensored, 
     * 1: left-censored, 
     * 2: right-censored, 
     * 3: missing
     */
}

parameters {
    /* the first point is fixed at the origin (0)
     * the next D points are stored in a lower triangular matrix (X1)
     * the remaining points are given by a N-D-1 x D matrix (X2)
     */
    cholesky_factor_cov[D] X1;
    matrix[D, N-D-1] X2;
    real<lower=0> sigma; // nuisance parameter
}

transformed parameters {
    // put all coordinates in a single matrix X
    matrix[D, N] X; // [0, X1', X2]
    X = append_col(rep_vector(0, D), append_col(X1', X2));
}

model {
    // try to keep sigma small
    sigma ~ exponential(1);
    // calculate the current distances (between the Xs)
    for ( e in 1:E ) {
        int i1; int i2;
        real dist;
        
        // compute the distance
        i1 = edges[e][1];
        i2 = edges[e][2];
        dist = distance(col(X, i1), col(X, i2));

        // likelihood
        if ( censoring[e] == 0 ) { // uncensored
            distances[e] ~ normal(dist, sigma);
        } else if ( censoring[e] == 1 ) { // left-censored
            target += normal_lcdf(distances[e] | dist, sigma);
        } else if ( censoring[e] == 2 ) { // right-censored
            target += normal_lccdf(distances[e] | dist, sigma);
        } else if ( censoring[e] == 3 ) { // missing
            // do nothing
        }
    }
}

generated quantities {
    vector[D] meanX; // the means of each MCMC sample
    matrix[D, N] Xc; // centered positions
    matrix[D, N] Xcr; // centered and rotated positions

    // center X
    for ( i in 1:D ) {
        meanX[i] = mean(row(X, i));
    }
    Xc = X - rep_matrix(meanX, N);
    // rotate Xc
    Xcr = eigenvectors_sym(tcrossprod(Xc))' * Xc;
}"""

m6_1s = SampleModel("m6.1s", bmds_model);