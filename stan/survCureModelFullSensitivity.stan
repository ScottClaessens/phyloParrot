// generated on 31.03.2022 using the rethinking R package (version 2.13)
functions{
    // ornstein-uhlenbeck covariance kernel function
    // for phylogenetic gaussian process covariance
    matrix cov_GPL1(matrix x, real sq_alpha, real sq_rho, real delta) {
        int N = dims(x)[1];
        matrix[N, N] K;
        for (i in 1:(N-1)) {
          K[i, i] = sq_alpha + delta;
          for (j in (i + 1):N) {
            // covariance = eta^2 * exp(-rho^2 * distance)
            K[i, j] = sq_alpha * exp(-sq_rho * x[i,j] );
            K[j, i] = K[i, j];
          }
        }
        K[N, N] = sq_alpha + delta;
        return K;
    }
}
data{
    int N;                // number of species
    int eventLit[N];      // observed tool use in literature? (no observed tool use in lit = 0, observed tool use in lit = 1)
    int eventVid[N];      // observed tool use in videos? (no observed tool use in vids = 0, observed tool use in vids = 1)
    int numLit[N];        // number of papers published in literature before tool use identified (or if no tool use, total number of papers published)
    int numVid[N];        // number of videos published on YouTube before tool use identified (or if no tool use, total number of videos published)
    int Species[N];       // species id
    real EQ[N];           // encephalisation quotient (standardised)
    int Feeding[N];       // feeding type (generalist = 1, specialist = 2)
    matrix[N,N] distMat;  // scaled phylogenetic pairwise distance matrix
}
parameters{
    vector[N] z;          // standardised normal parameters for non-centered gaussian process 
    real aLit;            // rate for literature exponential distribution (survival function)
    real aVid;            // rate for video exponential distribution (survival function)
    vector[2] aP;         // vector of intercepts for the probability NOT tool user (one for generalists, one for specialists)
    real bEQ;             // coefficient for EQ predicting the probability NOT tool user
    real<lower=0> etasq;  // maximum covariance for gaussian process
    real<lower=0> rhosq;  // rate for gaussian process
}
transformed parameters{
    vector[N] k;          // varying species-specific intercepts for probability NOT tool user
    matrix[N,N] L_SIGMA;  // cholesky decomposition of covariance matrix
    matrix[N,N] SIGMA;    // covariance matrix
    // non-centered parameterisation
    SIGMA = cov_GPL1(distMat, etasq, rhosq, 0.01);
    L_SIGMA = cholesky_decompose(SIGMA);
    k = L_SIGMA * z;
}
model{
    // declare lambdaLit, lambdaVid, and p
    real lambdaLit;   // 1 / exp(aLit)
    real lambdaVid;   // 1 / exp(aVid)
    vector[N] p;      // probability NOT tool user, for each species
    // priors for gaussian process
    rhosq ~ exponential( 0.5 );
    etasq ~ exponential( 0.5 );
    // priors for linear model on probability NOT tool user
    bEQ ~ normal( 0 , 1 );
    ////// modify intercept prior from normal(0, 1) to normal(1.78507, 2)
    ////// which is 0.14 on the probability scale (actual proportion of observed tool users)
    aP ~ normal( 1.78507 , 2 );
    // prior for survival functions
    aLit ~ normal( 0 , 1 );
    aVid ~ normal( 0 , 1 );
    // standardised normal priors for non-centered parameterisation
    z ~ normal( 0 , 1 );
    // linear model on probability NOT tool user (with logit link function)
    for ( i in 1:N ) {
        p[i] = aP[Feeding[i]] + k[Species[i]] + bEQ * EQ[i];
        p[i] = inv_logit(p[i]);
    }
    // survival functions
    lambdaLit = 1 / exp(aLit);
    lambdaVid = 1 / exp(aVid);
    // survival cure likelihood (literature)
    // if no observed tool use in lit/vids, then species is either not a tool user (p)
    // or a tool user but hasn't been observed in lit/vids yet (1 - p)
    // else if observed tool use in lit/vids, then species is a tool user (1 - p) that has been observed in lit/vids
    // Pr(y | event == 0) = p + (1 - p) * ExponentialCCDF( lambda )
    // Pr(y | event == 1) = (1 - p) * Exponential( lambda )
    for ( i in 1:N )
        if ( eventLit[i] == 0 ) target += log_mix(p[i], 0, exponential_lccdf(numLit[i] | lambdaLit));
    for ( i in 1:N )
        if ( eventLit[i] == 1 ) target += log1m(p[i]) + exponential_lpdf(numLit[i] | lambdaLit);
    for ( i in 1:N )
        if ( eventVid[i] == 0 ) target += log_mix(p[i], 0, exponential_lccdf(numVid[i] | lambdaVid));
    for ( i in 1:N )
        if ( eventVid[i] == 1 ) target += log1m(p[i]) + exponential_lpdf(numVid[i] | lambdaVid);
}
