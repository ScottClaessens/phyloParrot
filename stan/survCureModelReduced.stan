data{
    int N;                // number of species
    int eventLit[N];      // observed tool use in literature? (no observed tool use in lit = 0, observed tool use in lit = 1)
    int eventVid[N];      // observed tool use in videos? (no observed tool use in vids = 0, observed tool use in vids = 1)
    int numLit[N];        // number of papers published in literature before tool use identified (or if no tool use, total number of papers published)
    int numVid[N];        // number of videos published on YouTube before tool use identified (or if no tool use, total number of videos published)
}
parameters{
    real aLit;            // rate for literature exponential distribution (survival function)
    real aVid;            // rate for video exponential distribution (survival function)
    real aP[N];           // probability that each species is NOT a tool user (logit scale)
}
model{
    // declare lambdaLit, lambdaVid, and p
    real lambdaLit;   // 1 / exp(aLit)
    real lambdaVid;   // 1 / exp(aVid)
    real p[N];        // probability that each species is not a tool user
    // priors on probability NOT tool user
    aP[1:N] ~ normal( 0 , 1 );
    // prior for survival functions
    aLit ~ normal( 0 , 1 );
    aVid ~ normal( 0 , 1 );
    // linear model on probability NOT tool user (with logit link function)
    for ( i in 1:N ) p[i] = inv_logit(aP[i]);
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
