# Hierarchical Functional Models

Hierarchical Functional Models, unified sampler with three types of priors
- Matrix Normal (Naive Bayesian Conjugate Prior)
- Strongly Separable Latent Factor Prior
- Non-Separable Latent Factor Prior

# Installation:
```r
devtools::install_github("Qian-Li/HFM")
```

# Updates (March 6th, 2018):
- [x] Unified implementation across all priors;
- [x] Modified separable simulation: simulated B-spline coefficients directly instead of GP observations;
- [x] Second-moment assessment (Random effects)
- [x] Simulation 1: Random missing observations on varying bs.df (num of spline coefs); Sep and Non-sep 500
- [x] Simulation 2: Random missing observations on varying spatial dependency (rho = .4/.6/.8); Sep and Non-sep 500
- [x] Simulation 3: Varying Sample size and SNR; (100);
- [x] Simulation 4: Robustness of noLF when true noLF unknown; (100);

# Problems: 
- [ ] Crash occasionally on `pinv()` or `chol()`;