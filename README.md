# Hierarchical Functional Models

Hierarchical Functional Models, unified sampler with three types of priors
- Matrix Normal (Naive Bayesian Conjugate Prior) (`hfm.NB`)
- Strongly Separable Latent Factor Prior (`hfm.LF` with two inputs)
- Non-Separable Latent Factor Prior (`hfm.LF` with one input)

# Installation:
```r
# For Mac Users
devtools::install_github("Qian-Li/HFM")
# For Windows and Linux Users, please try install the binary version
```

## Updates (March 6th, 2018):
- [x] Unified implementation across all priors;
- [x] Modified separable simulation: simulated B-spline coefficients directly instead of GP observations;
- [x] Second-moment assessment (Random effects)
- [x] Simulation 1: Random missing observations on varying bs.df (num of spline coefs); Sep and Non-sep 500
- [x] Simulation 2: Random missing observations on varying spatial dependency (rho = .4/.6/.8); Sep and Non-sep 500
- [x] Simulation 3: Varying Sample size and SNR; (100);
- [x] Simulation 4: Robustness of noLF when true noLF unknown; (100);

## Updates (March 29th, 2018):
- [x] Data Analysis: Alpha and Gamma band regressed on Groups;
- [x] Data Analysis: Alpha and Gamma band regressed on Groups and Age;
- [x] Data Analysis: Alpha and Gamma band regressed on Age and VDQ;

## Updates (May 3rd, 2018):
- [x] 1.0 Published: all issues fixed and sampling progrssion bar added :-)
- [ ] Coming soon: A vignette and short working example of how to use the package

## Debug: 
- [x] Crash occasionally on `pinv()` or `chol()`;