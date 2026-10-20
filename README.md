# Discriminating between the Weibull and log-normal distributions with a test based on the third log-cumulant

Repository for the paper *"Discriminating between the Weibull and log-normal distributions with a test based on the third log-cumulant"*, by **Caio C. F. Santos** (Universidade Federal de Pernambuco / CASTLab) and **Raydonal Ospina** (Universidade Federal da Bahia / UFPE).

## Overview

The Weibull and log-normal lifetime models are widely used and often fit the same data well, which makes the choice between them hard in small samples.

We propose a test of the log-normal against the Weibull based on the **third log-cumulant**, which is

- **zero** for every log-normal distribution, and
- **−2ζ(3)/β³ < 0** for every Weibull distribution with shape parameter β.

The sign of this log-cumulant has a geometric reading in the log-cumulant diagram. The statistic standardises the bias-corrected sample third log-cumulant and, by invariance under the group of transformations x ↦ b·xᶜ, has a parameter-free distribution under both families.

## Highlights

- **Size:** with the bias correction, the standard normal critical value preserves the chosen significance level (α = 5%) in Monte Carlo experiments.
- **Power:** reaches 0.95 at n ≈ 86.
- **Comparison:** benchmarked against the ratio of maximised likelihoods (RML) rule and Kappenman's (1988) rule in terms of probability of correct selection (PCS). Under the log-normal the test errs with probability close to α for every n; under the Weibull its PCS exceeds 0.90 from n = 69 onwards.
- **Applications:** three classical data sets — the test does not reject the log-normal in two and rejects it in the third, in agreement with the other two rules.

## Keywords

Mellin transform · log-cumulants · skewness test · model discrimination · invariant tests · reliability

## Repository contents

```
.
├── README.md
├── paper/            # LaTeX source and compiled PDF of the article
├── code/             # Test statistic, bias correction, Monte Carlo experiments
├── data/             # Classical life data sets used in the applications
└── figures/          # Log-cumulant diagram, power and PCS plots
```

## How to cite

```bibtex
@article{santos2026weibulllognormal,
  title   = {Discriminating between the Weibull and log-normal distributions
             with a test based on the third log-cumulant},
  author  = {Santos, Caio C. F. and Ospina, Raydonal},
  year    = {2026},
  note    = {Preprint}
}
```

## Authors

- **Caio C. F. Santos** — [ORCID 0009-0003-6410-8854](https://orcid.org/0009-0003-6410-8854) — contact: [caio.franca@ufpe.br](mailto:caio.franca@ufpe.br)
- **Raydonal Ospina** — [ORCID 0000-0002-9884-9090](https://orcid.org/0000-0002-9884-9090)

## License

To be chosen by the authors (e.g., MIT for the code; CC BY 4.0 for the paper).

