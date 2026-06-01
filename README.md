# Computational-QCD

# Computational QCD: Computational QCD framework for extracting the running strong coupling constant from ATLAS multi-jet measurements.

## Overview

This repository contains a computational framework for the determination of the strong coupling constant ((\alpha_s)) from multi-jet measurements using data from the ATLAS experiment at CERN.

The analysis combines experimental measurements, perturbative Quantum Chromodynamics (QCD) predictions, parton distribution functions (PDFs), systematic uncertainty propagation, and statistical inference techniques to study the energy-scale dependence of (\alpha_s). The framework implements chi-square minimization, nuisance-parameter profiling, PDF uncertainty evaluation, and scale variation studies for precision QCD measurements.

---

## Physics Motivation

The strong coupling constant ((\alpha_s)) is a fundamental parameter of the Standard Model that determines the strength of the strong interaction. Quantum Chromodynamics predicts that (\alpha_s) varies with energy scale, a phenomenon known as asymptotic freedom.

Multi-jet cross-section ratios provide a powerful probe of this behavior because many experimental and theoretical uncertainties partially cancel in the ratio. This repository implements tools used to extract (\alpha_s) and investigate its running using ATLAS Run-2 measurements.

---

## Features

* Extraction of the strong coupling constant ((\alpha_s))
* Chi-square minimization using ROOT and TMinuit
* Profile likelihood fits with nuisance parameters
* MSHT20 PDF uncertainty evaluation
* CT18 PDF uncertainty support
* Renormalization and factorization scale variation studies
* Statistical and systematic uncertainty propagation
* Non-perturbative correction handling
* Theory-to-data comparison tools
* Publication-quality ROOT plotting utilities

---

## Repository Contents

### `alphaS_Running_Extraction_MSHT.C`

Main fitting framework used to determine the running strong coupling constant across multiple energy scales.

### `chi2Min.C`

Core statistical fitting engine implementing:

* Chi-square minimization
* Nuisance parameter treatment
* PDF uncertainty evaluation
* Scale uncertainty calculations
* Theory prediction generation

### `Style.C` / `Style.h`

ATLAS-style plotting configuration for ROOT.

### `Utils.h`

Utility functions for graph manipulation, uncertainty bands, ratio calculations, and ROOT-based visualization.

---

## Software Requirements

* ROOT
* C++
* TMinuit
* MSHT20 PDF parameterizations
* Input ROOT files containing measured and theoretical distributions

---

## Scientific Applications

## Analysis Workflow

```text
ATLAS Multi-Jet Data
        |
        v
HEPData Input
        |
        v
Experimental Covariance
        |
        v
Theory Predictions + PDF/Scale Variations
        |
        v
χ² Minimization with TMinuit/MIGRAD
        |
        v
αs Extraction
        |
        v
Running of αs(Q)



## Acknowledgments and Disclaimer

This repository contains analysis software developed in the context of research performed within the ATLAS Collaboration. The framework includes original code developed by the author together with components inspired by, adapted from, or built upon software tools and methodologies commonly used within the high-energy physics community.

Experimental measurements, detector information, and published physics results referenced by this work remain the property of the ATLAS Collaboration and their respective authors. This repository is intended for scientific, educational, and research purposes only.

The software and views presented here are those of the author and do not constitute an official ATLAS Collaboration software release, recommendation, or publication.

---

## Author

**Zahra Farazpay**

Ph.D. in Physics

Computational Physics | Quantum Chromodynamics | Scientific Computing | Statistical Inference | Machine Learning

Former ATLAS Collaboration Author

