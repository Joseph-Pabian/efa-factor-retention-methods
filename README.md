# Factor Retention in Exploratory Factor Analysis  
### A Multi-Method Psychometric Decision Framework

## Overview

Determining the number of factors to retain is one of the most consequential decisions in Exploratory Factor Analysis (EFA). Over- or under-extraction can distort dimensional interpretation, bias parameter estimates, and compromise construct validity.

This project demonstrates a **triangulated, multi-method approach to factor retention** using the LibQUAL service quality dataset. Rather than relying on a single heuristic (e.g., Kaiser’s eigenvalue > 1 rule), multiple classical, regression-based, simulation-based, residual-based, and model fit criteria are compared to arrive at a defensible factor solution.

The goal is to illustrate principled statistical decision-making under methodological disagreement — a critical skill in psychometrics and quantitative research.

---

## Research Question

How many latent dimensions best represent the LibQUAL service quality items when evaluated across complementary retention criteria?

---

## Methods Implemented

### 1. Classical / Visual Methods
- Scree Plot
- Kaiser Criterion (Eigenvalues > 1)
- Objective Scree Criteria:
  - Optimal Coordinates
  - Acceleration Factor
  - `nFactors` package procedures
  ## Scree Plot

![Scree Plot](images/Screeplot.png)


### 2. Regression-Based Criterion
- Zoski & Jurs Regression Approach

### 3. Model Fit-Based Evaluation
- Maximum Likelihood Chi-Square Tests across competing factor solutions

### 4. Simulation-Based Methods
- Horn’s Parallel Analysis
- Revised Parallel Analysis (custom implementation included)
- Comparative Data Method (custom implementation included)

### 5. Residual-Based Methods
- Velicer’s Minimum Average Partial (MAP)
- Very Simple Structure (VSS)
- Residual Correlation Matrix Diagnostics

---

## Results Summary

Across multiple complementary criteria:

## Findings Overview

![findings overview](images/findings_overview.png)


