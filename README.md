# EFA Factor Retention Methods

This repository demonstrates a multi-method approach to determining the number of factors to retain in Exploratory Factor Analysis (EFA) using the LibQUAL service quality dataset.

Rather than relying on a single heuristic (e.g., eigenvalues > 1), this project compares classical, regression-based, simulation-based, and residual-based retention techniques to produce a defensible and triangulated factor decision.

---

## 🧠 Methods Implemented

### Classical / Visual
- Scree plot
- Eigenvalues > 1

### Objective Scree Criteria
- Optimal coordinates
- Acceleration factor
- `nFactors` package methods

### Regression-Based
- Zoski & Jurs regression

### Model Fit-Based
- Maximum likelihood chi-square tests across competing solutions

### Simulation-Based
- Horn’s parallel analysis
- **Revised Parallel Analysis** (custom function file included)
- **Comparative Data Method** (custom function file included)

### Residual-Based
- Velicer’s Minimum Average Partial (MAP)
- Very Simple Structure (VSS)
- Residual correlation matrix diagnostics

---

## 🛠 Packages Used

- `foreign`
- `psych`
- `nFactors`

---

## 📊 Key Takeaway

Across multiple complementary criteria, a **3-factor solution** is supported for these data, though several methods suggest alternative possibilities. This demonstrates why triangulation across retention rules is preferable to relying on a single statistical criterion.

---

