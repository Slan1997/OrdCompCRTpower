# OrdCompCRTpower

**OrdCompCRTpower** is an interactive R Shiny application for power and sample size calculations for clinical trials with **ordinal and ordinal composite outcomes**, with a primary focus on cluster randomized trials (CRTs).

The application provides a graphical interface for specifying trial-design assumptions and automatically performs the corresponding power and sample size calculations. It is intended to make methods for ordinal trial design more accessible to statisticians, clinical trialists, and applied researchers without requiring users to directly implement the underlying formulas in R.

## Features

OrdCompCRTpower integrates tools for:

- deriving ordinal composite endpoint distributions from correlated binary components;
- power and sample size calculation for individually randomized trials (IRTs) with ordinal outcomes;
- power and sample size calculation for CRTs with three-category ordinal outcomes;
- design-stage calculations for CRTs with ordinal outcomes containing four or more categories;
- comparison of GEE-Ordinal, Whitehead design-effect, and GEE-Binary approaches;
- sensitivity analyses across outcome distributions, treatment effects, cluster sizes, numbers of clusters, and within-cluster correlation assumptions; and
- interactive visualization and downloadable power results.

## Application Modules

The application contains four main modules.

### 1. Composite

The **Composite** module derives the distribution of an ordinal composite endpoint from its underlying binary components.

Users specify:

- marginal prevalences of the binary components;
- the rule used to construct the ordinal composite endpoint; and
- assumptions about dependence among the binary components.

The dependence structure is modeled using vine copulas. The module supports:

- **Rank by severity**, where binary components are ordered according to clinical severity;
- **Summation**, where the ordinal outcome is defined by the number of positive components; and
- a built-in **COVID-19 ordinal scale example** illustrating a more complex composition rule.

For user-defined rank-by-severity and summation endpoints, dependence can be specified using either simplified or fully manual inputs.

The module derives the resulting ordinal endpoint distribution and provides a copyable version of the distribution that can be transferred directly into the IRT or CRT modules for subsequent power and sample size calculations.

### 2. IRT

The **IRT** module performs power and sample size calculations for individually randomized trials with ordinal outcomes using **Whitehead's method**.

Users can specify:

- the ordinal outcome distribution;
- treatment effect;
- Type I error rate;
- target power; and
- range of total sample sizes.

The module displays:

- the resolved treatment and control outcome distributions;
- an interactive power curve;
- the minimum total sample size required to achieve the target power; and
- a downloadable power table.

The IRT module supports ordinal outcomes with an arbitrary number of ordered categories.

### 3. CRT (3 Cat.)

The **CRT (3 Cat.)** module performs power and sample size calculations for CRTs with three-category ordinal outcomes.

The primary method is the proposed **GEE-Ordinal** approach, which directly incorporates the matrix-valued within-cluster correlation structure of the ordinal outcome.

The module also allows comparison with:

- **Whitehead-DE**, which applies Whitehead's ordinal method with a scalar-ICC design-effect adjustment;
- **GEE-Binary**, based on dichotomizations of the ordinal outcome;
- **GEE-Independence**; and
- Whitehead's method ignoring within-cluster correlation.

Users can generate:

1. power versus the total number of clusters for a fixed cluster size;
2. power versus cluster size for a fixed total number of clusters; or
3. the required total number of clusters versus cluster size for a specified target power.

Within-cluster dependence for the ordinal outcome is specified through

```text
rho11, rho12, rho22
```

corresponding to the 2 x 2 correlation matrix for the first two category indicators.

The application includes validity checks for these correlation parameters and displays allowable ranges conditional on the other specified correlation values.

### 4. CRT (4+ Cat.)

The **CRT (4+ Cat.)** module supports CRT design when the ordinal outcome contains four or more categories.

Because direct specification of the full within-cluster correlation matrix becomes increasingly difficult as the number of categories increases, this module implements a pragmatic design-stage strategy in which adjacent ordinal categories are collapsed to form a three-category outcome.

Users can specify one or more clinically meaningful three-category collapsing schemes and compare the resulting power and sample size calculations.

Binary dichotomizations can be specified separately from the three-category collapsing schemes, allowing GEE-Binary calculations to correspond to clinically meaningful event definitions.

The original ordinal outcome does **not** need to be collapsed in the final trial analysis; collapsing is used here as a design-stage approximation for power and sample size calculation.

## Statistical Methods

The application brings together methods developed for two related trial-design problems.

### Ordinal Composite Endpoints

For an ordinal composite endpoint constructed from correlated binary components, the app uses **vine copulas** to derive the joint distribution of the binary components from their marginal prevalences and dependence assumptions.

The joint distribution is then mapped to the ordinal composite endpoint according to the specified composition rule.

The resulting ordinal distribution can subsequently be used for power and sample size calculations.

### Individually Randomized Trials

For IRTs, power and sample size are calculated using **Whitehead's method** for ordinal outcomes under the proportional odds model.

### Cluster Randomized Trials

For CRTs with three-category ordinal outcomes, the application implements closed-form **GEE-based power and sample size formulas** that account for the multidimensional within-cluster correlation structure induced by ordinal outcomes.

The app also provides Whitehead design-effect and binary-collapse approaches for comparison and sensitivity analysis.

For outcomes with four or more categories, the proposed three-category GEE framework can be applied after clinically meaningful adjacent-category collapsing.

## Who Should Use This App?

OrdCompCRTpower is intended for:

- statisticians;
- clinical trialists;
- biostatisticians;
- clinical investigators; and
- applied researchers

who are designing studies with ordinal or ordinal composite outcomes.

The application is particularly useful for exploring how power and required sample size change under different assumptions about:

- ordinal outcome distributions;
- treatment effect sizes;
- number of clusters;
- cluster sizes;
- within-cluster correlations;
- binary dichotomizations; and
- ordinal category collapsing schemes.

The interactive interface is also intended to facilitate communication between statisticians and clinical collaborators by making design assumptions and their consequences easier to visualize.

## Running the Application Locally

Clone this repository:

```bash
git clone https://github.com/Slan1997/OrdCompCRTpower.git
cd OrdCompCRTpower
```

Open the R project in RStudio and install any required packages that are not already installed.

Then launch the application from R using:

```r
shiny::runApp()
```

Alternatively, open the main application file in RStudio and click **Run App**.

## Typical Workflow

A typical workflow depends on the study design.

For an **ordinal composite endpoint**:

```text
Binary component information
        |
        v
Composite module
        |
        v
Derived ordinal distribution
        |
        +---------> IRT module
        |
        +---------> CRT module
```

For an ordinal outcome whose distribution is already known, users can proceed directly to the appropriate IRT or CRT module.

Within the CRT modules, users can vary cluster size, number of clusters, correlation assumptions, and outcome definitions to perform sensitivity analyses and assess alternative trial designs.

## Methodological Background

The statistical methodology implemented in OrdCompCRTpower was developed as part of:

> **Statistical Advancements in Power and Sample Size Determination for Clinical Trials with Ordinal Outcomes**  
> Lan Shi, PhD Dissertation, Vanderbilt University, 2026.

The application integrates two main methodological developments:

1. a vine copula-based framework for deriving distributions of ordinal composite endpoints from correlated binary components; and
2. closed-form GEE-based power and sample size methods for cluster randomized trials with ordinal outcomes.

## Current Limitations

The current version focuses primarily on:

- two-arm trials;
- equal allocation;
- proportional-odds treatment effects; and
- selected within-cluster correlation structures.

The application does not currently support direct upload of pilot datasets for estimating design parameters.

The Composite module currently supports rank-by-severity and summation compositions, together with a built-in COVID-19 ordinal scale example. Arbitrary user-defined mappings from binary component patterns to ordinal categories are not yet implemented.

Potential future extensions include unequal allocation, additional small-sample corrections, direct pilot-data input, more flexible ordinal composite definitions, and more complex CRT designs.

## Source Code

The application is implemented in **R** using **Shiny**.

Repository:

https://github.com/Slan1997/OrdCompCRTpower

## Bug Reports and Questions

If you encounter any bugs or have questions about the application, please feel free to contact me at slannist97@gmail.com or open an issue in this GitHub repository.

## Citation

If you use OrdCompCRTpower in your work, please cite:

> Shi L. *Statistical Advancements in Power and Sample Size Determination for Clinical Trials with Ordinal Outcomes*. PhD Dissertation, Vanderbilt University; 2026.

Additional citations for the methodological papers underlying the Composite and CRT modules will be added as they become available.

## Author

**Lan Shi**  
Department of Biostatistics  
Vanderbilt University School of Medicine

## Disclaimer

OrdCompCRTpower is intended to support study design and sensitivity analysis. Power and sample size calculations depend on the validity of the user-specified assumptions, including outcome distributions, treatment effects, and dependence parameters. Investigators should evaluate multiple plausible assumptions and use clinical and statistical judgment when selecting a final study design.
