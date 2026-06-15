# ipeval: Interventional Prediction Evaluation

Provides methods to evaluate predictive performance of models that
estimate risks under hypothetical intervention scenarios
(interventional/causal/counterfactual predictions) with observational
data subject to treatment-outcome confounding. Inverse probability of
treatment weighting (IPTW) is used to construct a pseudopopulation in
which all individuals receive a specified intervention, enabling
assessment of agreement between predicted risks under the intervention
and observed outcomes in the pseudo-population corresponding to that
intervention. Supports interventions with binary or categorical
treatment levels, applied either at a single time point or as
longitudinal treatment strategies with sequential treatment decisions.
Performance measures supported are AUC (Area Under the receiving
operating characteristic Curve), Brier score, observed-expected ratio,
and calibration plots. Methods implemented in this package are based on
work by Keogh and Van Geloven (2024)
[doi:10.1097/EDE.0000000000001713](https://doi.org/10.1097/EDE.0000000000001713)
.

## See also

Useful links:

- <https://github.com/survival-lumc/ipeval>

- <https://survival-lumc.github.io/ipeval/>

- Report bugs at <https://github.com/survival-lumc/ipeval/issues>

## Author

**Maintainer**: Jasper van Egeraat <j.w.a.van_egeraat@lumc.nl>

Authors:

- Nan van Geloven \[copyright holder\]

- Ruth Keogh \[copyright holder\]

Other contributors:

- Leiden University Medical Center \[funder\]
