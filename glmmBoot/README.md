Mixed Poisson Implementation
================

This code to illustrates the implementation of the bootstrapping methods
described in [(Flores-Agreda and Cantoni, Under
Review)](https://www.researchgate.net/publication/315768128_Bootstrapping_Generalized_Linear_Mixed_Models_via_a_Weighted_Laplace_Approximation)
in a real-data context. Data comes from the  example found in
[Molenberghs & Verbeke
(2006)](https://www.springer.com/gp/book/9780387251448) initially
reported by [Faught et. al.
(1996)](https://www.ncbi.nlm.nih.gov/pubmed/8649570).

The aim of the study was to verify the effects of an new anti-epileptic
drug (AED) compared to a placebo on the number of seizures experienced
by patients during the study.
<!-- % After that period, 45 patients were assigned to the placebo group, 44 to the active (new) treatment group.  -->
<!-- 89 Patients were assigned to either group (placebo or treatment), measured on a weekly basis and followed during 16 weeks, after which they were entered into a longer-term study with the number of visits ranging between 2 and 27 weeks, hence making clusters variable in size. -->
<!-- % The out- come of interest is the number of epileptic seizures experienced during the last week, i.e., since the last time the outcome was measured. The key re- search question is whether or not the additional new treatment reduces the number of epileptic seizures. -->
To do this, consider a mixed Poisson model for the outcome containing
two potentially correlated random effects: one for a random intercept
and another one for the visit time i.e. \[alt
text\]\[<https://github.com/dfloresag/glmmBoot/blob/master/glmmBoot/img/eq01.png>\]
- \(T_{ij}\) represents the effect of the treatment and - \(t_{ij}\) the
visit time.

The variance-covariance structure of the vector of Normal random effects
\(\mathbf{u}_i = [u_{i1}, u_{i2}]^T\), comprises a correlation
coefficient \(\rho\),
i.e.

\[\mathbf{\Delta}(\mathbf{\sigma}) = \left[\begin{array}{c c} \sigma_{1}^2 & \rho\sigma_{1}\sigma_{2} \\  \rho\sigma_{1}\sigma_{2} & \sigma_{2}^2 \end{array}\right]\].

Analysis has been performed using TMB to provide Laplace-Approximated
estimates and their Standard Errors obtained via numerical
approximations of the required derivatives of the LALL, while RWLB and
PB replicates of the model parameters were generated using the
previously described algorithms. The more accurate AGQ estimate is
unavailable in software for models with a random slope, therefore, only
our LAML is available.

The Bootstrap distributions were then used to construct (i) Estimates of
the Parameters, by averaging over the replicates (ii) Estimates of the
Standard Errors by computing the standard deviations for each
distribution and (iii) percentile-based Confidence Intervals (CI) with a
level of 95%.