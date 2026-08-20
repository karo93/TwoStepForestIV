# TwoStepForestIV

Replication materials for Two-Step Forest-Based Subgroup Discovery with Instrumental Variables, which studies heterogeneous treatment effects and subgroup discovery under imperfect compliance. Building on the two-step framework of Bayesian Causal Forest with Instrumental Variables (BCF-IV), the repository implements a model-agnostic first-stage discovery procedure in which heterogeneous complier treatment effects can be estimated using alternative machine learning methods.

The repository provides implementations of two non-Bayesian forest-based approaches: DRRF-IV, a debiased transformed-outcome regression-forest estimator, and GRF-IV, an instrumental-variable adaptation based on the generalized random forest framework. Their performance is compared with BCF-IV in simulation studies and an empirical application.

# Overview

Instrumental variable methods provide a way to estimate causal effects when treatment assignment is endogenous or different from the actual received treatment. In settings with imperfect compliance, the treatment effect of interest is the conditional complier average causal effect (CCACE) rather than the standard CATE. This repository studies a two-step procedure for identifying interpretable subgroups with heterogeneous complier treatment effects. In the first step, a machine learning estimator is used to estimate individual-level heterogeneous treatment effects. These estimated effects are subsequently used as outcomes in a CART procedure to identify a tree partition capturing treatment-effect heterogeneity. In the second step, the discovered tree structure is held fixed and subgroup-specific causal effects are estimated using instrumental-variable regressions within the terminal nodes.

The main contribution is to make the discovery step model-agnostic, allowing different machine learning estimators to be incorporated into the two-step framework. The repository focuses on two forest-based implementations:

* DRRF-IV, which uses estimated nuisance functions to construct a debiased transformed outcome and estimates heterogeneous complier effects with a regression forest.
* GRF-IV, which uses the instrumental-variable functionality of generalized random forests to estimate heterogeneous complier effects.

The empirical application revisits an instrumental-variable study of prompt intensive care unit admission and 28-day mortality across UK National Health Service hospitals, with bed availability serving as the instrument.
