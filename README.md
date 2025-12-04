<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->
  
# Signature Recontextualization

The signature recontextualization problem describes a simple goal in computational biology:
Given gene signature X of a genetic or chemical perturbation in model organism Y, what is the corresponding gene signature of the same perturbation in model organism Z?

Any method (gene-regulatory network or deep-learning based) that is able to input a ranked list of genes from one biological context and output another ranked list of genes for another context performs this task of 'signature recontextualization'.

This repository contains benchmarking tasks and data for evaluating these signature recontextualization methods.

## Installation:
```
BiocManager::install("montilab/sigrecon", dependencies=TRUE)
```
