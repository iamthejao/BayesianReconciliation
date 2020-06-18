# Reconciling hierarchical forecasts via Bayes' rule

This package implements  the Bayesian reconciliation of hierarchical forecasts.

In particular the function `hierRecBayesianExperiment` reconciles forecasts for hierachical / grouped time series.
The reconciliation algorithm is implemented in `bayesRecon` in reconciliationMethods.R


## Required packages for reproducing the paper experiments

Required packages should be automatically installed by the commands in the beggining of the scripts.

* `forecast` (produces base forecasts with either auto.arima or ets)
* `hts`   (algorithms for reconciling hierarchical time series, including minT)
* `huge`, `SHIP`, `corpcor` (covariance matrix estimation)
* `matlib`, `matrixStats`, `wordspace` (collection of matrix functions and statistics)
* `MASS`, `tmvtnorm`, `truncnorm` (truncated normal and other density functions)
* `Rcpp` (R and C++ interface)
* `optparse` (argument option parser)


## Choice of data sets and base forecasts
Two real data sets can be used: `infantgts` (available from `hts`) or `tourism` (raw data available from [https://robjhyndman.com/publications/mint/](https://robjhyndman.com/publications/mint/). The csv file of tourism is available in this repository. When the   `tourism` data set is selected, function `loadTourism.R` reads the csv file and builds the hierarchy.

Alternatively, the data can be synthetically generated:
* from a hierarchy having two bottom time series and one top time series  (`synthetic`).
* from a hierarchy having four bottom time series, two intermediate and one top time series  (`syntheticLarge`).

The base forecasts can be created using either `auto.arima` or `ets`, both available from `forecast`.

## Examples

```R
 hierRecBayesianExperiment(dset="infantgts", fmethod="ets", h=1, iTest=1) 
 hierRecBayesianExperiment(dset="infantgts", fmethod="arima", h=1, iTest=2) 
 hierRecBayesianExperiment(dset="tourism", fmethod="ets", h=2, iTest=1)
 hierRecBayesianExperiment(dset="tourism", fmethod="arima", h=1) 
```


Examples with generated data sets:
```R
 hierRecBayesianExperiment(dset="syntheticLarge", h=3, synth_n=300)  # hierarchy is 4-2-1 
```
Arguments:

* `dset` : can be either `infantgts`, `tourism`, `synthetic` or `syntheticLarge` 

* `fmethod` : method for generating the base forecasts: it can be either `arima` or `ets`. Default: 'ets'.

* `h`: forecast horizon for which we reconcile the forecasts. We use between 1 and 4 in the experiments of the paper.
Default: 1.

* `iTest`: controls how the split the data between train and test. The training set contains the data from the first observation up to the observation in position (length(timeSeries) - h - (iTest - 1)). This is useful for parallelizing the experiments. Admissible values are between 1 and 45. Default: 1.

An additional set of parameters applies to set of experiments with synthetic time series.

* `seed` : seed (default:0)

* `synth_n` : length of the generated time series. Default: 100.

* `synthCorrel` : correlation of the noise affecting the bottom time series. Applies only to the `synthetic` case; default: 0.5.
For the `syntheticLarge` case, the covariance matrix of the noise is set as in the MinT paper (Wickramasuriya et al., 2019)

* `savePredictions` : save predictions on prediction folder. default: TRUE

* `saveSamples` : save samples on sample folder. It can be heavy according to sampleSize x hierarchy size. default: FALSE

* `enforceKhOne` : kh=1 (TRUE) or kh=h (FALSE) experiment. default: FALSE

* `sampleSize` : Sample size for sampled mean, median, ES, etc. default: 100000

The reconciliation results are saved scattered in the folder `results` according to the batch split used with the name following the pattern `{h}_{dset}_{fmethod}_{from_split}to{to_split}.csv` as in `2_tourism_ets_28to35`
The file shows the mae, rmse and energy score for many combinations of groups of the hierarchy (All, Upper, Bottom, etc) and reconciliation schemes (BottomUp, pMinT, LG) and the Base forecasts.

## Running batch experiment

The experiments can be ran all at once (or more easily parallelized) by running the batch experiment script which implements the calls to the function `hierRecBayesianExperiment`.

```R
 Rscript batchHier.R -d infantgts -m arima -p "results/" -k true
```

For help with the arguments:

```R
Rscript batchHier.R --help
```

## pMinT and LG reconciliation function

The algorithm is implemented in `bayesRecon` inside `reconciliationMethods.R`. It can be used as following:

```R
bayesRecon(basePredictions, sumMatrix, residualCovarianceMatrix, reconType="pmint")
bayesRecon(basePredictions, sumMatrix, residualCovarianceMatrix, reconType="lg")
```

* `preds` : are the h-step predictions of each node

* `mSumMatrix` : is the hierarchy sum matrix

* `mCovar` : is the hierarchy residuals (estimated) covariance matrix

* `sampleSize` : sample size. default: 100000

* `method` : can be either `pmint` for pMinT or `lg` for Linear Gaussian (LG).

* `kh`: should be either 1 or h.

## Analyzing the results
The previous functions save the raw results in a csv file within the `results/` directory. The results are then analysed using a python notebook provided along with the R code.s


