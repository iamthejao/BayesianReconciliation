# BayesianReconciliation

Package dependence should be handled automatically inside.

Batch run:

The experiments can be run in batch by using the script batchHier.R

Rscript batchHier.R -d infantgts -m arima -p "results/" -k true

For help with the arguments:
Rscript batchHier.R --help


Single run:

Most important function to run the experiments is hierRecBayesianExperiment.R

Arguments:

dset: (infantgts, tourism, synthetic, syntheticLarge)
h: steps ahead prediction
fmethod: (arima, ets)
iTest: used for train/test split
Seed
synth_n: size in timesteps of synthetic large dataset
synthCorrel: correlation for synthetic dataset
testProbability: frequency of MinT check
savePredictions
saveSamples
runPositive
enforceKhOne: kh=1 or kh=h
sampleSize
