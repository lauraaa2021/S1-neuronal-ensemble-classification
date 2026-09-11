# S1 Neuronal Ensemble Classification

This repository contains the Igor Pro code developed as part of the Master's Thesis:

**Design, Development and Application of a Methodology for Neuronal Ensemble Dynamics Analysis Underlying Sensory Discrimination**

The methodology was developed for the semi-automated classification of calcium-imaging responses recorded from the mouse primary somatosensory cortex (S1) during repeated tactile stimulation.

## Requirements

- Igor Pro 9 (64-bit)
- SARFIA
- NeuroMatic

## Main functions

- `GetCalciumEvents_PreStimWindows1`
- `GetCalciumEvents_PreStimWindows_TwoThresholds`
- `ClassifyByStim2`
- `ClassifyROIEnsembles`
- `OverlayCaPN_WithThresholds`

## Method overview

The pipeline:

1. Calculates local pre-stimulus baselines and standard deviations.
2. Defines positive and negative event-detection thresholds.
3. Converts calcium traces into a signed discrete representation (+1, 0, -1).
4. Classifies individual stimulus responses as ON, OFF, non-responsive, or mixed.
5. Integrates responses across repeated stimuli to obtain the final ON, OFF, or NON classification.
6. Allows visual researcher supervision and correction of ambiguous responses.

## Data availability

The experimental calcium-imaging datasets are not included in this repository.

## Author

Laura Blas Zúmel  
Master of Science in Neurotechnology  
Universidad Politécnica de Madrid  
2026
