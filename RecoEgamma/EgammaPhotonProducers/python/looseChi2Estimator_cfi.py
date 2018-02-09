import FWCore.ParameterSet.Config as cms

import TrackingTools.KalmanUpdators.Chi2MeasurementEstimatorESProducer_cfi
chi2CutForConversionTrajectoryBuilder = TrackingTools.KalmanUpdators.Chi2MeasurementEstimatorESProducer_cfi.Chi2MeasurementEstimator.clone()
chi2CutForConversionTrajectoryBuilder.ComponentName = 'eleLooseChi2'
chi2CutForConversionTrajectoryBuilder.MaxChi2 = 100000.
chi2CutForConversionTrajectoryBuilder.nSigma = 3.

