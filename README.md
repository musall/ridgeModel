# ridgeModel

Linear regression model, using ridge regression to fit a combination of task- or movement-events to high-dimensional neural data sets.

The script 'tutorial_linearModel.m' shows how to build a large design matrix and fit it to neural data. To use it, you need the [widefield data set](http://repository.cshl.edu/38599/) from our paper ‘Single-trial neural dynamics are dominated by richly varied movements’ by Musall, Kaufman et al., 2019.
The tutorial will use an example widefield recording and extract different events (e.g. stimulus onset, whisking, etc) that are used to model the neural data. Have a look at the [paper](https://www.nature.com/articles/s41593-019-0502-4) for a more detailed explanation of the model.

Each event is transformed into a 'design matrix' which is needed to estimate neural responses to each moment in time, before and after the event occurred. This is done by the function 'makeDesignMatrix' which turns an event matrix where each column is a trace of binary events into a larger design matrix and also returns an index that indicates which regressors in the new matrix correspond to the initial event trace.
In the tutorial, we are using three different event types that influence the range of regressors in the created design matrix.
'Whole-trial' matrices (eventType = 1) include every sample after the occurrence of the event until the end of the trial. We use this to estimate trials of a certain type, e.g. trials that were rewarded.
'Post-event' matrices (eventType = 2) include every sample after event occurrence until time 'sPostTime' has been reached. We use this for example to fit a stimulus-response.
'Peri-event' matrices (eventType = 3) include samples before the event, defined by 'mPreTime', and samples after the event, defined by 'mPostTime'. We usually use this for movements.
The model can also include continuous regressors. Just add them in as additional columns to the design matrix. In general, it is good practice to keep all regressors zero-mean and normalize their magnitude by dividing them by their standard deviation. The fitting function 'ridgeMML' will also do this internally by default.

After the design matrix ('fullR') is created, it is fitted against a neural data set ('Vc'), using the 'ridgeMML' function. Here, neural data is a matrix with each column being a recording of some kind (e.g. a single neuron or in the tutorial a widefield dimension).
RidgeMML estimates the ridge penalty lambda for each column in Vc seperately ('ridgeVals') and also returns the corresponding beta weights for each regressor in fullR ('dimBeta').

The line 

```Vm = (fullR * dimBeta)'```

creates a reconstructed data set 'Vm' based on the model fits.

The beta weights can then be explored to assess the neural responses to different events in the model.
Lastly, the tutorial shows how to run cross-validation using the 'crossValModel' function to test the model's predictive performance on left-out data.


# Other code, related to the paper

The repo also includes other functions that were either used in the paper or somewhat simplified to provide an intuition of the neural data in general. They are located in the 'delayedDetection' folder.

'delayDec_rebuildAvg' creates average responses of the widefield data for different times in the behavioral paradigm. It also creates d' maps for choice and sensory modality.

'delayDec_regressModel' is the model that was used in the paper. All its functionality is also included in the tutorial.

'delayDec_compareRegressors' compares cvR² and dR² for different model variables, creates plots that show increased spatial specificity for dR² over cvR² and computes dR² for each time point in the trial.

'delayDec_compare2pRegressors' is mostly similar to 'delayDec_compareRegressors' but is used for 2-photon data.

'delayDec_Neuropixels' fits a simplified version of the linear model to Neuropixels recordings in awake mice during visual stimulation. Subsequently, use 'delayDec_NeuropixelsResults' to visualize the results.