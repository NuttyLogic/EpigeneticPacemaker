# EpigeneticPacemaker

<h3>Parameters</h3>
---

**iter_limit: int, default=100**
> The maximum number of iterations

**error_tolerance: float, default=0.00001**
> The tolerance for optimization between iterations, if the error between two consecutive iterations is below the 
tolerance threshold the model will stop fitting

---
 
<h3>Attributes</h3>

**EPM: Dict[model_paramaters]**

- MC_error - initial model error after first model iteration
- MC_rates - initial site rates after first model iteration
- MC_intercepts - initial sites intercepts after first model iteration
- MC_rss - initial residual sum of squares
- EPM_error - observed fit model error 
- EPM_rates - fit model rates
- EPM_intercepts - fit model intercepts
- EPM_iter - number of iteration used to fit model 
- EPM_rss - fit model residual of sum squares
- chi2 - log(MC_rss / EPM_rss) * number_sites * number_samples
- pval - chi2 p 
---
<h3>Methods</h3>

**fit(meth_array: numpy.array, states: numpy.array)**
> Fit an EPM model with an m x n methylation array and n methylation states, where m is the number of methylation sites 
and n is the number of samples. The meth_array and states should be passed an numpy arrays or an assertion error will be thrown. 

**predict(meth_array: numpy.array) -> predicted epigenetic states: numpy.array**
> Predict n epigenetic states given an m x n methylation array

**score(meth_array: numpy.array, states: numpy.array) -> Pearson R: tuple(significance, R value)**
> Predict n epigenetic states given an m x n methylation array and compare against known values. *Note*, predicted values may 
not be linearly dependent on input values.  

***

# EpigeneticPacemakerCV

<h3>Parameters</h3>
---

**iter_limit: int, default=100**
> The maximum number of iterations

**error_tolerance: float, default=0.00001**
> The tolerance for optimization between iterations, if the error between two consecutive iterations is below the 
tolerance threshold the model will stop fitting

**verbose:  bool, default=False**
> Verbose model fitting

**cv_folds: int, default=3**
> Number of cross validation folds 

**randomize_order: bool, default=True**
> Randomize the test folds, otherwise step through the CV folds in order

---
 
<h3>Attributes</h3>

**EPM: Dict[model_paramaters]**

- EPM_rates - CV model rates
- EPM_intercepts - CV model intercepts

**models: Dict[cv models]**
> Cross validation models over every fold, holds all model attributes

**predicted_states: numpy.array**
> Predicted epigenetic state for test samples, ie. samples left out of model fitting for a particular fold

---
<h3>Methods</h3>

**fit(meth_array: numpy.array, states: numpy.array)**
> Fit an EPM model with an \(m \times n\) methylation array and n methylation states, where m is the number of methylation sites 
and n is the number of samples. The meth_array and states should be passed an numpy arrays or an assertion error will be thrown. 

**predict(meth_array: numpy.array) -> predicted epigentic states: numpy.array**
> Predict n epigenetic states given an m x n methylation array

**score(meth_array: numpy.array, states: numpy.array) -> Pearson R: tuple(significance, R value)**
> Predict n epigenetic states given an m x n methylation array and compare against known values. *Note*, predicted values may 
not be linearly dependent on input values.  
