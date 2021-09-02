import numpy as np


def pearson_correlation(X: np.array, Y: np.array, copy_input=True) -> np.array:
    """Vectorized fucntion to calculate pearson correlation coefficient between
        a traits of interest X and a methylation / target matrix Y
        Params:
            * *X (np.ndarray)*: matrix of $$m$$ samples and $$n$$ features
            * *Y (np.ndarray)*: matrix of $$i$$ features and $$j$$ samples
            * *copy_input (bool) = True*: copy input data, if False data may be reshapped
        Returns:
            * *PCC Matrix (np.ndarray)*: $$n$$ by $$i$$ matrix of PCC coefficients
                                         between each $$i$$ features and $$n$$ features
    """
    # reshape X if 1-d
    X_fit = X if not copy_input else np.copy(X, order='k')
    if len(X_fit.shape) == 1:
        X_fit = X_fit.reshape((-1, 1))
    Y_fit = Y if not copy_input else np.copy(Y, order='k')
    if len(Y_fit.shape) == 1:
        Y_fit = Y_fit.reshape((1, -1))
    # calculate mean for each row and phenotype mean
    Y_means = np.mean(Y_fit, axis=1)
    X_means = np.mean(X_fit, axis=0)

    # subtract means from observed values
    Y_trans = Y_fit - Y_means.reshape((-1, 1))
    X_trans = X_fit - X_means

    # calculate covariance
    covariance = np.array([np.sum(Y_trans * X_trans[:, count], axis=1) for count in range(X_fit.shape[1])])
    Y_var = np.sqrt(np.sum(Y_trans ** 2, axis=1))
    Y_var = np.array([Y_var for _ in range(X_fit.shape[1])])
    X_var = np.sqrt(np.sum(X_trans ** 2, axis=0))
    return covariance / (Y_var * X_var.reshape((-1, 1)))
