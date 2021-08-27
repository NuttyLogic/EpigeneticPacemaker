import numpy as np


class EPMScaler:
    """Simple min, max scaler implementation for array with m samples and n features"""

    def __init__(self):
        """
        attributes:
            self.X_min: np.ndarray = min of m rows
            self.X_max: np.ndarray = max of m rows
        """
        self.X_min = None
        self.X_max = None

    def fit(self, X: np.ndarray):
        """Set min and max values"""
        self.X_min = np.min(X, axis=0)
        self.X_max = np.max(X, axis=0)

    def transform(self, X: np.ndarray) -> np.ndarray:
        """Scale input array to reference distribution range"""
        # Scale input value to 0 - 1
        X_std = (X - np.min(X, axis=0)) / (np.max(X, axis=0) - np.min(X, axis=0))
        # Scale values to original reference range
        return X_std * (self.X_max - self.X_min) + self.X_min