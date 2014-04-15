#!/usr/bin/env python

import numpy as np
from scipy.stats import rankdata as rd
from ..utils import find_nearest_idx


class Analog(object):
    """Various methods to produce deterministic/probabilistic forecasts via analogs."""

    def __init__(self, forecast, train, ):
        """
        Initialize the analog object.

        Parameters
        ----------
        fcst : NumPy array
            Either a 2-d (lat,lon) or 3-d (n_vars,lat,lon) numpy array of "current" forecast data.
        train : NumPy array
            Either a 3-d (time,lat,lon) or 4-d (n_vars,time,lat,lon) numpy array of "past"/training data.
        """


        # --- Set up the class attributes
        self.forecast = forecast
        self.train = train

        return self

    def analog_params(self, window=3, comp_method=['Rank'], weights=[]):
        """
        Set up some parameters for whatever analog method.

        Parameters
        ----------
        window : integer
            +/- number of grid points (n/s, e/w) around each forecast grid point to include in local domain pattern\n
            matching. For example, if window = 3, we use 49 grid points to calculate differences between forecast/training\n
            data (3*2 + 1)*(3+2 +1).
        comp_method : list
            Method with which to pattern match. Options: ['Rank','MAE','RMSE','Corr'].
            if n_vars > 1, can use different methods to compare each variable.
        """

        # --- Here, we check to make sure everything is copacetic.
        if (len(self.forecast.shape) > 2) and (self.forecast.shape[0] > 1):
            self.n_vars = self.forecast.shape[0]
            # --- First, make sure we're dealing with same number of vars between fcst/train
            if self.forecast.shape[0] != self.train.shape[0]:
                raise ValueError("Different number of variables between forecast and training data arrays!")
            # --- Make sure there aren't more weights than variables
            if len(weights) > self.forecast.shape[0]:
                raise ValueError("More weights than variables!")
            if len(comp_method) > self.forecast.shape[0]:
                raise ValueError("More methods than variables!")
        else:
            self.n_vars = 1

        # --- Everything cool? Ok, let's do this.
        self.window = window
        self.comp_method = comp_method
        self.weights = weights

        return self

    def domain(self, lat_bounds=[], lon_bounds=[], forecast_lats=[], forecast_lons=[], lat_inc=1, lon_inc=1):
        """
        Main and forecast domain setup.

        Parameters
        ----------
        lat_bounds/lon_bounds : list or NumPy array
            A list/numpy array of all latitude within the forecast and training data.\n
            len(allLats) = forecastData[0] or forecastData[1] if forecastData.shape == 3.
        forecast_last/forecast_lons : list or NumPy array
            Either a list/numpy array of min/max forecast lats/lons (defining a domain) or a single point.
        lat_inc/lon_inc : float
            The change in latitude/longitude each grid point (e.g. dx or dy).
        """
        # --- Here, we check to make sure everything is copacetic.
        if len(forecast_lats) != len(forecast_lons):
            raise ValueError("len(lat_bounds) != len(lon_bounds)\nWe either need to forecast for a point or domain.")

        # --- Everything cool? Ok, let's do this.
        if len(forecast_lats) == 1:
            self.point_fcst = True
        else:
            self.point_fcst = False

        self.allLats = np.arange(lat_bounds[0], lat_bounds[-1] + 1, lat_inc)
        self.allLons = np.arange(lon_bounds[0], lon_bounds[-1] + 1, lon_inc)
        self.fcstLats = forecast_lats
        self.fcstLons = forecast_lons

        return self

    def compare(self, fcst, train, idxs):
        """
        Does the heavy lifting in finding analogs.

        Parameters
        ----------
        self : object
            Analog object.
        fcst : NumPy array
            Local domain forecast data, comparing this to past forecasts.
        train : NumPy array
            Local domain training data, "past" forecasts, reanalysis, whatever.
        idxs : NumPy array
            Empty array, used for indices for closest matching patters.

        Returns
        -------
        idxs : NumPy array
            Array of indices for closest matching patters.
        """

        if self.n_vars == 1:
            if self.comp_method == "Rank":
                temp_ranks = np.zeros((self.train[0]+1,self.train[0,...].reshape(-1).shape[0]))

    def analog_single_date(self):
        """
        Used to find analogs for a single forecast domain.

        returns:
            analog_idxs - indices of closest analogs, from best pattern match to worst.
        """

        # --- Pre-generating analog indices array, this should be faster.
        if len(self.fcstLats) <= 1:
            analog_idxs = np.ones(self.train.shape[0]) * -9999.9
        else:
            analog_idxs = np.ones(self.train.shape) * -9999.9 if self.n_vars == 1 else \
                analog_idxs = np.ones((self.train.shape[0], self.train.shape[2], self.train.shape[3])) * -9999.9

        # --- Ok, let's do this...
        if self.point_fcst:
            # --- We only need one local analog for this, so find closest lat/lon grid point and form from there.
            nearest_lat_idx = find_nearest_idx(self.allLats[:], self.fcstLats[0])[0]
            nearest_lon_idx = find_nearest_idx(self.allLons[:], self.fcstLons[0])[0]

            # --- We then form a domain around the point in question...
            temp_forecast = self.forecast[int(nearest_lat_idx - self.window):int(nearest_lat_idx - self.window) + 1, \
                            int(nearest_lon_idx - self.window):int(nearest_lon_idx - self.window) + 1]

            temp_train = self.train[..., int(nearest_lat_idx - self.window):int(nearest_lat_idx - self.window) + 1, \
                         int(nearest_lon_idx - self.window):int(nearest_lon_idx - self.window) + 1]

            # --- We then find the closest analogs
            idxs = compare(self,temp_forecast,temp_train)




