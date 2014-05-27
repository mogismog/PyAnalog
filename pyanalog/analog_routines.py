#!/usr/bin/env python

import warnings

import numpy as np

from .utils import find_nearest_idx

import fortran_routines as fort

# --- First, some exception classes
class AnalogError(Exception):
    pass
class DomainError(Exception):
    pass

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

        # --- Let's also set some defaults here used by other functions
        self.window = 3
        self.comp_method = []
        self.weights = []

    def analog_params(self, window=3, weights=[], comp_methods=['Rank']):
        """
        Set up some parameters for whatever analog method.

        Parameters
        ----------
        window : integer
            +/- number of grid points (n/s, e/w) around each forecast grid point to include in local domain pattern\n
            matching. For example, if window = 3, we use 49 grid points to calculate differences between forecast/training\n
            data (3*2 + 1)*(3+2 +1). Default = 3
        comp_method : list
            Method with which to pattern match. Options: ['Rank','MAE','RMSE','Corr'].
            if n_vars > 1, can use different methods to compare each variable.
            Default = ['Rank']
        weights : list
            If we are using more than one variable, then one field may be more important to the predictor than others.
            As such, we can weigh each variable by multiplying the pattern matching results by some value n, where 0<n<1.
            IMPORTANT: Be sure that sum(n) == 1.
            Default = None
        """

        # --- Here, we check to make sure everything is copacetic.
        if (len(self.forecast.shape) > 2) and (self.forecast.shape[0] > 1):
            self.n_vars = self.forecast.shape[0]
            # --- First, make sure we're dealing with same number of vars between fcst/train
            if self.forecast.shape[0] != self.train.shape[0]:
                raise AnalogError("Different number of variables between forecast and training data arrays!")
            # --- Make sure there aren't more weights than variables
            if len(weights) > self.forecast.shape[0]:
                raise AnalogError("More weights than variables!")
            if len(comp_methods) > self.forecast.shape[0]:
                    raise AnalogError("More methods than variables!")
        else:
            self.n_vars = 1

        # --- Now, to make sure the method selected is available, and if it is then add
        available_methods = ['rank','rmse','mae','corr']
        self.comp_method = []
        for mthd in comp_methods:
            if any(mthd.lower() in avail for avail in available_methods):
                self.comp_method.append(mthd.lower())
            else:
                warnings.warn('Method {} not an option! Reverting to rank method'.format(mthd.lower()))
                self.comp_method.append('rank')

        # --- Everything cool? Ok, let's do this.
        self.window = window
        self.weights = weights

        return self

    def domain(self, lat_bounds=[], lon_bounds=[], forecast_lats=[], forecast_lons=[], lat_inc=1, lon_inc=1):
        """
        Sets up the forecast domain.

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
            raise DomainError("len(lat_bounds) != len(lon_bounds)\nWe either need to forecast for a point or domain.")

        # --- Everything cool? Ok, let's do this.
        if len(forecast_lats) == 1:
            self.point_fcst = True
        else:
            self.point_fcst = False

        self.allLats = np.arange(lat_bounds[0], lat_bounds[-1] + 1, lat_inc)
        self.allLons = np.arange(lon_bounds[0], lon_bounds[-1] + 1, lon_inc)
        # --- If we're dealing with a point forecast, we will find the closest lat/lon grid points to the fcst point.
        if len(forecast_lats) > 1:
            self.closest_lat = self.allLats[find_nearest_idx(self.allLats[:], self.fcstLats[0])[0]]
            self.closest_lon = self.allLons[find_nearest_idx(self.allLons[:], self.fcstLons[0])[0]]
        else:
            self.fcstLats = forecast_lats
            self.fcstLons = forecast_lons

        return self

    def closest_dates_1fcst(self):
        """
        Used to find analogs for a single forecast domain.

        returns:
            analog_idxs - indices of closest analogs, from best pattern match to worst.
        """

        # --- Pre-generating analog indices array, this should be faster.
        if len(self.fcstLats) <= 1:
                analog_idxs = np.zeros(self.train.shape[0])
        else:
            if self.n_vars == 1:
                analog_idxs = np.zeros(self.train.shape)
            else:
                analog_idxs = np.zeros((self.train.shape[0], self.train.shape[2], self.train.shape[3])) * -9999.9

        # --- Ok, let's do this...
        if self.point_fcst:
            if self.n_vars > 1:
                for nv in self.n_vars:
                    if self.comp_method[nv] == 'rank':
                        analog_idxs[...] = fort.rank_analog_point(self.train, self.forecast, self.train.shape[0],\
                                                                  self.allLats.shape[0],self.allLons.shape[0], \
                                                                  self.closest_lat, self.closest_lon, self.window)
                    if self.comp_method[nv] == 'rmse':
                        analog_idxs[...] = fort.rmse_analog_point(self.train, self.forecast, self.train.shape[0],\
                                                                  self.allLats.shape[0],self.allLons.shape[0], \
                                                                  self.closest_lat, self.closest_lon, self.window)
                    if self.comp_method[nv] == 'mae':
                        analog_idxs[...] = fort.mae_analog_point(self.train, self.forecast, self.train.shape[0],\
                                                                  self.allLats.shape[0],self.allLons.shape[0], \
                                                                  self.closest_lat, self.closest_lon, self.window)
                    if self.comp_method[nv] == 'corr':
                        analog_idxs[...] = fort.corr_analog_point(self.train, self.forecast, self.train.shape[0],\
                                                                  self.allLats.shape[0],self.allLons.shape[0], \
                                                                  self.closest_lat, self.closest_lon, self.window)


