#!/usr/bin/env python

import numpy as np


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

        # --- Make sure things are copacetic prior to setting things...

        # --- Here, we check to make sure everything is copacetic.
        if (len(self.forecast.shape) > 2) and (self.forecast.shape[0] > 1):
            # --- First, make sure we're dealing with same number of vars between fcst/train
            if self.forecast.shape[0] != self.train.shape[0]:
                raise ValueError("Different number of variables between forecast and training data arrays!")
            # --- Make sure there aren't more weights than variables
            if len(weights) > self.forecast.shape[0]:
                raise ValueError("More weights than variables!")
            if len(comp_method) > self.forecast.shape[0]:
                raise ValueError("More methods than variables!")

        # --- Everything cool? Ok, let's do this.
        self.window = window
        self.comp_method = comp_method
        self.weights = weights

        return self

    def domain(self, lat_bounds=[], lon_bounds=[], forecast_lats=[], forecast_lons=[], lat_inc=1, lon_inc=1):
        """
        Domain setup.

        Parameters
        ----------
        lat_bounds/lon_bounds : list or NumPy array
            A list/numpy array of all latitude within the forecast and training data.\n
            len(allLats) = forecastData[0] or forecastData[1] if forecastData.shape == 3.
        forecast_last/forecast_lons : list or NumPy array
            Either a list/numpy array of min/max forecast lats/lons (defining a domain) or a single point.


        comp_method : list
            Method with which to pattern match. Options: ['Rank','MAE','RMSE','Corr'].
            if n_vars > 1, can use different methods to compare each variable.
        """

        # --- Make sure things are copacetic prior to setting things...

        # --- Here, we check to make sure everything is copacetic.
        if (len(self.forecast.shape) > 2) and (self.forecast.shape[0] > 1):
            # --- First, make sure we're dealing with same number of vars between fcst/train
            if self.forecast.shape[0] != self.train.shape[0]:
                raise ValueError("Different number of variables between forecast and training data arrays!")
            # --- Make sure there aren't more weights than variables
            if len(weights) > self.forecast.shape[0]:
                raise ValueError("More weights than variables!")
            if len(comp_method) > self.forecast.shape[0]:
                raise ValueError("More methods than variables!")

        # --- Everything cool? Ok, let's do this.
        self.window = window
        self.comp_method = comp_method
        self.weights = weights

        return self


    def compare(fcst, train, latidxs, lonidxs, weight, method):
        """


        :param weight:
        :param method:
        :param fcst:
        :param train:
        :param latidxs:
        :param lonidxs:
        """
        return fcst

    def analog_single_date(forecastData, trainingData, forecastLats, forecastLons, allLats, allLons, \
                           window, **kwargs):
        """
        Method to pattern match and find analogous dates at either single location or multiple grid points.

        Used for a single forecast date.


        forecastData - Either a 2-d (lat,lon) or 3-d (n_vars,lat,lon) numpy array of "current" forecast data.
        trainingData - Either a 3-d (time,lat,lon) or 4-d (n_vars,time,lat,lon) numpy array of "past"/training data.
        forecastLats - Either a list/numpy array of min/max forecast lats (defining a domain) or a single point.
        forecastLons - Either a list/numpy array of min/max forecast lons (defining a domain) or a single point.
        allLats - A list/numpy array of all latitudes within the forecast and training data.\n
            len(allLats) = forecastData[0] or forecastData[1] if forecastData.shape == 3.
        allLons - A list/numpy array of all longitudes within the forecast and training data.\n
            len(allLons) = forecastData[0] or forecastData[1] if forecastData.shape == 3.
        window - +/- number of grid points (n/s, e/w) around each forecast grid point to include in local domain pattern\n
            matching. For example, if window = 3, we use 49 grid points to calculate differences between forecast/training\n
            data (3*2 + 1)*(3+2 +1). If one_dom = True, this is ignored.


        optional:
            method - Method with which to pattern match. Options: ['Rank','MAE','RMSE','Corr'].
                if n_vars > 1, can use different methods to compare each variable.
            one_dom - If we're forecasting for multiple grid points, this can set it to pattern match one big domain.
            weights - If interested in an additive model, weights can be set based on the importance of each item.
                note: sum(weights) == 1.


        returns:
            analog_idxs - indices of closest analogs, from best pattern match to worst.
        """
        meths = kwargs.get('method', ['Rank'])
        one_dom = kwargs.get('one_dom', False)
        weights = kwargs.get('weights', [])
        n_proc = kwargs.get('n_proc', 1)

        # --- Here, we check to make sure everything is copacetic.
        if (forecastData.shape > 2) and (forecastData.shape[0] > 1):
            # --- First, make sure we're dealing with same number of vars between fcst/train
            if forecastData.shape[0] != trainingData.shape[0]:
                raise ValueError("Different number of variables between forecast and training data arrays!")
            # --- Make sure there aren't more weights than variables
            if len(weights) > forecastData.shape[0]:
                raise ValueError("More weights than variables, exiting...")
            if len(meths) > forecastData.shape[0]:
                raise ValueError("More methods than variables, exiting...")


        # --- Analog method if we're pattern matching one large domain
        #if one_dom:



