""" 
A script to output the values for problem 2c: 
Hubble Time and Lookback Time at z=2 for various values of H_0.

Astro 534 Winter 2015 U. Mich. 
Written by Tom Rice

"""

from __future__ import division

import numpy as np
import astropy.units as u
import astropy.constants as c

def hubble_time(redshift, hubble_constant):
    """ 
    Returns Hubble time at redshift z using hubble constant H_0.

    Assumes H_0 given in km s^-1 Mpc^-1.

    """

    z = redshift
    H_0 = u.Quantity(hubble_constant, u.km / u.s / u.Mpc)

    t_H = 2/(3*H_0) * (1 + z)**(-3/2)

    return t_H.to('Gyr')

def lookback_time(redshift, hubble_constant):

    t_H_0 = hubble_time(0, hubble_constant)
    t_H_z = hubble_time(redshift, hubble_constant)

    return t_H_0 - t_H_z



def calculate_values_for_2c(redshift=2, list_of_H0s=[50,70,100]):

    z = redshift

    print "Hubble Time at z={0}".format(z)
    for H_0 in list_of_H0s:

        print "    H_0={0} km s^-1 Mpc^-1:".format(H_0)
        print "   {0:.2f}".format(hubble_time(z, H_0))

    print "Lookback Time at z={0}".format(z)
    for H_0 in list_of_H0s:

        print "    H_0={0} km s^-1 Mpc^-1:".format(H_0)
        print "   {0:.2f}".format(lookback_time(z, H_0))




