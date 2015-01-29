""" 
A script to create plots for problem 4: 
Hubble Time and angular diameter distance for various cosmologies,
from z=0 to z=10.

Astro 534 Winter 2015 U. Mich. 
Written by Tom Rice

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as c

from scipy.integrate import quad


def _hubble_time_inner_function(z, omega_matter, omega_lambda):
    return 1 / ((1 + z) * (omega_matter * (1 + z)**3 + omega_lambda)**(1/2))


def hubble_time_general(redshift, hubble_constant, omega_matter, omega_lambda):
    """
    Returns Hubble time at redshift z using hubble constant H_0, omega_matter, omega_lambda.

    Assumes H_0 given in km s^-1 Mpc^-1.

    """

    z = redshift
    H_0 = u.Quantity(hubble_constant, u.km / u.s / u.Mpc)

    hubble_time_integrand = lambda z: _hubble_time_inner_function(
        z, omega_matter, omega_lambda)

    t_H_z = 1 / H_0 * quad(hubble_time_integrand, z, np.inf)[0]

    return t_H_z.to(u.Gyr)


def _angular_diameter_inner_function(z, omega_matter, omega_lambda):
	return 1 / (omega_matter * (1+z)**3 + omega_lambda)**(1/2)


def angular_diameter_distance(redshift, hubble_constant, omega_matter, omega_lambda):
    """
    Returns angular diameter distance at redshift z using H_0, omega_matter, omega_lambda.

    Assumes H_0 given in km s^-1 Mpc^-1.

    """

    z = redshift
    H_0 = u.Quantity(hubble_constant, u.km / u.s / u.Mpc)

    angular_diameter_integrand = lambda z: _angular_diameter_inner_function(
        z, omega_matter, omega_lambda)

    d_A_z = c.c / H_0 * 1/(1+z) *  quad(angular_diameter_integrand, 0, z)[0]

    return d_A_z.to(u.Mpc)

def make_plot_4ab():

	fig = plt.figure()

	return fig

def make_plot_4c():

	fig = plt.figure()

	return fig
