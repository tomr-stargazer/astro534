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

from hw1_2c import hubble_time


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


def angular_diameter_distance(redshift, hubble_constant, omega_matter, omega_lambda, redshift_of_observer=0):
    """
    Returns angular diameter distance at redshift z using H_0, omega_matter, omega_lambda.

    Assumes H_0 given in km s^-1 Mpc^-1.

    """

    z = redshift
    z_0 = redshift_of_observer
    a_0 = 1/(1+z_0)

    H_0 = u.Quantity(hubble_constant, u.km / u.s / u.Mpc)

    angular_diameter_integrand = lambda z: _angular_diameter_inner_function(
        z, omega_matter, omega_lambda)

    d_A_z = a_0 * c.c / H_0 * 1/(1+z) *  quad(angular_diameter_integrand, z_0, z)[0]

    return d_A_z.to(u.Mpc)

def analytic_desitter_distance(redshift, hubble_constant):

    z = redshift
    H_0 = u.Quantity(hubble_constant, u.km / u.s / u.Mpc)

    d_A_z = 2*c.c / H_0 * ( (1+z)**(-1) - (1+z)**(-3/2))

    return d_A_z.to(u.Mpc)


def make_plot_4ab(n_steps=50):

    EdS_cosmology = {'name': 'Einstein deSitter', 'omega_matter': 1, 'omega_lambda': 0, 'H_0': 50}
    LCDM_cosmology = {'name': r'$\Lambda$CDM', 'omega_matter': 0.27, 'omega_lambda': 0.73, 'H_0': 70}

    z_array = np.linspace(0, 10, n_steps)

    fig = plt.figure()

    for cosmology in [EdS_cosmology, LCDM_cosmology]:

        tH_array = np.zeros_like(z_array)

        for i, z in enumerate(z_array):

            tH_array[i] = hubble_time_general(
                z, cosmology['H_0'], cosmology['omega_matter'], 
                cosmology['omega_lambda']).value

        plt.plot(z_array, tH_array, 'o', ms=2, label=cosmology['name'])

    analytic_tH_array = hubble_time(z_array, 50)

    plt.plot(z_array, analytic_tH_array, 'r-', zorder=-1, label='Analytic Einstein deSitter')

    plt.legend()

    plt.xlabel("Redshift $z$", fontsize=18)
    plt.ylabel("Hubble Time $t_H$ (Gyr)", fontsize=18)
    return fig

def make_plot_4c(n_steps=50):

    EdS_cosmology = {'name': 'Einstein deSitter', 'omega_matter': 1, 'omega_lambda': 0, 'H_0': 50}
    LCDM_cosmology = {'name': r'$\Lambda$CDM', 'omega_matter': 0.27, 'omega_lambda': 0.73, 'H_0': 70}

    z_array = np.linspace(0, 10, n_steps)

    fig = plt.figure()

    for cosmology in [EdS_cosmology, LCDM_cosmology]:

        dA_array = np.zeros_like(z_array)

        for i, z in enumerate(z_array):

            dA_array[i] = angular_diameter_distance(
                z, cosmology['H_0'], cosmology['omega_matter'], 
                cosmology['omega_lambda']).value

        plt.plot(z_array, dA_array, 'o', ms=2, label=cosmology['name'])

    analytic_dA_array = analytic_desitter_distance(z_array, 50)

    plt.plot(z_array, analytic_dA_array, 'r-', zorder=-1, label='Analytic Einstein deSitter')

    plt.legend()

    plt.xlabel("Redshift $z$", fontsize=18)
    plt.ylabel("Distance $d_A$ (Mpc)", fontsize=18)

    return fig
