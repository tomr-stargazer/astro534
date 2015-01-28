""" 
A script to create plots for problem 4: 
Hubble Time and angular diameter distance for various cosmologies,
from z=0 to z=10.

Astro 534 Winter 2015 U. Mich. 
Written by Tom Rice

"""

from __future__ import division

import numpy as np
import astropy.units as u
import astropy.constants as c

from scipy.integrate import quad


def hubble_time_integrand(z, omega_matter, omega_lambda):

	return 1 / ( (1+z) * (omega_matter*(1+z)**3 + omega_lambda)**(1/2) )

def hubble_time_general(redshift, hubble_constant, omega_matter, omega_lambda):
	"""

	"""

	z = redshift
	H_0 = u.Quantity(hubble_constant, u.km / u.s / u.Mpc)

	t_H_z = quad(hubble_time_integrand, z, np.inf)

	return t_H_z