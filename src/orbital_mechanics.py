from constants import *
import numpy as np

# two body motion for small orbiting object
def acceleration(mu, r):
    return -mu * r / np.magnitude(r) ** 2.0

# potential energy of orbit relative to earth
def potential_energy(mu, r, mass):
    return -mu * mass / r

# total specific energy (epsilon)
def total_specific_energy(mu, r, v):
    return v ** 2 / 2.0 - mu / r

# total specific energy (epsilon)
def total_specific_energy(mu, a):
    return - mu / (2 * a)

# specific angular momentum given radius and velocity vector
def specific_angular_momentum(r, v):
    return np.cross(r,v)

# specific angular momentum given radius, velocity, and flight path angle
def specific_angular_momentum(r, v, phi):
    return r * v * np.cos(phi)

# radius in two body problem at
def radius(phi, e, theta):
    return p / (1 + e * np.cos(theta))

# mass ratio
def mass_ratio(delv, Isp, g0 = g0_SI):
    return np.exp(delv/(g0 * Isp))

# inert mass fraction
def inert_mass_fraction(mass_inert, mass_prop):
    return mass_inert / (mass_inert + mass_prop)

# propellant mass fraction
def propellant_mass_fraction(mass_prop, mass_inert):
    return mass_prop / (mass_prop + mass_inert)

def mass_propellant(mass_payload, f_inert, mr):
    return mass_payload * (1 - f_inert) * (mr - 1) / (1 - f_inert * mr)

# mass propellant
# def mass_propellant(delv, Isp, mass_payload, f_inert, g0 = g0_SI):
#     mr = mass_ratio(delv, Isp, g0)
#     return mass_propellant(mass_payload, f_inert, mr)

# compute imf, mr from common starting known masses
def generic_mass_sizing(mass_initial, mass_prop, mass_payload):
    # calcs
    mass_inert = mass_initial - mass_prop - mass_payload
    imf = intert_mass_fraction(mass_inert,mass_prop)
    mr = mass_initial / (mass_initial - mass_prop)

    return imf, mr, mass_inert

