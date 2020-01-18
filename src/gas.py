from constants import *
import numpy as np

#pressure (psia) from altitude relationship from Appendix 2 Sutton 9th Edition
def pressure_at_alt(alt):
    # strip units b/c correlation
    alt = alt.magnitude
    if(alt < 83000):
        return (-4.272981E-14*alt**3 + 0.000000008060081*alt**2 
            - 0.0005482655*alt + 14.69241)
    else:
        return 0.0

#temperature(R) from altitude relationship from Appendix 2 Sutton 9th Edition
def temperature_at_alt(alt):
    # strip units b/c correlation
    alt = alt.magnitude
    if(alt < 32809):
        return (-0.0036*alt + 518)
    else:
        return 399

#density(lbm/ft^3) from altitude relationship from Appendix 2 Sutton 9th Edition
def density_at_alt(alt):
    # strip units b/c correlation
    alt = alt.magnitude
    if(alt < 82000):
        return ((0.00000000001255)*alt ** 2 - (0.0000019453)*alt
                + 0.07579)
    else:
        return 0

class Gas:

    # initialize gas with name
    def __init__(self,name,altitude):
        self.name = name
        self.pressure = pressure_at_alt(altitude) * units("psi") 
        self.temperature = temperature_at_alt(altitude) * units("R")
        self.density = density_at_alt(altitude) * units("lb/ft^3")

        if(name == "air"):
            self.k = k_air
            self.speed_of_sound = np.sqrt(self.k * R_air_EN * self.temperature)
        else:
            print("Warning: k not set for gas \"%s\"!" % name)
            self.k = 0.0 * units("")

        # convert to current project units
        self.speed_of_sound.ito("ft/s")


