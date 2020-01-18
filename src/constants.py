from pint import UnitRegistry

#L[i]brary to assign units to var[i]ables
units = UnitRegistry()
# units.default_system = 'imperial'
# units.auto_reduce_dimensions = True

#custom definition of slug since it's not included
units.define("slug = 1 * pound * s**2 / ft = slugs")
units.define("lbmol = 453.59237 * mol")

# conversions
lbfToLbm = 32.2 * units["lb * ft / (lbf * s^2)"]

#misc constants
g0_EN = 32.2 * units["ft / s**2"]
g0_SI = 9.81 * units["m/s ** 2"]
p_atm = 14.69241 * units["psi"]
R_universal_EN = 1544.0 * units["lbf * ft/(lbmol * R)"]
R_universal_SI = 8.3144598 * units["J / mol * 1 / K"]
MW_air_SI = 28.97 * units["g / (kmol)"]
MW_air_EN = 28.97 * units["lb / (lbmol)"]
R_air_SI = 287.05 * units["J * 1 / kg * 1 / K"]
R_air_EN = R_universal_EN / MW_air_EN * lbfToLbm 
k_air = 1.4 * units[""]
density_water_SI = 997.0 * units["kg/m^3"]