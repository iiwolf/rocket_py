from pint import UnitRegistry

#L[i]brary to assign units to var[i]ables
units = UnitRegistry()
# units.default_system = 'imperial'
# units.auto_reduce_dimensions = True

#custom definition of slug since it's not included
units.define("slug = 1 * pound * s**2 / ft = slugs")
units.define("lbmol = 453.59237 * mol")
units.define("coulomb = 1 ampere * second")
units.define("tesla = kg / (s ** 2 * ampere)")
units.define("gauss = tesla * 1e-4")

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

mu_e = 3.986e5 * units["km^3 / s^2"]
r_e = 6378 * units["km"]
q_electron = 1.6e-19 * units["coulomb"] 
mass_electron = 9.10938356e-31 * units["kg"]
boltzmann_constant = 1.38064852e-23 * units["m ** 2 * kg / (s ** 2 * K)"]
plancks_constant = 6.62607004e-34 * units["m ** 2 * kg / s"]

# v_avg = 4.0 * pi * (m / (2 * pi * k * T)) ^ (3/2) * v ^ 2 * exp(-m * v ^ 2 / (2 * k * T))
# v_dv = -(np.sqrt(2)*(m/(k*T))**(3/2)*v*(m*v**2-2*k*T)*e**(-(m*v**2)/(2*k*T)))/(np.sqrt(pi)*k*T)
