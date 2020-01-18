import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import time
import os

# custom functions
from nozzle_calcs import CF2, Plot3_7_PR09, CFmin
from constants import *
from gas import Gas
from rocket_io import *

# An internal balLstics program that wll calculate the steady-state chamber pressure, mass of
#   propellant, and thrust of the soLd rocket as a function of time and altitude (Including Nozzle
#   Erosion). 

# area of annulus for some inner radius r and outer radius ro
def area_of_annulus(r_inner, r_outer):
    return np.pi * (r_outer ** 2 - r_inner ** 2)

# area of burn surface as a function of web (for N cyLndrcal soLd grains wth single perforation)
def burn_surface_area(L, r_inner, r_outer, N):
    return N * 2 * (np.pi * r_inner * L + area_of_annulus(r_inner,r_outer))

# mass of propellant as a function of web
def mass_of_propellant(L, r, rf, rho_p, N):
    return N * rho_p * area_of_annulus(r,rf) * L

# steady-state chamber pressure including temperature sensitivity terms
def chamber_pressure(A_b, A_t, rho_p, T_bi, T_b0, a0, sigma_p, c_star, n):
    return (a0 * rho_p * c_star * np.exp(sigma_p * (T_bi - T_b0)) * A_b / A_t / lbfToLbm) ** (1 / (1 - n))

# propellant burning rate including temperature sensitivity terms
def burning_rate(T_bi, T_b0, a0, sigma_p, p1, n):
    return a0 * np.exp(sigma_p * (T_bi - T_b0)) * p1 ** n

# equations used for determining the throat area and throat erosion. 
def calc_w_max(r0,rf,L0):
    return min(rf - r0, L0/2)

#CD from Mach from figure 4.3 in Sutton 9th Edition
def drag_coefficient_from_mach(M):

    # strip units b/c correlation
    M = M.magnitude

    if(M < 0.6):
        return 0.15
    elif(M < 1.2):
        return -0.12 + 0.45 * M
    elif(M < 1.8):
        return 0.76 - 0.283 * M
    elif(M < 4):
        return 0.311 - 0.034 * M
    else:
        return 0.175

# Computes trajectory from number of grains, area ratio, and ballast weight
def run_rocket_config(ri, rf, L0, ngrain, Ae_At_0, At_0, mass_ballast_0, T_bi):

    # givens/constants
    n = 0.35
    a0 = 0.030 * units["(in / s) * (lbf / in^2) ^ -0.35"]
    sigma_p = 0.001 * units["1/F"]
    c_star = 5210 * units["ft/s"]
    rho_p = 0.065 * units["lb/in^3"]
    T_b0 = 70.0 * units["F"]
    # T_bi = 70.0 * units["F"]
    w_step = (0.01 * units["in"])
    k_exhaust = 1.3
    grain_spacing = 0.125 * units["in"]
    mass_ballast = mass_ballast_0
    mass_structure = 40.0 * units["lb"]
    mass_per_length = 0.25 * units["lb/in"]

    # maxes
    MAX_PRESSURE = 1000.0 * units["psi"]
    MAX_ACCELERATION = 15.0 * units[""]
    MAX_CASE_LENGTH = 34.0 * units["in"]

    # L[i]sts / dynamic vars
    r = []                      # radius of propellant
    L = []                      # length of propellant
    At = []                     # throat area
    Ae_At = []                  # area ratio
    w = []                      # web distance
    time = []                   # time 
    delDt = []                  # change in diameter at the throat
    CF = []                     # coefficient of thrust
    thrust = []                 # thrust
    It = []                     # total impulse
    Is = []                     # specific impulse
    rate = []                   # burn rate
    throat_diameter = []        # diameter at throat
    mass_propellant = []        # mass of propellant
    burn_area = []                 # area of burn
    p1 = []                     # chamber pressure
    ambient_gas = []              # list of gas properties (Gas class)
    Is_compute = []                     # specific impulse

    # pr4 / phyics variables
    altitude = []
    velocity = []
    acceleration = []
    mach = []
    drag_coeff = []
    drag = []
    F_M = []                    # force over mass term for force balance
    D_M = []                    # drag over mass term for force balance
    
    # vehicle characteristics
    vehicle_mass = []
    vehicle_diameter = 6.19 * units("in")
    vehicle_cx_area = np.pi * vehicle_diameter ** 2 / 4.0

    # initial values
    r.append(ri)    
    L.append(L0)
    At.append(At_0)
    Ae_At.append(Ae_At_0)
    w.append(0.0 * units["in"])
    time.append(0.0 * units["s"])
    It.append(0.0 * units["lbf * sec"])
    Is.append(0.0 * units["sec"])
    Is_compute.append(0.0 * units["sec"])
    throat_diameter.append(np.sqrt(4.0 * At[0] / np.pi))
    altitude.append(0 * units("ft"))                      # altitude of H untsville, AL

    # initial physics
    velocity.append(0.0 * units("ft/s"))
    
    # misc initial calculations
    Ae = Ae_At[0] * At[0]
    w_max = calc_w_max(r[0],rf,L[0])
    mass_case = 4.0 * mass_per_length * (L[0] + grain_spacing)

    burned_out = False
    output=[]
    physics_output=[]
    i = 0
    burnout_it = -1
    last_it = -1
    It_total = It[0]

    # pre-run design constraints
    Aport_At0 = np.pi * r[0] ** 2 / At[0]
    case_length = (L0 + grain_spacing) * ngrain 
    if(Aport_At0 < 2.0):
        print("Invalid Configuration: A_p,0 / A_t,0 > 2.0 (%s)" % Aport_At0)
        return 0.0        
    if(case_length > MAX_CASE_LENGTH):
        print("Invalid Configuration: case length > 34.0 in (%s)" % case_length)
        return 0.0

    # loop while propellant isn't burned out
    while(True):

        #time dependent values that rely on previous value
        if(i > 0):
            velocity.append(velocity[i-1] + acceleration[i-1] * (time[i] - time[i-1]))
            altitude.append(altitude[i-1] + (velocity[i] + velocity[i-1]) / 2 * (time[i] - time[i-1]))
            
            if(altitude[i] < 0.0):
                time.pop()
                altitude.pop()
                velocity.pop()
                break
        elif(i > 1000):
            print("Trajectory timeout!")
            break
            
        # update gas state
        ambient_gas.append(Gas("air",altitude[i]))

        #with T[i] calculate Mach number and the CD
        mach.append(velocity[i] / ambient_gas[i].speed_of_sound)
        drag_coeff.append(drag_coefficient_from_mach(mach[i]) * units[""])
        drag.append(0.5 * ambient_gas[i].density * drag_coeff[i] * velocity[i] ** 2 * vehicle_cx_area)

        if(velocity[i-1] < 0):
            time.pop()
            altitude.pop()
            velocity.pop()
            last_it = i
            break
            drag[i] *= -1
        
        # calculate area of current burn surface
        burn_area.append(burn_surface_area(L[i], r[i], rf, ngrain))

        # if propellant is still burning
        if(not burned_out):

            p1.append(chamber_pressure(burn_area[i], At[i], rho_p, T_bi, T_b0, a0, sigma_p, c_star, n))
            rate.append(burning_rate(T_bi, T_b0, a0, sigma_p, p1[i], n))

            # calculate force stuff
            # print("%s / %s = %s" % (p1[i], ambient_gas[i].pressure, p1[i] / ambient_gas[i].pressure))
            CF.append(CF2(k_exhaust,Ae_At[i], p1[i] / ambient_gas[i].pressure))
            thrust.append(CF[i] * p1[i] * At[i])
            mass_propellant.append(mass_of_propellant(L[i], r[i], rf, rho_p, ngrain))
            
            # after one iteration, we can calculate impulse
            if(i > 0):
                It.append(0.5*(thrust[i] + thrust[i - 1])*(time[i] - time[i - 1]))
                Is.append(It[i] * lbfToLbm / ((mass_propellant[i - 1] - mass_propellant[i]) * g0_EN))

            # check this or next step will burnout, if so, do half step instead
            if(w[i] >= w_max):
                burned_out = True
            else:

                w.append(w[i] + w_step)         # current web distance
                L.append(L[0] - 2 * w[i + 1])   # current grain length (decreases w[i]th burn)
                r.append(r[0] + w[i + 1])       # current radius (increases w[i]th burn)

                # calculate new time step from burn rate
                time.append(time[i] + w_step / rate[i]) 

                # calculate change in diameter at throat
                delDt.append(0.000087 * units["in ^ 3 / s / lbf"] * p1[i] * (time[i + 1] - time[i]))

                if(w[i + 1] + w_step > w_max):
                    w_step = 0.005 * units["in"]
                    burnout_it = i


        # else if propellant is all burned up
        else:
            It.append(0.0 * It[0])
            Is.append(0.0 * Is[0])
            p1.append(0.0 * p1[0])
            rate.append(0.0 * rate[0])
            CF.append(0.0 * CF[0])
            thrust.append(0.0 * thrust[0])

        # same as else above, but we want it to happen the first time burn is turned off
        if(burned_out):  

            mass_propellant.append(0 * mass_propellant[0])

            w.append(w[i])
            L.append(L[0] - 2 * w[i + 1])   # current grain length (decreases w[i]th burn)
            r.append(r[0] + w[i + 1])       # current radius (increases w[i]th burn)

            # switch to 0.1 time step
            time.append(time[i] + 0.1 * units["s"])
            delDt.append(0.0 * delDt[0])

        vehicle_mass.append(mass_ballast + mass_propellant[i] + mass_structure + mass_case)

        F_M.append(thrust[i] / vehicle_mass[i] * lbfToLbm)
        D_M.append(drag[i] / vehicle_mass[i])
        acceleration.append(F_M[i] - D_M[i] - g0_EN)

        ## New Step ## 
        # The following calcs had initial value and utilize the above values #
        # i.e. a lot of [i + 1]'s #

        # get new area at throat
        throat_diameter.append(throat_diameter[i] + delDt[i])
        At.append(np.pi * (throat_diameter[i + 1] / 2) ** 2)
        Ae_At.append(Ae / At[i + 1])
        # print("%s\t%s\t%s\t%s\n" % (time[i],vehicle_mass[i], w[i], velocity[i]))

        It_total += It[i]

        # check if iteration violates any condition
        if(acceleration[i] / g0_EN >= MAX_ACCELERATION):
            print("Invalid Configuration: exceeded max acceleration")
            print("--> %s >= %s" % (acceleration[i] / g0_EN, MAX_ACCELERATION))
            return 0.0
        if(p1[i] > MAX_PRESSURE):
            print("Invalid Configuration: exceeded max pressure")
            return 0.0
        if(CF[i] < CFmin(Ae_At[i]) and not burned_out):
            print("Invalid Configuration: below min CF")
            # pressure_p = np.array([val.magnitude for val in p1]) 
            # print("%f, %f, %f" %(np.min(pressure_p), np.mean(pressure_p), np.max(pressure_p)))

            return 0.0

        # output
        p1[i].ito("psi")
        output.append([w[i],burn_area[i],mass_propellant[i] ,time[i],throat_diameter[i],At[i] 
            ,p1[i] ,CF[i],thrust[i],It[i],rate[i],delDt[i],Ae_At[i]])
        physics_output.append([time[i],vehicle_mass[i],velocity[i],altitude[i],
            ambient_gas[i].pressure,ambient_gas[i].temperature,ambient_gas[i].speed_of_sound,mach[i],
            drag_coeff[i],ambient_gas[i].density,F_M[i],D_M[i],acceleration[i]])

        #increment i
        i += 1

    # for i in range(0,len(It)):
    #     print("%s \t %s \t %s" % (It[i], mass_propellant[i], Is[i]))
    # pressure_p = np.array([val.magnitude for val in p1]) 
    # print(np.mean(pressure_p))
    # return max(altitude).magnitude

    headers = ["web burned", "area", "mass", "time", "Dt ", "At", "P1_e", "C_F", "F_e", "I", "r", "del d", "Ae/At"]
    physics_headers = ["time", "Mveh", "velocity", "altitude", "p3", "temperature", "a", "Mach", "C_d", "rho_3", "F/M", "D/M", "acceleration"]
    # print_output_table(headers,output)
    # print_output_table(physics_headers,physics_output)
    # write_output_table(headers,output,"param_data.txt")
    # write_output_table(physics_headers,physics_output,"physics_data.txt")

    if(burnout_it == -1):
        # print("Warning: no burnout!")
        burnout_it = last_it        

    # convert pint arrays to numpy for plotting/processing
    Is_sub = Is[0:burnout_it]
    Is_p = np.array([val.to_base_units().magnitude for val in Is_sub]) 
    It_p = np.array([val.to_base_units().magnitude for val in It]) 
    time_p = np.array([val.to_base_units().magnitude for val in time]) 
    mp_p = np.array([val.magnitude for val in mass_propellant])
    p1_p = np.array([val.magnitude for val in p1]) 
    F_p = np.array([val.magnitude for val in thrust]) 
    CF_p = np.array([val.magnitude for val in CF]) 
    Ae_At_p = np.array([val.magnitude for val in Ae_At]) 
    acceleration_p = np.array([val.magnitude for val in acceleration]) 
    velocity_p = np.array([val.magnitude for val in velocity]) 
    altitude_p = np.array([val.magnitude for val in altitude]) 

    # Misc post calculations
    Is_ave = (It_total / mass_propellant[0].magnitude)
    max_altitude = max(altitude).magnitude

    config_folder = "results/final/%d" % int(max_altitude)
    if(not os.path.isdir(config_folder)):
        os.mkdir(config_folder)
    matplotlib.rcParams.update({'font.size': 14})
    plt.plot(time_p[:burnout_it],mp_p[:burnout_it])
    plt.xlabel("time [s]")
    plt.ylabel("mass of propellant [kg]")
    plt.savefig("results/final/%d/mass_propellant.png" % int(max_altitude))
    plt.clf()

    plt.plot(time_p[:burnout_it],p1_p[:burnout_it])
    plt.xlabel("time [s]")
    plt.ylabel("chamber pressure [kg]")
    plt.savefig("results/final/%d/chamber_pressure.png" % int(max_altitude))
    plt.clf()

    plt.plot(time_p[:burnout_it],F_p[:burnout_it])
    plt.xlabel("time [s]")
    plt.ylabel("thrust [kg]")
    plt.savefig("results/final/%d/thrust.png" % int(max_altitude))
    plt.clf()

    plt.plot(time_p,acceleration_p)
    plt.xlabel("time [s]")
    plt.ylabel("acceleration [ft/s^2]")
    plt.savefig("results/final/%d/acceleration.png" % int(max_altitude))
    plt.clf()

    plt.plot(time_p,velocity_p)
    plt.xlabel("time [s]")
    plt.ylabel("velocity [ft/s]")
    plt.savefig("results/final/%d/velocity.png" % int(max_altitude))
    plt.clf()

    plt.plot(time_p,altitude_p)
    plt.xlabel("time [s]")
    plt.ylabel("altitude [ft]")
    plt.savefig("results/final/%d/altitude.png" % int(max_altitude))
    plt.clf()

    Plot3_7_PR09(Ae_At_p[:burnout_it],CF_p[:burnout_it]).savefig("results/final/%d/thrust_coefficient.png" % int(max_altitude))
    plt.clf()
    # print("")
    # print("It_total: %s" % np.sum(It_p))
    # print("Is_ave: %s" % Is_ave)
    # print("P1_max %s" % max(p1))
    # print("F_max %s" % max(thrust))
    # print("Aport / At_0: %s" % Aport_At0)

    # print("A_max %s" % max(altitude))
    burnout_it = burnout_it - 1
    print("%.3f | %.3f | %.3f | %.3f |%.3f | %.3f | %.3f" % (mp_p[0],burn_area[0].magnitude,time[burnout_it].magnitude,
            Is_ave.magnitude, max(p1).magnitude,max(thrust).magnitude, Aport_At0.magnitude))
    print("%.3f | %.3f | %.3f | %.3f |%.3f | %.3f |%.3f | %.3f | %.3f" % (mass_case.magnitude, 
        vehicle_mass[0].magnitude, altitude[burnout_it].magnitude, velocity[burnout_it].magnitude, acceleration[burnout_it].magnitude,
        max_altitude, max(velocity).magnitude, max(acceleration).magnitude, max(acceleration).magnitude / g0_EN.magnitude))
    # print(" = %f ft" % max_altitude)
    return (mp_p[0],burn_area[0].magnitude,time[burnout_it].magnitude,
            Is_ave.magnitude, max(p1).magnitude,max(thrust).magnitude, Aport_At0.magnitude, 
            max_altitude, max(velocity).magnitude, max(acceleration).magnitude)
