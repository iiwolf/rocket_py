import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import time
import itertools

from constants import *
from rocket_calc import run_rocket_config
from rocket_io import PR09_output_plots

MAX_FINAL_RADIUS = 2.375
MIN_FINAL_RADIUS = 1.000
MAX_INITIAL_RADIUS = 1.8
MIN_INITIAL_RADIUS = 0.5
MAX_LENGTH = (34.0 - 0.125 * 3.0) / 4.0
MIN_LENGTH = 2.0    
MIN_THROAT_AREA = 0.1
MAX_THROAT_AREA = 4.0
MIN_AREA_RATIO = 1.0
MAX_AREA_RATIO = 8.0
MIN_BALLAST_WEIGHT = 0.0
MAX_BALLAST_WEIGHT = 1.0

# check if altitdue falls within any target altitudes to a certain tolerance
def is_valid_altitude(altitude, target_altitudes, tol):
    for alt in target_altitudes:
        if(altitude >= alt - tol and altitude <= alt + tol):
            return True
    return False

def define_design_space():
    print("Defining design space...")
    perm = itertools.permutations([MIN_INITIAL_RADIUS, MAX_INITIAL_RADIUS, MIN_FINAL_RADIUS, MAX_FINAL_RADIUS,
            MIN_LENGTH, MAX_LENGTH, MIN_THROAT_AREA, MAX_THROAT_AREA, MIN_AREA_RATIO, 
            MAX_AREA_RATIO, MIN_BALLAST_WEIGHT, MAX_BALLAST_WEIGHT])
    print(len(list(perm)))

# run range of rocket configs and dump data
def raw_rocket_runs():
    
    # target altitdues and tolerance
    target_altitudes = [5000.00, 10000.00, 15000.00]
    tol = 500.00

    # variable parameters
    ngrains = np.arange(8,30)
    throat_areas = np.arange(1.0, 1.51, 0.1)
    area_ratios = np.arange(2.30, 2.81, 0.1)
    mass_ballasts = np.arange(1.0, 1.1, 0.5)
    initial_radii = np.arange(1.20, 1.51, 0.1)
    final_radii = np.arange(1.70, 2.61, 0.1)
    lengths = np.arange(2.0,6.1,2.0)
    it = 0    

    # configurations = np.zeros(shape=(num_iterations,5))
    f = open("output_0.dat","w")
    for throat_area in throat_areas:
        for area_ratio in area_ratios:
            for ri in initial_radii:
                for rf in final_radii:
                    for L0 in lengths:
                        for mass_ballast in mass_ballasts:
                            valid_count = 0
                            for ngrain in ngrains:
                                it+=1
                                try:
                                    max_altitude = run_rocket_config(ri * units["in"], rf * units["in"], L0 * units["in"],
                                                ngrain, area_ratio * units[""], throat_area * units["in^2"], mass_ballast * units["lb"])

                                    print("\t%f %f %f %f %f %f %f %f\n" % (ri, rf, L0, ngrain, area_ratio, throat_area, mass_ballast, max_altitude))

                                    # plt.scatter(ngrain,max_altitude)
                                    # plt.draw()
                                    # plt.pause(0.01)
                                    if(is_valid_altitude(max_altitude, target_altitudes, tol)):
                                        if(valid_count == 0):               
                                            print("Nozzle Configuration: A_t = %f \tAb_At = %f\n" % (throat_area, area_ratio))
                                            f.write("Nozzle Configuration: A_t = %f \tAb_At = %f\n" % (throat_area, area_ratio))
                                        elif(valid_count == 3):
                                            print("VALID DESIGN!")
                                            f.write("VALID DESIGN!")
                                            
                                        f.write("\t%f %f %f %f %f %f %f %f\n" % (ri, rf, L0, ngrain, area_ratio, throat_area, mass_ballast, max_altitude))
                                        valid_count+=1

                                    if(it % 1000 == 0):
                                        print("%d" % it)
                                    elif(it % 10000 == 0):
                                        f.close()
                                        f = open("output_%d.dat" % it)

                                except Exception as e:
                                    print("Error: " + str(e))
                                    pass
                
    f.close()                            

def plot_rocket_runs():

    configurations = []
    legend_labels = []
    i = 0
           
    # variable parameters
    ngrains = np.arange(1,5)
    throat_areas = np.arange(1.3, 1.31, 0.3)
    area_ratios = np.arange(2.30, 2.31, 0.3)
    mass_ballasts = np.arange(0.00, 1.0, 0.1)
    initial_radii = np.arange(1.1, 1.11, 0.1)
    final_radii = np.arange(1.75, 1.751, 0.25)
    lengths = np.arange(6.0,6.1,2.0)

    # configurations = np.zeros(shape=(num_iterations,5))
    for throat_area in throat_areas:
        # plt.clf()
        color_count = 0
        for area_ratio in area_ratios:
            configurations = []
            for ri in initial_radii:
                for rf in final_radii:
                    for L0 in lengths:
                        for mass_ballast in mass_ballasts:
                            for ngrain in ngrains:

                                max_altitude = run_rocket_config(ri * units["in"], rf * units["in"], L0 * units["in"],
                                            ngrain, area_ratio * units[""], throat_area * units["in^2"], mass_ballast * units["lb"])
                                print("%f, %f, %f, %f, %f, %f, %f, %f" % (ri, rf, L0, ngrain, area_ratio, throat_area, mass_ballast, max_altitude))

                                if(max_altitude != 0.0):
                                    configurations.append(np.array([ngrain, throat_area, area_ratio, mass_ballast, max_altitude]))
                                    color = "C" + str(ngrain)
                                    i += 1

            configs = np.array(configurations)
            # legend_labels.append("ngrain: %d" % ngrain)
            if(configs.size > 2):
                color_count+=1
                plt.plot(configs[:,0],configs[:,4])
                plt.legend(legend_labels)

                plt.title("At %f" % throat_area)
                plt.xlabel("ngrain [-]")
                plt.ylabel("altitude [ft]")
                plt.savefig("results/config_throat_area_%f.png" % throat_area)
                plt.draw()
                plt.pause(0.01)

def plot_initial_temp_comparisons():
    temps = np.arange(30,120,10)
    results = []

    for temp in temps:
        # temperature sensitivity
        results.append(run_rocket_config(1.0 * units["in"], 1.75 * units["in"], 6.0 * units["in"],
                        3.0, 2.0 * units[""], 2.165 * units["in^2"], 0.00492999 * units["lb"], 
                        temp * units["F"]))

    results = np.array(results)
    alts = results[:,7]
    vel = results[:,8]
    acc = results[:,9]

    plt.plot(temps, alts)
    plt.scatter(temps,alts)
    plt.xlabel("initial propellant temperature [F]")
    plt.ylabel("altitude [ft]")
    plt.savefig("temp_vs_alt.png")
    plt.clf()

    plt.plot(temps, vel)
    plt.scatter(temps, vel)
    plt.xlabel("initial propellant temperature [F]")
    plt.ylabel("velocity [ft/s]")
    plt.savefig("temp_vs_velocity.png")
    plt.clf()

    plt.plot(temps, alts)
    plt.scatter(temps, alts)
    plt.xlabel("initial propellant temperature [F]")
    plt.ylabel("acceleration [ft/s^2]")
    plt.savefig("temp_vs_acceleration.png")
    plt.clf()
  

def run_all_configs():
    my_vars = np.zeros(shape=10)

    print(my_vars)
    my_vars += run_rocket_config(1.0 * units["in"], 1.75 * units["in"], 6.0 * units["in"],
                    2.0, 2.3 * units[""], 1.3 * units["in^2"], 0.7866799 * units["lb"], 70.0 * units["F"])
    # my_vars += run_rocket_config(1.0 * units["in"], 2.375 * units["in"], 8.0 * units["in"],
                    # 4.0, 4.0 * units[""], 1.0 * units["in^2"], 1.0 * units["lb"], 70.0 * units["F"])
    print(my_vars)
    my_vars += run_rocket_config(1.0 * units["in"], 1.75 * units["in"], 6.0 * units["in"],
                    2.0, 2.3 * units[""], 1.3 * units["in^2"], 0.7866799 * units["lb"], 70.0 * units["F"])
    print(my_vars)
    my_vars += run_rocket_config(1.0 * units["in"], 1.75 * units["in"], 6.0 * units["in"],
                    3.0, 2.3 * units[""], 2.165 * units["in^2"], 0.00492999 * units["lb"], 70.0 * units["F"])
    print(my_vars)
    my_vars += run_rocket_config(1.0 * units["in"], 1.75 * units["in"], 6.0 * units["in"],
                    4.0, 2.0 * units[""], 2.5 * units["in^2"], 0.138799 * units["lb"], 70.0 * units["F"])
    print(my_vars)

def achieve_target_altitudes():
    achieve_target_altitude(1.0 * units["in"], 1.75 * units["in"], 6.0 * units["in"],
                    3.0, 2.0 * units[""], 2.165 * units["in^2"], 0.00492999 * units["lb"], 10000.0)

    achieve_target_altitude(1.0 * units["in"], 1.75 * units["in"], 6.0 * units["in"],
                    4.0, 2.0 * units[""], 2.5 * units["in^2"], 0.138799 * units["lb"], 15000.0)

def achieve_target_altitude(ri, rf, L0, ngrain, Ae_At_0, At_0, mass_ballast_0, T_bi, target_alt):
    diff = 1e6
    step = 0.1
    min_diff = 1e6
    step_count = 0
    tol = 0.01
    while(abs(diff) > tol):
        max_altitude = run_rocket_config(ri, rf, L0 , ngrain, 
                Ae_At_0, At_0, mass_ballast_0, T_bi)
        # print("%s : %f" %(mass_ballast_0, max_altitude))

        diff = max_altitude - target_alt
        # if(abs(diff) > 2000):
            # return diff
        min_diff = min(min_diff, abs(diff))

        if(diff < 0):
            mass_ballast_0 -= step * units["lb"]
            step = step / 10.0
            step_count += 1
            if(step_count > 3):
                break
        else:
            step_count = 0
            mass_ballast_0 += step * units["lb"]

        if(mass_ballast_0.magnitude > 1.0 or mass_ballast_0.magnitude < 0.0):
            print("Invalid configuration: invalid mass ballast")
            break
    
    return min_diff


start = time.time()
p = [1.2 * units["in"], 2.22 * units["in"], 3.5 * units["in"]]


# grain plots
# PR09_output_plots(1.0,2.375,8.0,4,0.125).savefig("baseline_cx.png")
# PR09_output_plots(p[0].magnitude,p[1].magnitude,p[2].magnitude,2,0.125).savefig("5000_cx.png")
# PR09_output_plots(p[0].magnitude,p[1].magnitude,p[2].magnitude,3,0.125).savefig("10000_cx.png")
# PR09_output_plots(p[0].magnitude,p[1].magnitude,p[2].magnitude,4,0.125).savefig("15000_cx.png")
# exit(0)

print("5000")
run_rocket_config(p[0],p[1],p[2], 2.0, 
    2.1 * units[""], 1.0 * units["in^2"], 0.28396 * units["lb"], 70.0 * units["F"])
print("10000")
run_rocket_config(p[0],p[1],p[2], 3.0, 
                        2.1 * units[""], 1.6 * units["in^2"], 0.47683 * units["lb"], 70.0 * units["F"])
print("15000")                        
run_rocket_config(p[0],p[1],p[2], 4.0, 
                        2.1 * units[""], 2.0 * units["in^2"], 0.10344 * units["lb"], 70.0 * units["F"])

exit(0)
print("5000")
achieve_target_altitude(p[0],p[1],p[2], 2.0, 
    2.1 * units[""], 1.0 * units["in^2"], 0.28396 * units["lb"], 70.0 * units["F"], 5000.0)
print("10000")
achieve_target_altitude(p[0],p[1],p[2], 3.0, 
                        2.1 * units[""], 1.6 * units["in^2"], 0.47683 * units["lb"], 70.0 * units["F"], 10000.0)
print("15000")                        
achieve_target_altitude(p[0],p[1],p[2], 4.0, 
                        2.1 * units[""], 2.0 * units["in^2"], 0.10344 * units["lb"], 70.0 * units["F"], 15000.0)

exit(0)
my_vars = np.zeros(shape=10)
throat_areas = np.arange(1.0, 2.00, 0.2)
area_ratios = np.arange(1.5, 4.5, 0.5)
initial_radii = np.arange(1.1, 1.11, 0.1)
final_radii = np.arange(1.5, 2.3, 0.1)
lengths = np.arange(4.0,6.1,0.5)
configs=[]
count = 0
super_min = 1e9
# nozzle config
for throat_area in throat_areas:
    for area_ratio in area_ratios:
        f = open("data/comp2/%f_%f.dat" % (throat_area, area_ratio),"w")
        print("%f_%f.dat" % (throat_area, area_ratio))
        # propellant config
        for ri in initial_radii:
            for rf in final_radii:
                for L0 in lengths:
                    tol = 200.0
                    valid_count = 0
                    min_diff_total = 0.0

                    # for each altitude, try to reach +-200 with only ngrain
                    for i, alt in enumerate([5000.0, 10000.0, 15000.0]):
                        min_diff = 1e6
                        for ngrain in range(2, 10):
                            min_diff = min(min_diff,abs(achieve_target_altitude(ri * units["in"], rf * units["in"], L0 * units["in"],
                                            ngrain, area_ratio * units[""], throat_area * units["in^2"], 0.0 * units["lb"],70.0 * units["F"], alt)))
                            configs.append([ri * units["in"], rf * units["in"], L0 * units["in"],
                                            ngrain * units[""], area_ratio * units[""], throat_area * units["in^2"], 0.0 * units["lb"],70.0 * units["F"], min_diff * units["ft"]])
                            # print("%f" % min_diff)        
                            # print("%s : %f" %(mass_ballast_0, max_altitude))
                            for val in configs[-1]:
                                f.write("{:.4f~P}".format(val).ljust(14,' '))
                                
                            if(abs(min_diff) < 200):
                                # print("Valid design!")
                                valid_count+=1
                                break

                        # print("")
                        min_diff_total += min_diff

                    # print("\n\t\t\t____Total min diff ", min_diff_total)
                    f.write("--%f\n" %(min_diff_total))
                    super_min = min(super_min, min_diff_total)
                    # plt.scatter(ri, min_diff_total)
                    # plt.show()
                    # plt.draw()
                    # plt.pause(0.01)
        f.close()                            
        if(valid_count >=3): 
            print("YEET!")
            for val in configs[-1]:
                print("{:.4f~P}".format(val).ljust(14,' '))

            exit(0)

f.close()
# achieve_target_altitude(propellant[0], propellant[1], propellant[2], 3.0,

#                         2.1 * units[""], 1.5 * units["in^2"], 0.6 * units["lb"], 70.0 * units["F"], 5000.0)
# achieve_target_altitude(1.0 * units["in"], 1.75 * units["in"], 4.3 * units["in"],
#                      3.0, 2.1 * units[""], 1.5 * units["in^2"], 0.0 * units["lb"], 70.0 * units["F"], 10000.0)
# achieve_target_altitude(1.0 * units["in"], 1.75 * units["in"], 4.3 * units["in"],
#                      5.0, 2.1 * units[""], 1.5 * units["in^2"], 0.0 * units["lb"], 70.0 * units["F"], 15000.0)

print(my_vars)
# run_all_configs()

#raw_rocket_runs()
#plot_rocket_runs()
# define_design_space()
end = time.time()
print("Finished: total time elapsed %fs" %(end - start))
exit(0)
max_altitude = run_rocket_config(1.5 * units["in"], 2.375 * units["in"], 8.0 * units["in"],
                1.0, 1.5 * units[""], 2.7 * units["in^2"], 1.0 * units["lb"])
print(max_altitude)
max_altitude = run_rocket_config(1.5 * units["in"], 2.375 * units["in"], 8.0 * units["in"],
                2.0, 1.5 * units[""], 2.7 * units["in^2"], 1.0 * units["lb"])
print(max_altitude)
max_altitude = run_rocket_config(1.5 * units["in"], 2.375 * units["in"], 8.0 * units["in"],
                3.0, 1.5 * units[""], 2.7 * units["in^2"], 1.0 * units["lb"])
print(max_altitude)
