import numpy as np
from constants import *
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def print_output_table(headers,output):
    for i in range(0,len(headers)):
        print("{}".format(headers[i]).ljust(14,' '),end="")
    print("")
    for i in range(0,len(output)):
        for val in output[i]:
            # print("\n%s\t" % val)
            print("{:.3f~P}".format(val).ljust(14,' '),end="")
            # print("{:10.2f}".format(val.magnitude),end="")
        print("")
    
def write_output_table(headers,output,filename):
    f = open(filename,"w")
    for i in range(0,len(headers)):
        f.write("{}".format(headers[i]).ljust(14,' '))
    f.write("\n")
    for i in range(0,len(output)):
        for val in output[i]:
            if(val.magnitude < 0.01):
                f.write("{:.6f~P}".format(val).ljust(14,' '))
            else:
                f.write("{:.4f~P}".format(val).ljust(14,' '))

            # f.write("{:10.2f}".format(val.magnitude),end="")
        f.write("\n")
    f.close()

def trajectory_output_plots():
    
    fig, ax_acc = plt.subplots()
    # ax_vel = ax_acc.twinx()
    ax_alt = ax_acc.twinx()
    ax_acc.plot(time_p,acceleration_p,color='b')
    ax_acc.plot(time_p,velocity_p,color='g')
    ax_alt.plot(time_p,altitude_p,color='r')
    fig.legend(["acceleration [ft/s**2]","velocity [ft/s]","altitude [m]"])
    plt.show()
    
def PR09_output_plots(ri, rf, L0, ngrain, grain_spacing):
    # plot 4 main plots
    units.setup_matplotlib()
    fig, ax = plt.subplots()
    # plt.rcParams["figure.figsize"] = (12.5,11)
    matplotlib.rcParams.update({'font.size': 14})

    # use rectangles to create cross section view
    # plt.subplot(2, 2, 1)
    total_propellant_length = L0 * 4 + grain_spacing * (4 - 1)
    single_propellant_length = L0 #total_propellant_length / ngrain
    initial_grain_width = rf - ri
    rect = patches.Rectangle((0,-rf),total_propellant_length,rf*2,linewidth=1,edgecolor='b',facecolor="none")

    # for each grain (top and bottom)
    for i in range(0,ngrain):
        for j in range(0,2):
            x = i * (single_propellant_length + grain_spacing)
            
            # "top"
            grain = patches.Rectangle((x,rf - initial_grain_width),single_propellant_length,initial_grain_width,linewidth=1,edgecolor='r',facecolor="none")
            plt.gca().add_patch(grain)

            # "bottom"
            grain = patches.Rectangle((x,-rf),single_propellant_length,initial_grain_width,linewidth=1,edgecolor='r',facecolor="none")
            plt.gca().add_patch(grain)

    plt.gca().add_patch(rect)

    plt.axis('equal')
    plt.xlabel("radial location - $r$ [$in$]")
    plt.ylabel("axial location - $z$ [$in$]")

    # plt.subplot(2, 2, 2)
    # plt.plot(time_p, mp_p)
    # plt.xlabel("time - $t$ [$s$]")
    # plt.ylabel("mass propellant - $m_p$ [$lb_m$]")

    # plt.subplot(2, 2, 3)
    # plt.plot(time_p, p1_p)
    # plt.xlabel("time - $t$ [$s$]")
    # plt.ylabel("chamber pressure - $p_1$ [$psi$]")

    # plt.subplot(2, 2, 4)
    # plt.plot(time_p, F_p)
    # plt.xlabel("time - $t$ [$s$]")
    # plt.ylabel("thrust - $F$ [$lb_f$]")
    return plt