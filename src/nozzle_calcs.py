#propulsion suite
import matplotlib
from matplotlib.ticker import ScalarFormatter
import matplotlib.pyplot as plt
import numpy as np
# matplotlib.use('tkagg')

# A2/At Ratio Function
def A2At(M2, k):
    return  (1.0 / M2) \
            * (2.0 / (k + 1.0) \
            * (1.0 + (k - 1.0)/2.0*(M2**2.0)))**((k + 1.0)/(2.0*(k-1.0)))

# Newton solver to get M2 for given A2At and k
def MachAtExit(k, A2At, M2_Guess):             
    StopCriteria = 0.000001                         #Percent Rel Error to Stop
    EA = StopCriteria * 1.1                         #Doctoring the Stopping Criteria
    M2 = M2_Guess                                   #Exit Mach Number Guess
    i = 0                                           #Setting the Iteration Counter
    while EA > StopCriteria and i < 100: 
        i = i + 1                         
        AFUN = (2 + (k - 1) * M2 ** 2) / (k + 1)
        BFUN = (k + 1) / (2 * (k - 1))
        CFUN = 1 / AFUN
        DFUN = 1 / M2 ** 2
        DERFUN = ((AFUN) ** BFUN) * (CFUN - DFUN)       #Derivative of Mach # Function
        FUNFUN = ((1 / M2) * AFUN ** BFUN) - A2At  #Mach Number Root Function
        Mold = M2                                      #Old Solution
        M2 = M2 - FUNFUN / DERFUN                       #New Solution via Newton Rhapson
        EA = abs((M2 - Mold) / M2) * 100               #Percent Relative Error
        
    return M2                                      

#thrust coefficient (eq 3-30) given k, area, and pressure ratios p2/p1 and p1/p3
def CF(k, A2_At, p2_p1, p1_p3): 
    p3_p1 = 1.0 / p1_p3                #flip p1/p3 ratio
    term1 = 2 * k ** 2 / (k - 1)
    term2 = 2 / (k + 1)
    term3 = (k + 1) / (k - 1)
    term4 = (k - 1) / k
    CF = (term1 * (term2 ** term3) * (1 - (p2_p1 ** term4))) ** 0.5 + (p2_p1 - p3_p1) * A2_At
    return CF

#approximation for minimum CF as prescribed by Dr. Frederick
def CFmin(A2At):
    Cfmin = -0.0445 * np.log(A2At) ** 2 + 0.5324 * np.log(A2At) + 0.1843
    return Cfmin

def TestFunctions():
    k = 1.3
    A2_At = 10
    #given results
    M2 = MachAtExit(k,A2_At,1.5)
    print(M2)
    p2_p1 = (1 + 0.5 * (k - 1) * M2 ** 2) ** (-k / (k - 1))
    print(p2_p1)
    print(CF(k,A2_At,p2_p1,50))

    #comparison with my A2At
    for AreaRatio in np.arange(2,10):
        M2 = MachAtExit(k,AreaRatio,1.5)
        # print(M2)
        # print(A2At(M2,k))

#CF for a given k, A2/At, and p1/p3
def CF2(k,A2_At,p1_p3):
    M2 = MachAtExit(k,A2_At,1.5)
    p2_p1 = (1 + 0.5 * (k - 1) * M2 ** 2) ** (-k / (k - 1))
    return CF(k,A2_At,p2_p1,p1_p3)
    
#T0/T: this term shows up a lot of places (3-12) (only Mt = 1.0?)
def T0_T(k,M):
    return (1 + 0.5 * (k - 1) * M ** 2)

#P0/P: just T0/T to k/k-1
def P0_P(k,M):
    return T0_T(k,M) ** (k / (k - 1))

#Ay/Ax for Mx,My and k (3-14)
def ExpansionAreaRatio(k,Mx,My):
    return (Mx / My) * np.sqrt( (T0_T(k,My)/T0_T(k,Mx) )
    ** ((k + 1)/(k - 1)))


#re-create plot3-1
def Plot3_1():

    #givens
    k1=1.2
    k2=1.3
    machRange = np.arange(0.10,10,0.01)

    #area curves
    areaCurve1 = ExpansionAreaRatio(k1,1.0,machRange)    #Mx = 1 Ax = At
    areaCurve2 = ExpansionAreaRatio(k2,1.0,machRange)    #Mx = 1 Ax = At

    #temperature curves
    tempCurve1 = 1.0/T0_T(k1,machRange)
    tempCurve2 = 1.0/T0_T(k2,machRange)

    #pressure curves
    pressureCurve1 = 1.0/P0_P(k1,machRange)
    pressureCurve2 = 1.0/P0_P(k2,machRange)

    #area plots
    fig, axPT = plt.subplots()

    #common to all axes
    axPT.set_xlim(0.10,10)
    axPT.set_xscale("log")
    axPT.set_xlabel("Mach Number")
    axPT.yaxis.grid(True, which='major')
    axPT.yaxis.grid(True, which='minor')
    axPT.xaxis.grid(True, which='major')
    axPT.xaxis.grid(True, which='minor')
    # axPT.set_xticks([0.1,1,10],('asdsada','1','10'))

    #pressure/temp plots
    axPT.plot(machRange,tempCurve1,linestyle='dashed')
    axPT.plot(machRange,tempCurve2)
    axPT.plot(machRange,pressureCurve1,linestyle='dashed',color='C0')
    axPT.plot(machRange,pressureCurve2,color='C1')
    axPT.set_ylim(0.01,2.0)
    # axPT.set_yticks([0.01,0.10,1.0],('0.01','0.10','1.0'))
    axPT.set_ylabel("Pressure ratio p/p0 and temperature ratio T/T0\n\n")
    axPT.yaxis.set_label_position("left")
    axPT.set_yscale("log")

    axArea = axPT.twinx()
    axArea.plot(machRange,areaCurve1,linestyle='dashed')
    axArea.plot(machRange,areaCurve2,)
    axArea.yaxis.set_label_position("right")
    axArea.set_ylim(1.0,500)
    axArea.set_yticks([1.0,10,100,500])
    axArea.set_ylabel("\n\nArea Ratio A/At")
    axArea.set_yscale("log")
    

    plt.title("Figure 3-1")
    # plt.grid(b=True, which='major', color='#999999', linestyle='-', alpha=0.2)
    # plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
    plt.legend(("k=1.2","k=1.3"))
    plt.show()



#re-create plot3-7
def Plot3_7():
    k = 1.3
    p1_p3 = [5,10,20,50,100,200,500,1000,2000,5000]
    A2_At = np.arange(1,1000,0.1)
    xlabels=[1,3,10,30,100,300,1000]
    it = 9
    fig = plt.figure()
    # ax = fig.add_axes(xlabels)
    # ax.set_ticks(xlabels)
    for p1p3 in p1_p3:

        it = it + 1
        cfArr = []
        areaArr = []
        for AR in A2_At:

            cfVal = CF2(k,AR,p1p3)

            if(cfVal < CFmin(AR)):
                break

            M2 = MachAtExit(k,AR,1.5)
            p2_p1 = (1 + 0.5 * (k - 1) * M2 ** 2) ** (-k / (k - 1))

            cfArr.append(cfVal)
            areaArr.append(AR)
        
        # print("AR %f" % AR)
        if(not cfArr):
            # print("Empty list!")
            cfArr = [0]*5
        if(len(cfArr) < 5):
            lastVal = cfArr[-1]
            for i in range(len(cfArr),5):
                cfArr.append(lastVal)
        # if((AR < 5.0 or it/10 in (10,20,50,100,200,400,800,1000))):
        #     print("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f" % (AR,M2,p2_p1,cfArr[0],cfArr[1],cfArr[2],cfArr[3],cfArr[4]))
        plt.plot(areaArr,cfArr)

    # exit
    plt.plot(A2_At,CFmin(A2_At))
    plt.xlim(1,1000)
    plt.ylim(0.6,1.9)
    plt.yticks(np.arange(0.6, 2.3, step=0.2))
    plt.xlabel("Area Ratio (A2/At)")
    plt.ylabel("CF")
    plt.title("Thrust Coefficient for k=1.3")
    plt.xscale("log")
    plt.grid()
    plt.legend(p1_p3)
    plt.show()


#re-create plot3-7 with pr09 addition
def Plot3_7_PR09(Ae_At,CF):
    k = 1.3
    p1_p3 = [5,10,20,50,100,200,500,1000,2000,5000]
    A2_At = np.arange(1,1000,0.1)
    xlabels=[1,3,10,30,100,300,1000]
    it = 9
    fig = plt.figure()

    for p1p3 in p1_p3:

        it = it + 1
        cfArr = []
        areaArr = []
        for AR in A2_At:

            cfVal = CF2(k,AR,p1p3)

            if(cfVal < CFmin(AR)):
                break

            M2 = MachAtExit(k,AR,1.5)
            p2_p1 = (1 + 0.5 * (k - 1) * M2 ** 2) ** (-k / (k - 1))

            cfArr.append(cfVal)
            areaArr.append(AR)
        
        # print("AR %f" % AR)
        if(not cfArr):
            # print("Empty list!")
            cfArr = [0]*5
        if(len(cfArr) < 5):
            lastVal = cfArr[-1]
            for i in range(len(cfArr),5):
                cfArr.append(lastVal)
        plt.plot(areaArr,cfArr)

    # exit
    plt.plot(Ae_At,CF,linewidth=6.0,c="black")
    plt.plot(A2_At,CFmin(A2_At))
    plt.xlim(1,1000)
    plt.ylim(0.6,1.9)
    plt.yticks(np.arange(0.6, 2.3, step=0.2))
    plt.xlabel("area ratio [-]")
    plt.ylabel("coefficient of thrust [-]")
    plt.xscale("log")
    plt.grid()
    plt.legend(p1_p3, loc="right")
    return plt
# Plot3_1()
# Plot3_7()