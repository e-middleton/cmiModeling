import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import colors
import pandas as pd
import json
import yaml
import matplotlib as mpl

# global attributes
coast = pd.read_csv("./data/coastline.csv")

### VECTOR STYLING ###
miniVecScale = 500 # mini vec scale used for shorter vectors
arrowWidth = 0.005
arrowWidth2 = 0.003
headLength = 2
headAxisLength=2
headWidth = 3

# colorscaling
colorbarFault = 6 # pos and neg afterslip colorbar limits


def slipDist(fault, cmi, fault_slip, cmi_slip, outputDir) :
    '''Plots the slip distribution for both the cmi alone with gps vectors
    and the fault connected to the cmi without gps vectors'''
    
    # plot for visualizing
    xmin = 122
    xmax = 148
    ymin = 33
    ymax = 45

    

    maxF = np.max(np.abs(fault_slip))
    maxH = np.max(np.abs(cmi_slip))
    max_mag = np.max([maxF, maxH])
    
    both = {}
    both["points"] = np.vstack((fault["points"], cmi["points"]))
    shift_val = len(fault["points"][:,0])
    both["verts"] = np.vstack((fault["verts"], cmi["verts"]+shift_val))

    fig, ax = plt.subplots( figsize=(6,6))
    rso = ax.tripcolor(both["points"][:,0],
                        both["points"][:,1], 
                        both["verts"],
                        facecolors=(np.vstack(((fault_slip, cmi_slip)))).flatten(), 
                        vmin=-max_mag, vmax=max_mag)

    cbar1 = fig.colorbar(rso, ax=ax, orientation='horizontal')
    ax.plot(coast.lon, coast.lat, color="k", linewidth=0.5)
    cbar1.set_label("Slip (m)")
    ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax), aspect='equal')
    ax.title.set_text("Dip Slip") #graph 1
    ax.set_ylabel("Latitude")
    ax.set_xlabel("Longitude")

    plt.savefig(outputDir + 'images/slipDist.pdf') # save the figure only runs if called by plotFig
    plt.close('all')

    # same plot, but colorscale from cmi slip
    
    fig, ax = plt.subplots( figsize=(6,6))
    rso = ax.tripcolor(both["points"][:,0],
                        both["points"][:,1], 
                        both["verts"],
                        facecolors=(np.vstack(((fault_slip, cmi_slip)))).flatten(), 
                        vmin=-np.max(np.abs(cmi_slip)), vmax=np.max(np.abs(cmi_slip))) # colorbar based on cmi slip 
    cbar1 = fig.colorbar(rso, ax=ax, orientation='horizontal')
    ax.plot(coast.lon, coast.lat, color="k", linewidth=0.5)
    cbar1.set_label("Slip (m)")
    ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax+2), aspect='equal')
    ax.title.set_text("Dip Slip") #graph 2
    ax.set_ylabel("Latitude")
    ax.set_xlabel("Longitude")

    plt.savefig(outputDir + 'images/slipDistCmi.pdf') # save the figure only runs if called by plotFig
    plt.close('all')

    return


def afterslip(fault_slip_mag, fault, outputDir):
    '''plots the afterslip distribution along the subduction zone,
    side by side with a plot showing the location of the peak afterslip'''
    # plot for visualizing
    xmin = 138
    xmax = 148
    ymin = 33
    ymax = 45


    # afterslip total magnitude 
    fig, ax = plt.subplots(figsize=(5,6))
    rso = ax.tripcolor(fault["points"][:,0],
                        fault["points"][:,1], 
                        fault["verts"],
                        facecolors=(fault_slip_mag).flatten(), 
                        vmin=-1*colorbarFault, vmax=colorbarFault)
    cbar1 = fig.colorbar(rso, ax=ax, orientation='vertical')
    ax.plot(coast.lon, coast.lat, color="k", linewidth=0.5)
    cbar1.set_label("Slip (m)")
    ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax), aspect='equal')
    ax.title.set_text("Fault Slip Magnitude") #graph 1
    ax.set_ylabel("Latitude")
    ax.set_xlabel("Longitude")
    plt.savefig(outputDir + "images/afterslip.pdf")
    plt.close('all')


def plotDipSlip(fault, fault_d_slip, outputDir):

    xmin = 138
    xmax = 148
    ymin = 33
    ymax = 45

    # dip slip only
    fig, ax = plt.subplots(figsize=(5,6))
    rso = ax.tripcolor(fault["points"][:,0],
                        fault["points"][:,1], 
                        fault["verts"],
                        facecolors=(fault_d_slip).flatten(), 
                        vmin=-1*colorbarFault, vmax=colorbarFault)
    cbar1 = fig.colorbar(rso, ax=ax, orientation='vertical')
    ax.plot(coast.lon, coast.lat, color="k", linewidth=0.5)
    cbar1.set_label("Dip slip (m)")
    ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax), aspect='equal')
    ax.title.set_text("Fault Dip Slip") #graph 1
    ax.set_ylabel("Latitude")
    ax.set_xlabel("Longitude")
    plt.savefig(outputDir + "images/dipSlip.pdf")
    return
    

def plotSlipLocation(fault, fault_slip_mag, outputDir):
    xmin = 138
    xmax = 148
    ymin = 32
    ymax = 48

    fig, ax = plt.subplots(figsize=(5,6))
    custom_colors = ['lightsteelblue', 'cornflowerblue', 'royalblue', 'mediumblue'] #, 'midnightblue']
    cmap = colors.ListedColormap(custom_colors)

    cmap = mpl.colormaps['Blues'].resampled(4) # 4 discrete bins for colorbar

    # norm = colors.BoundaryNorm(bounds, cmap.N)
    rso2 = ax.tripcolor(fault["points"][:,0], fault["points"][:,1], fault["verts"], facecolors=(fault_slip_mag).flatten(), cmap=cmap, edgecolors='lightgray', linewidth=0.01)
    cbar2 = fig.colorbar(rso2, ax=ax,orientation='vertical')
    cbar2.set_label('slip magnitude (m)')
    ax.plot(coast.lon, coast.lat, color="k", linewidth=0.5)
    ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax), aspect='equal')
    ax.title.set_text("Max Slip Location") #graph 2

    plt.savefig(outputDir + "images/maxSlipLocation.pdf")
    return


def plotDispByComponent(dispMat, allElemBegin, estSlip, gps, outputDir) :
    '''plot displacements, including observed, predicted, and displacements separated
    by the component (i.e., displacement from CMI and displacement from fault)'''

    # calculate displacements by which components
    # calc disp from cmi, beginning from the cmi elements to the end of the cmi elements
    cmi_disp = dispMat[:, allElemBegin[1]:allElemBegin[2]].dot(estSlip[allElemBegin[1]:allElemBegin[2]]) 

    # fault disp
    fault_disp = dispMat[:, allElemBegin[0]:allElemBegin[1]].dot(estSlip[allElemBegin[0]:allElemBegin[1]])

    fig, ax = plt.subplots(1, 2, figsize=(10,5))
    # coastline plotting
    ax[0].plot(coast.lon, coast.lat, color="k", linewidth=0.5)
    ax[1].plot(coast.lon, coast.lat, color="k", linewidth=0.5)

    # plot the component of horizontal displacement from afterslip
    ax[0].quiver(gps.lon, gps.lat, fault_disp[0::3], fault_disp[1::3], 
                 scale=miniVecScale, color='b', label="fault", 
                 headlength=headLength, headaxislength=headAxisLength, headwidth=headWidth, width=arrowWidth, 
                 lw=0.2, ec='k',)
    
    # plot the component of horizontal displacement from VER / cmi slip
    Q2 = ax[0].quiver(gps.lon, gps.lat, cmi_disp[0::3], cmi_disp[1::3], 
                      scale=miniVecScale, color='yellow', label="cmi", 
                      headlength=headLength, headaxislength=headAxisLength, headwidth=headWidth, width=arrowWidth2, 
                      lw=0.2, ec='k',)
    # quiver based off of cmi disp, but same units as afterslip so it should be alright
    ax[0].quiverkey(Q2, X=0.3, Y=0.8, U=30, label="30 cm", labelpos='N', color='b')
    ax[0].set_ylim([34, 42])
    ax[0].set_xlim([136, 144])
    ax[0].set_title("Horizontal contributions")

    # zero vector for the starting vertical position
    ogVec = np.zeros_like(fault_disp[0::3])

    # plot the component of vertical displacement from afterslip
    # note the vector scale is dividied by 5, so shorter displacements appear larger (inverse)
    ax[1].quiver(gps.lon, gps.lat, ogVec, fault_disp[2::3], 
                 scale=miniVecScale/5, color='b', label="fault", 
                 headlength=headLength, headaxislength=headAxisLength, headwidth=headWidth, width=arrowWidth, 
                 lw=0.2, ec='k',)
    
    # plot the component of vertical displacement from VER / cmi
    Q3 = ax[1].quiver(gps.lon, gps.lat, ogVec, cmi_disp[2::3], 
                      scale=miniVecScale/5, color='yellow', label="cmi", 
                      headlength=headLength, headaxislength=headAxisLength, headwidth=headWidth, width=arrowWidth2, 
                      lw=0.2, ec='k',)
    ax[1].quiverkey(Q3, X=0.3, Y=0.8, U=10, label="10 cm", labelpos='N', color='b')
    ax[1].set_ylim([34, 42])
    ax[1].set_xlim([136, 144])
    ax[1].set_title("Vertical contribution")
    ax[1].legend()

    plt.savefig(outputDir + 'images/dispByComponent.pdf')
    plt.close('all')

    fig, ax = plt.subplots( figsize=(6,5))
    ax.plot(coast.lon, coast.lat, color="k", linewidth=0.5)

    # plot the component of horizontal displacement from afterslip against total
    
    # plot the component of horizontal displacement from observed
    Q2 = ax.quiver(gps.lon, gps.lat, fault_disp[0::3], fault_disp[1::3], 
                      scale=miniVecScale*2, color='b', label="fault", 
                      headlength=headLength+2, headaxislength=headAxisLength+1, headwidth=headWidth+2, width=arrowWidth-0.001,
                      minshaft=1.5
                      )
    ax.quiver(gps.lon, gps.lat, cmi_disp[0::3], cmi_disp[1::3], 
                 scale=miniVecScale*2, color='r', label="cmi", 
                 headlength=headLength+2, headaxislength=headAxisLength+1, headwidth=headWidth+2, width=arrowWidth-0.001,
                 minshaft=1.5
                 )
    
    # quiver based off of cmi disp, but same units as afterslip so it should be alright
    ax.quiverkey(Q2, X=0.8, Y=0.1, U=50, label="50 cm", labelpos='N', color='r')
    ax.set_ylim([36, 40])
    ax.set_xlim([138.5, 144.5])
    ax.legend()

    plt.savefig(outputDir + 'images/dispByComponent_agata.pdf')
    plt.close('all')

    fig, ax = plt.subplots(1, 2, figsize=(10,5))
    ax[0].plot(coast.lon, coast.lat, color="k", linewidth=0.5)
    ax[1].plot(coast.lon, coast.lat, color="k", linewidth=0.5)

    # plot the component of horizontal displacement from afterslip against total
    
    # plot the component of horizontal displacement from observed
    Q2 = ax[0].quiver(gps.lon, gps.lat, gps.east_vel, gps.north_vel, 
                      scale=miniVecScale*2, color='r', label="observed", 
                      headlength=headLength, headaxislength=headAxisLength, headwidth=headWidth, width=arrowWidth, 
                      lw=0.2, ec='k',)
    ax[0].quiver(gps.lon, gps.lat, fault_disp[0::3], fault_disp[1::3], 
                 scale=miniVecScale*2, color='cyan', label="fault", 
                 headlength=headLength, headaxislength=headAxisLength, headwidth=headWidth-1, width=arrowWidth-0.0005, 
                 lw=0.2, ec='k',)
    
    # quiver based off of cmi disp, but same units as afterslip so it should be alright
    ax[0].quiverkey(Q2, X=0.3, Y=0.8, U=100, label="100 cm", labelpos='N', color='r')
    ax[0].set_ylim([34, 42])
    ax[0].set_xlim([136, 144])
    ax[0].set_title("Horizontal afterslip contributions")

    Q3 = ax[1].quiver(gps.lon, gps.lat, gps.east_vel, gps.north_vel, 
                      scale=miniVecScale*2, color='r', label="observed", 
                      headlength=headLength, headaxislength=headAxisLength, headwidth=headWidth, width=arrowWidth, 
                      lw=0.2, ec='k',)
    ax[1].quiver(gps.lon, gps.lat, cmi_disp[0::3], cmi_disp[1::3], 
                 scale=miniVecScale*2, color='yellow', label="cmi", 
                 headlength=headLength, headaxislength=headAxisLength, headwidth=headWidth-1, width=arrowWidth-0.0005, 
                 lw=0.2, ec='k',)
    
    # quiver based off of cmi disp, but same units as afterslip so it should be alright
    ax[1].quiverkey(Q3, X=0.3, Y=0.8, U=100, label="100 cm", labelpos='N', color='r')
    ax[1].set_ylim([34, 42])
    ax[1].set_xlim([136, 144])
    ax[1].set_title("Horizontal cmi contributions")

    plt.savefig(outputDir + 'images/dispByComponent_FreedHoriz.pdf')
    plt.close('all')
    return

def plotTotalDisp(predDisp, gps, vecScale, outputDir):
    fig, ax = plt.subplots(1, 2, figsize=(10,5))
    ax[0].plot(coast.lon, coast.lat, color="k", linewidth=0.5)
    ax[1].plot(coast.lon, coast.lat, color="k", linewidth=0.5)
    Q= ax[0].quiver(gps.lon, gps.lat, gps.east_vel, gps.north_vel, scale=vecScale, color='darkgreen', label='observed', headlength=headLength, headaxislength=headAxisLength, headwidth=headWidth, width=arrowWidth, 
                 lw=0.2, ec='k',)
    ax[0].quiver(gps.lon, gps.lat, predDisp[0::3], predDisp[1::3], scale=vecScale, color='orange', label='predicted', headlength=headLength, headaxislength=headAxisLength, headwidth=headWidth-1, width=arrowWidth2, 
                 lw=0.2, ec='k',)
    ax[0].quiverkey(Q, X = 0.3, Y=0.7, U=100, label='100 cm',labelpos='N', color='darkgreen')
    ax[0].set_title("Observed and Predicted Horizontal (cm)")
    ax[0].set_ylim([34, 42])
    ax[0].set_xlim([136, 144])
    
    Q0 = ax[1].quiver(gps.lon, gps.lat, np.zeros_like(gps.north_vel), gps.up_vel, scale=vecScale/10, color='darkgreen', label='observed', headlength=headLength, headaxislength=headAxisLength, headwidth=headWidth, width=arrowWidth, 
                 lw=0.2, ec='k',)
    ax[1].quiver(gps.lon, gps.lat, np.zeros_like(predDisp[0::3]), predDisp[2::3], scale=vecScale/10, color='orange', label='predicted', headlength=headLength, headaxislength=headAxisLength, headwidth=headWidth-1, width=arrowWidth2, 
                 lw=0.2, ec='k',)
    ax[1].quiverkey(Q0, X = 0.3, Y=0.7, U=10, label='10 cm',labelpos='N', color='darkgreen')
    ax[1].set_title("Observed and Predicted Vertical (cm)")
    ax[1].set_ylim([34, 42])
    ax[1].set_xlim([136, 144])
    ax[1].legend()

    plt.savefig(outputDir + 'images/totalDisp.pdf')
    plt.close('all')
    return

# plot the ratio of displacement due to cmi vs total displacements
def plotRatio(gps, totalCmiDisp, totalDisp, saveFigures, ratioFig, outputDir) :
    ratio = totalCmiDisp / totalDisp

    plt.close('all')
    fig, ax = plt.subplots()
    dots = ax.scatter(gps.lon, gps.lat, c=ratio) # color by ratio of cmi_disp / total disp
    ax.plot(coast.lon, coast.lat, color="k", linewidth=0.5)
    ax.set_xlim(130, 145)
    ax.set_ylim(33, 45)
    # Add a customized colorbar
    cbar = fig.colorbar(dots, ax=ax, orientation='vertical', shrink=0.7,
                        label='Ratio of cmi:total', extend='both')

    if saveFigures and ratioFig:
        plt.savefig(outputDir + "images/ratioCmiToTotal.pdf")
    
    plt.close('all')

    return


def residualPlot(gps, predDisp, vecScale, outputDir) :
    # calculate gps displacement residuals
    # residual = actual - predicted

    residuals = np.empty((len(gps.lon), 3))
    actual = np.hstack((np.array(gps.east_vel).reshape(-1,1), np.array(gps.north_vel).reshape(-1,1), np.array(gps.up_vel).reshape(-1,1)))

    # predicted / calculated values
    calc_east = predDisp[0::3].flatten()
    calc_north = predDisp[1::3].flatten()
    calc_up = predDisp[2::3].flatten()

    residuals[:,0] = actual[:,0] - calc_east
    residuals[:,1] = actual[:,1] - calc_north
    residuals[:,2] = actual[:,2] - calc_up

    plt.close('all')
    fig, ax = plt.subplots(1, 2, figsize=(10,5))
    Q = ax[0].quiver(gps.lon, gps.lat, residuals[:,0], residuals[:,1], scale=vecScale/5, color='r', label="residuals", 
                     headlength=headLength, headaxislength=headAxisLength, headwidth=headWidth, width=arrowWidth, 
                      lw=0.1, ec='k',)
    ax[0].quiverkey(Q, X=0.3, Y=0.8, U=20, label="20 cm", labelpos='N', color='r')
    ax[0].set_ylim([34, 42])
    ax[0].set_xlim([136, 144])
    ax[0].plot(coast.lon, coast.lat, color="k", linewidth=0.5)
    ax[0].set_title("residual displacements (horizontal)")

    Q1 = ax[1].quiver(gps.lon, gps.lat, 0, residuals[:,2], scale=vecScale/6, color='r', label="residuals", 
                      headlength=headLength, headaxislength=headAxisLength, headwidth=headWidth, width=arrowWidth, 
                      lw=0.1, ec='k',)
    ax[1].quiverkey(Q1, X=0.3, Y=0.8, U=10, label='10 cm', labelpos='N', color='r')
    ax[1].plot(coast.lon, coast.lat, color="k", linewidth=0.5)
    ax[1].set_ylim([34, 42])
    ax[1].set_xlim([136, 144])
    plt.title("residual displacements (vertical)")

    plt.savefig(outputDir + "images/residuals.pdf")
    plt.close('all')
    return

def plotFigs(config, fault, cmi, gps, estSlip, predDisp, dispMat, vecScale, allElemBegin):
    ''' use config file to determine which figures need to be plotted '''

    outputDir = config["results"]["outputDir"]
    # fault strike slip
    fault_s_slip = np.array(estSlip[0:allElemBegin[1]:2])/100 # cm -> m
    # fault dip slip
    fault_d_slip = np.array(estSlip[1:allElemBegin[1]:2])/100 # cm -> m

    # fault slip magnitude
    fault_slip_mag = np.sqrt(np.square(fault_s_slip) + np.square(fault_d_slip))

    # cmi strike slip 
    cmi_s_slip = np.array(estSlip[allElemBegin[1]::3])/100

    # cmi dip slip     # negated here to match geographic convention for the fault for visualizations, east = pos, west = neg
    cmi_d_slip = -1*np.array(estSlip[1+allElemBegin[1]::3])/100

    # cmi horizontal magnitude 
    cmi_slip_mag = np.sqrt(np.square(cmi_s_slip) + np.square(cmi_d_slip))


    if (config["results"]["saveFigures"]):
        if (config["results"]["afterslip"]):
            afterslip(fault_slip_mag, fault, outputDir)    # afterslip

        if (config["results"]["dipSlip"]):
            plotDipSlip(fault, fault_d_slip, outputDir)  # dip slip

        if (config["results"]["slipLocation"]):
            plotSlipLocation(fault, fault_slip_mag, outputDir) # afterslip location
        
        if (config["results"]["slipDist"]):          # total slip (fault and cmi)
            slipDist(fault, cmi, fault_slip_mag, cmi_slip_mag, outputDir)
        
        if (config["results"]["residuals"]):
            residualPlot(gps, predDisp, vecScale, outputDir)
        
        if (config["results"]["dispByComponent"]):
            plotDispByComponent(dispMat, allElemBegin, estSlip, gps, outputDir)
        
        if (config["results"]["totalDisp"]):
            plotTotalDisp(predDisp, gps, vecScale, outputDir)


def maxSlipMag(fault_slip, cmi_slip):
    ''' gets the maximum magnitude of the total slip magnitude per element of both meshes'''
    maxFaultMag = np.max(np.abs(fault_slip))
    maxCmiMag = np.max(np.abs(cmi_slip))

    return maxFaultMag, maxCmiMag

def calcRMSE(predDisp, gps):
    residuals = np.empty((len(gps.lon), 3))
    actual = np.hstack((np.array(gps.east_vel).reshape(-1,1), np.array(gps.north_vel).reshape(-1,1), np.array(gps.up_vel).reshape(-1,1)))

    # predicted / calculated values
    calc_east = predDisp[0::3].flatten()
    calc_north = predDisp[1::3].flatten()
    calc_up = predDisp[2::3].flatten()

    residuals[:,0] = actual[:,0] - calc_east
    residuals[:,1] = actual[:,1] - calc_north
    residuals[:,2] = actual[:,2] - calc_up

    # calc root mean square residuals
    resid_sqr = np.square(residuals) # square all residuals

    # Calculate the Mean Squared Error (MSE)
    mse = np.mean(resid_sqr)
    
    # Calculate the Root Mean Square Error (RMSE)
    rmse = np.sqrt(mse)

    return rmse

def calcMoment(faultSlip, cmiSlip, fault, cmi):
    '''calculate the moment for slip on fault vs slip on cmi
        moment = rigidity x area x slip'''
    '''calculate the moment for slip on fault vs slip on cmi
        moment = rigidity x area x slip'''

    rigidity = 3e10 #N / m^2 (is this the same though for the CMI at depth?)

    faultMoment = rigidity * (fault["area"]) * faultSlip.flatten() # flatten fault slip to match fault["area"] dimensions
    faultTotalMoment = np.sum(faultMoment) # in Nm

    cmiMoment = rigidity * (cmi["area"]) * cmiSlip.flatten()
    cmiTotalMoment = np.sum(cmiMoment)

    return faultTotalMoment, cmiTotalMoment


def numericalData(estSlip, predDisp, dispMat, gps, allElemBegin, fault, cmi, config):
    saveData = config["results"]["saveData"]
    outputDir = config["results"]["outputDir"]

    # fault strike slip
    fault_s_slip = np.array(estSlip[0:allElemBegin[1]:2])/100 # cm -> m
    # fault dip slip
    fault_d_slip = np.array(estSlip[1:allElemBegin[1]:2])/100 # cm -> m

    # fault slip magnitude
    fault_slip_mag = np.sqrt(np.square(fault_s_slip) + np.square(fault_d_slip))

    # cmi strike slip 
    cmi_s_slip = np.array(estSlip[allElemBegin[1]::3])/100 # cm -> m
    # cmi dip slip     
    cmi_d_slip = np.array(estSlip[1+allElemBegin[1]::3])/100 # cm -> m

    # cmi horizontal magnitude 
    cmi_slip_mag = np.sqrt(np.square(cmi_s_slip) + np.square(cmi_d_slip))

    maxFaultMag, maxCmiMag = maxSlipMag(fault_slip = fault_slip_mag, cmi_slip=cmi_slip_mag)
    rmse = calcRMSE(predDisp, gps)
    faultMoment, cmiMoment = calcMoment(faultSlip=fault_slip_mag, cmiSlip=cmi_slip_mag, fault=fault, cmi=cmi)
    avgFaultRake, avgCmiRake = rakeCalc(fault_d_slip=fault_d_slip, fault_s_slip=fault_s_slip, fault_slip_mag=fault_slip_mag, cmi_s_slip=cmi_s_slip, cmi_d_slip=cmi_d_slip, cmi_slip_mag=cmi_slip_mag)
    percentFault, percentCmi, weightedFault, weightedCmi = contributionsPercent(predDisp=predDisp, dispMat=dispMat, estSlip=estSlip, allElemBegin=allElemBegin, outputDir=outputDir)

    results = {}
    results["faultMaxMag (m)"] = maxFaultMag
    results["cmiMaxMag (m)"] = maxCmiMag
    results["faultMoment"] = faultMoment
    results["cmiMoment"] = cmiMoment
    results["rmse (cm)"] = rmse
    results["avgFaultRake (deg)"] = avgFaultRake
    results["avgCmiRake (deg)"] = avgCmiRake
    results["percentCmi"] = percentCmi
    results["percentFault"] = percentFault
    results["weightedPercentFault"] = weightedFault
    results["weightedPercentCmi"] = weightedCmi    
    results["percentCmi"] = percentCmi
    results["percentFault"] = percentFault
    results["weightedPercentFault"] = weightedFault
    results["weightedPercentCmi"] = weightedCmi

    if (saveData) :
        with open(outputDir + "numericalResults.txt", "w") as file:
            file.write(json.dumps(results, indent=4, sort_keys=True))

    return


# save current config settings for later reference
def saveConfig(config):
    with open(config["results"]["outputDir"] + 'configSettings.yaml', 'w') as file:
        yaml.dump(config, file, default_flow_style=False)
    return


def rakeCalc(fault_d_slip, fault_s_slip, fault_slip_mag, cmi_d_slip, cmi_s_slip, cmi_slip_mag) :
    '''calculates the average rake value for the cmi and fault
    using the slip magnitude to weight the rake value of the given element'''

    # cmiRake = np.rad2deg(np.arctan2(cmi_d_slip.flatten(), cmi_s_slip.flatten())) # (left lateral = 0 rake)
    faultRake = np.rad2deg(np.arctan2(fault_d_slip.flatten(), fault_s_slip.flatten()))
    cmiRake = np.rad2deg(np.arctan2(cmi_d_slip.flatten(), cmi_s_slip.flatten())) # (left lateral = 0 rake)

    weightedCmiRake = np.average(cmiRake, weights=cmi_slip_mag.flatten())
    weightedFaultRake = np.average(faultRake, weights=fault_slip_mag.flatten())

    return weightedFaultRake, weightedCmiRake


def contributionsPercent(predDisp, dispMat, estSlip, allElemBegin, outputDir) :
    ''' function for returning the percentage of total displacement that the cmi is responsible for
    compared to total percentage of displacements that the afterslip is contributing
    the percentages are calculated component (e/n/u) wise because trying to calc percentages 
    using the total displacement vector magnitude compounds rounding differences from np.float64'''

    # calc disp from cmi, beginning from the cmi elements to the end of the cmi elements
    cmi_disp = dispMat[:, allElemBegin[1]:allElemBegin[2]].dot(estSlip[allElemBegin[1]:allElemBegin[2]]) 
   
    # fault disp
    fault_disp = dispMat[:, allElemBegin[0]:allElemBegin[1]].dot(estSlip[allElemBegin[0]:allElemBegin[1]])

    with open(outputDir + 'numpy/cmiDisp.npy', 'wb') as f:
        np.save(f, cmi_disp)
    with open(outputDir + 'numpy/faultDisp.npy', 'wb') as f:
        np.save(f, fault_disp)

    # for each station, calculate the percentage of cmi and the percentage of fault contributions

    # use element wise division
    cmiPercentage = (cmi_disp / predDisp) * 100
    faultPercentage = (fault_disp / predDisp) * 100

    # for each given station, find the percent contributed by cmi and the percent from the fault
    # then, average each mechanisms contribution across all the stations
    totalCmiPercent = np.sum(cmiPercentage) / cmiPercentage.size
    totalFaultPercent = np.sum(faultPercentage) / faultPercentage.size

    weightedFaultPercent = np.average(faultPercentage, weights=np.abs(predDisp))
    weightedCmiPercent = np.average(cmiPercentage, weights=np.abs(predDisp))

    return totalFaultPercent, totalCmiPercent, weightedFaultPercent, weightedCmiPercent



##### ----- currently unused ----- #####


# method to output plots similar in style to Diao et al.
# observed displacement, cal_afterslip, cal_cmi, residual (horiz)
def plotLikeDiao(gps, predDisp, vectorLength, vecScale, dispMat, estSlip, allElemBegin, numDays, config) :
    saveFigures = config["results"]["saveFigures"]
    ratioFig = config["results"]["ratioFig"]
    outputDir = config["results"]["outputDir"]

    residuals = np.empty((len(gps.lon), 3))
    actual = np.hstack((np.array(gps.east_vel).reshape(-1,1), np.array(gps.north_vel).reshape(-1,1), np.array(gps.up_vel).reshape(-1,1)))

    # predicted / calculated values
    calc_east = predDisp[0::3].flatten()
    calc_north = predDisp[1::3].flatten()
    calc_up = predDisp[2::3].flatten()

    residuals[:,0] = actual[:,0] - calc_east
    residuals[:,1] = actual[:,1] - calc_north
    residuals[:,2] = actual[:,2] - calc_up

    # square the components
    east_disp = np.square(predDisp[0::3]).reshape(1,-1)
    north_disp = np.square(predDisp[1::3]).reshape(1, -1)

    # add together north and east disp, sum down column, take square root for magnitude of horiz disp
    totalDisp = np.sqrt(np.sum(np.vstack((east_disp, north_disp)), axis=0))

    # calc disp from cmi, beginning from the cmi elements to the end of the cmi elements
    cmi_disp = dispMat[:, allElemBegin[1]:allElemBegin[2]].dot(estSlip[allElemBegin[1]:allElemBegin[2]]) 
    cmi_e = np.square(cmi_disp[0::3]).reshape(1,-1) # square components
    cmi_n = np.square(cmi_disp[1::3]).reshape(1,-1)
    totalCmiDisp = np.sqrt(np.sum(np.vstack((cmi_e, cmi_n)), axis=0))

    # fault disp
    fault_disp = dispMat[:, allElemBegin[0]:allElemBegin[1]].dot(estSlip[allElemBegin[0]:allElemBegin[1]])

    plotRatio(gps, totalCmiDisp, totalDisp, saveFigures, ratioFig)

    fig, ax = plt.subplots(nrows=1, ncols=4, figsize=(22, 6), sharex=True, sharey=True)

    arrowWidth = 0.004

    # OBSERVED DISPLACEMENTS
    Q= ax[0].quiver(gps.lon, gps.lat, gps.east_vel/numDays, gps.north_vel/numDays, scale_units="xy", scale=vecScale, color='k', label='observed', width=arrowWidth)
    ax[0].quiverkey(Q, X = 0.3, Y=0.8, U=vectorLength, label='2 mm d-1',labelpos='N', color='k')
    ax[0].plot(coast.lon, coast.lat, color="k", linewidth=0.5) # coastline
    ax[0].set_title("Observed displacements (cm)")
    ax[0].set_ylim([35, 41])
    ax[0].set_xlim([135, 143])
    ax[0].set_aspect('equal', adjustable='box')

    # AFTERSLIP CONTRIBUTION # 

    Q2 = ax[1].quiver(gps.lon, gps.lat, fault_disp[0::3]/numDays, fault_disp[1::3]/numDays, scale_units="xy", scale=vecScale, color='g', label="displacements from afterslip", width=arrowWidth)
    ax[1].quiverkey(Q2, X=0.3, Y=0.8, U=vectorLength, label="2 mm d-1", labelpos='N', color='g')
    ax[1].plot(coast.lon, coast.lat, color="k", linewidth=0.5) # coastline
    ax[1].set_ylim([35, 41])
    ax[1].set_xlim([135, 143])
    ax[1].set_aspect('equal', adjustable='box')
    ax[1].set_title("Calculated Afterslip")

    # CMI CONTRIBUTION # 

    Q3 = ax[2].quiver(gps.lon, gps.lat, cmi_disp[0::3]/numDays, cmi_disp[1::3]/numDays, scale_units="xy", scale=vecScale, color='b', label="displacements from cmi slip", width=arrowWidth)
    ax[2].quiverkey(Q3, X=0.3, Y=0.8, U=vectorLength, label="2 mm d-1", labelpos='N', color='b')
    ax[2].plot(coast.lon, coast.lat, color="k", linewidth=0.5) # coastline
    ax[2].set_ylim([35, 41])
    ax[2].set_xlim([135, 143])
    ax[2].set_aspect('equal', adjustable='box')
    ax[2].set_title("Calculated CMI slip")

    
    # RESIDUALS #

    Q4 = ax[3].quiver(gps.lon, gps.lat, residuals[:,0]/numDays, residuals[:,1]/numDays, scale_units="xy", scale=vecScale, color='r', label="residuals", width=arrowWidth)
    ax[3].quiverkey(Q4, X=0.3, Y=0.8, U=vectorLength, label="2 mm d-1", labelpos='N', color='r')
    ax[3].plot(coast.lon, coast.lat, color="k", linewidth=0.5) # coastline
    ax[3].set_ylim([35, 41])
    ax[3].set_xlim([135, 143])
    ax[3].set_aspect('equal', adjustable='box')
    ax[3].set_title("Residual Displacements (horizontal)")
    plt.xticks([136, 137, 138, 139, 140, 141, 142])  # Set label locations.

    fig.subplots_adjust(wspace=0, hspace=0)

    plt.savefig(outputDir + "diaoFormattedDisplacements.pdf")
    plt.close('all')

    return


def observedCalculated(predDisp, gps, vecScale, estSlip, fault, config):
    outputDir = config["results"]["outputDir"]
    xmin = 130
    xmax = 144
    ymin = 32
    ymax = 43

    # dip slip on the CMI represents east-west motion, where east is negative and west is positive
    slip_type = 1 # 0 = strike slip 1 = dip slip, only for CMI, 2 (manual) is vertical/tensile
    end_idx = 2* len(fault["lon1"]) #end of fault elem beginning of cmi elem
    slip_vals = [estSlip[slip_type:end_idx:2]/100, estSlip[slip_type+end_idx::3]/100] # dip slip values for fault and CMI, converted from cm to m

    maxMag = np.max(np.abs(slip_vals[0]))

    fig, ax = plt.subplots( figsize=(8, 6))
    rso = ax.tripcolor(fault["points"][:,0],
                        fault["points"][:,1], 
                        fault["verts"],
                        facecolors=(slip_vals[0]).flatten(), 
                        vmin=-maxMag, vmax=maxMag)
    cbar1 = fig.colorbar(rso, ax=ax, orientation='vertical')
    cbar1.set_label("Dip Slip (m)")
    ax.plot(coast.lon, coast.lat, color="k", linewidth=0.5)
    Q = ax.quiver(gps.lon, gps.lat, gps.east_vel, gps.north_vel, scale=vecScale, color='k', label="observed")
    ax.quiverkey(Q, X=0.3, Y=0.8, U=50, label="50 cm", labelpos='N', color='k')
    ax.quiver(gps.lon, gps.lat, predDisp[0::3], predDisp[1::3], scale=vecScale, color='r', label="predicted")
    ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax), aspect='equal')
    ax.set_ylabel("Latitude")
    ax.set_xlabel("Longitude")
    plt.legend()
    plt.savefig(outputDir + "calcObs2.pdf")
    plt.close("all")
    return

def plotPercent(gps, predDisp, cmiDisp, outputDir):
    percent = 100*(cmiDisp/predDisp) #element wise division
    percent_up = percent[2::3]
    percent_east = percent[0::3]
    percent_north = percent[1::3]

    avgPercent = (percent_up + percent_north + percent_east)/3

    maxVal = 100
    minVal = -100

    # create masks for where the values are between -50-50%, these dots = more transparent, outliers = opaque plotting
    # the mask returns true for percentages either greater than 50 or less than -50
    avgMask = np.abs(avgPercent) > 50
    upMask = np.abs(percent_up) > 50
    eastMask = np.abs(percent_east) > 50
    northMask = np.abs(percent_up) > 50

    # convert to numpy arrays to use mask
    gpsLon = np.array(gps.lon).reshape(-1,1)
    gpsLat = np.array(gps.lat).reshape(-1,1)

    plt.close('all')
    fig, ax = plt.subplots(1,4, figsize=(20,5))
    dots = ax[0].scatter(gpsLon[avgMask], gpsLat[avgMask], c=avgPercent[avgMask], vmin=minVal, vmax=maxVal, cmap="PiYG") # color by ratio of cmi_disp / total disp
    ax[0].scatter(gpsLon[~avgMask], gpsLat[~avgMask], c=avgPercent[~avgMask], vmin=minVal, vmax=maxVal, cmap="PiYG", alpha=0.3) # flip the boolean mask to get the complement
    ax[0].plot(coast.lon, coast.lat, color="k", linewidth=0.5)
    ax[0].set_xlim(130, 145)
    ax[0].set_ylim(33, 45)
    ax[0].set_title("average percent")
    # Add a customized colorbar
    fig.colorbar(dots, ax=ax[0], orientation='vertical', shrink=0.7, extend='both')
    
    dots1 = ax[1].scatter(gpsLon[upMask], gpsLat[upMask], c=percent_up[upMask], vmin=minVal, vmax=maxVal, cmap="PiYG") # color by ratio of cmi_disp / total disp
    ax[1].scatter(gpsLon[~upMask], gpsLat[~upMask], c=percent_up[~upMask], vmin=minVal, vmax=maxVal, cmap="PiYG", alpha=0.3)
    ax[1].plot(coast.lon, coast.lat, color="k", linewidth=0.5)
    ax[1].set_xlim(130, 145)
    ax[1].set_ylim(33, 45)
    ax[1].set_title("vertical percent")
    fig.colorbar(dots1, ax=ax[1], orientation='vertical', shrink=0.7, extend='both')
    
    dots2 = ax[2].scatter(gpsLon[eastMask], gpsLat[eastMask], c=percent_east[eastMask], vmin=minVal, vmax=maxVal, cmap="PiYG") # color by ratio of cmi_disp / total disp
    ax[2].scatter(gpsLon[~eastMask], gpsLat[~eastMask], c=percent_east[~eastMask], vmin=minVal, vmax=maxVal, cmap="PiYG", alpha=0.3) 
    ax[2].plot(coast.lon, coast.lat, color="k", linewidth=0.5)
    ax[2].set_xlim(130, 145)
    ax[2].set_ylim(33, 45)
    ax[2].set_title("east percent")
    fig.colorbar(dots2, ax=ax[2], orientation='vertical', shrink=0.7, extend='both')
    
    dots3 = ax[3].scatter(gpsLon[northMask], gpsLat[northMask], c=percent_north[northMask], vmin=minVal, vmax=maxVal, cmap="PiYG") # color by ratio of cmi_disp / total disp
    ax[3].scatter(gpsLon[~northMask], gpsLat[~northMask], c=percent_north[~northMask], vmin=minVal, vmax=maxVal, cmap="PiYG", alpha=0.3)
    ax[3].plot(coast.lon, coast.lat, color="k", linewidth=0.5)
    ax[3].set_xlim(130, 145)
    ax[3].set_ylim(33, 45)
    ax[3].set_title("north percent")
    fig.colorbar(dots3, ax=ax[3], orientation='vertical', shrink=0.7, 
                        label='percent disp from cmi', extend='both')

    plt.savefig(outputDir + "percentCmi.pdf")
    
    plt.close('all')
    return



# freed match for horiz
    # fig, ax = plt.subplots(1, 2, figsize=(10,5))
    # ax[0].plot(coast.lon, coast.lat, color="k", linewidth=0.5)
    # ax[1].plot(coast.lon, coast.lat, color="k", linewidth=0.5)

    # # plot the component of horizontal displacement from afterslip against total
    
    # # plot the component of horizontal displacement from observed
    # Q2 = ax[0].quiver(gps.lon, gps.lat, gps.east_vel, gps.north_vel, 
    #                   scale=miniVecScale*2, color='r', label="observed", 
    #                   headlength=headLength, headaxislength=headAxisLength, headwidth=headWidth, width=arrowWidth, 
    #                   lw=0.2, ec='k',)
    # ax[0].quiver(gps.lon, gps.lat, fault_disp[0::3], fault_disp[1::3], 
    #              scale=miniVecScale*2, color='cyan', label="fault", 
    #              headlength=headLength, headaxislength=headAxisLength, headwidth=headWidth-1, width=arrowWidth-0.0005, 
    #              lw=0.2, ec='k',)
    
    # # quiver based off of cmi disp, but same units as afterslip so it should be alright
    # ax[0].quiverkey(Q2, X=0.3, Y=0.8, U=100, label="100 cm", labelpos='N', color='r')
    # ax[0].set_ylim([34, 42])
    # ax[0].set_xlim([136, 144])
    # ax[0].set_title("Horizontal afterslip contributions")

    # Q3 = ax[1].quiver(gps.lon, gps.lat, gps.east_vel, gps.north_vel, 
    #                   scale=miniVecScale*2, color='r', label="observed", 
    #                   headlength=headLength, headaxislength=headAxisLength, headwidth=headWidth, width=arrowWidth, 
    #                   lw=0.2, ec='k',)
    # ax[1].quiver(gps.lon, gps.lat, cmi_disp[0::3], cmi_disp[1::3], 
    #              scale=miniVecScale*2, color='yellow', label="cmi", 
    #              headlength=headLength, headaxislength=headAxisLength, headwidth=headWidth-1, width=arrowWidth-0.0005, 
    #              lw=0.2, ec='k',)
    
    # # quiver based off of cmi disp, but same units as afterslip so it should be alright
    # ax[1].quiverkey(Q3, X=0.3, Y=0.8, U=100, label="100 cm", labelpos='N', color='r')
    # ax[1].set_ylim([34, 42])
    # ax[1].set_xlim([136, 144])
    # ax[1].set_title("Horizontal cmi contributions")

    # plt.savefig(outputDir + 'images/dispByComponent_FreedHoriz.pdf')
    # plt.close('all')