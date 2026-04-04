import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import pandas as pd
import json
import yaml

# global
coast = pd.read_csv("coastline.csv")
lon_corr = 1

def slipDist(estSlip, gps, fault, cmi, vecScale, slipDist=False, saveFigures=False) :
    '''Plots the slip distribution for both the cmi alone with gps vectors
    and the fault connected to the cmi without gps vectors'''
    
    # plot for visualizing
    xmin = 122
    xmax = 148
    ymin = 32
    ymax = 48
    plt.close('all')

    # dip slip on the CMI represents east-west motion, where east is negative and west is positive
    slip_type = 1 # 0 = strike slip 1 = dip slip, only for CMI, 2 (manual) is vertical/tensile
    end_idx = 2* len(fault["lon1"]) #end of fault elem beginning of cmi elem
    slip_vals = [estSlip[slip_type:end_idx:2]/100, estSlip[slip_type+end_idx::3]/100] # dip slip values for fault and CMI, converted from cm to m

    max_mag_f = np.max(np.abs(slip_vals[0]))
    max_mag_h = np.max(np.abs(slip_vals[1]))
    if max_mag_f > max_mag_h:
        max_mag = max_mag_f
    else:
        max_mag = max_mag_h

    both = {}
    both["points"] = np.vstack((fault["points"], cmi["points"]))
    shift_val = len(fault["points"][:,0])
    both["verts"] = np.vstack((fault["verts"], cmi["verts"]+shift_val))

    fig, ax = plt.subplots(1, 2, figsize=(10,6))
    rso = ax[0].tripcolor(both["points"][:,0],
                        both["points"][:,1], 
                        both["verts"],
                        facecolors=(np.vstack(((slip_vals[0], slip_vals[1])))).flatten(), 
                        vmin=-max_mag, vmax=max_mag)

    cbar1 = fig.colorbar(rso, ax=ax[0], orientation='horizontal')
    ax[0].plot(coast.lon+360*(1-lon_corr), coast.lat, color="k", linewidth=0.5)
    cbar1.set_label("Slip (m)")
    ax[0].set(xlim=(xmin-2, xmax), ylim=(ymin, ymax), aspect='equal')
    ax[0].title.set_text("Fault Dip Slip") #graph 1
    ax[0].set_ylabel("Latitude")
    ax[0].set_xlabel("Longitude")

    rso = ax[1].tripcolor(cmi["points"][:,0], cmi["points"][:,1], cmi["verts"], facecolors=(slip_vals[1]).flatten(), vmin=-max_mag_h, vmax=max_mag_h)
    cbar1 = fig.colorbar(rso, ax=ax[1], orientation='horizontal')
    cbar1.set_label("Slip (m)")
    ax[1].quiver(gps.lon, gps.lat, gps.east_vel, gps.north_vel, scale=vecScale, color='k', label='observed')
    ax[1].set(xlim=(xmin-2, xmax), ylim=(ymin, ymax), aspect='equal')
    ax[1].title.set_text("CMI Dip Slip") #graph 1
    ax[1].set_ylabel("Latitude")
    ax[1].set_xlabel("Longitude")

    if saveFigures and slipDist:
        plt.savefig('slip_dist.pdf') # save the figure
    
    plt.close('all')

    return

def observedCalculated(predDisp, gps, vecScale, estSlip, fault):
    xmin = 130
    xmax = 144
    ymin = 32
    ymax = 43

    coast = pd.read_csv("coastline.csv")
    lon_corr = 1

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
    ax.plot(coast.lon+360*(1-lon_corr), coast.lat, color="k", linewidth=0.5)
    Q = ax.quiver(gps.lon, gps.lat, gps.east_vel, gps.north_vel, scale=vecScale, color='k', label="observed")
    ax.quiverkey(Q, X=0.3, Y=0.8, U=50, label="50 cm", labelpos='N', color='k')
    ax.quiver(gps.lon, gps.lat, predDisp[0::3], predDisp[1::3], scale=vecScale, color='r', label="predicted")
    ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax), aspect='equal')
    ax.set_ylabel("Latitude")
    ax.set_xlabel("Longitude")
    plt.legend()
    plt.savefig("calcObs2.pdf")
    plt.close("all")

    return


def afterslip(estSlip, fault):
    '''plots the afterslip distribution along the subduction zone,
    side by side with a plot showing the location of the peak afterslip'''

    # plot for visualizing
    xmin = 140
    xmax = 148
    ymin = 32
    ymax = 48
    plt.close('all')

    coast = pd.read_csv("coastline.csv")
    lon_corr = 1
    # dip slip on the CMI represents east-west motion, where east is negative and west is positive
    slip_type = 1 # 0 = strike slip 1 = dip slip, only for CMI, 2 (manual) is vertical/tensile
    end_idx = 2*len(fault["lon1"]) #end of fault elem beginning of cmi elem
    
    dip_slip_vals = np.array(estSlip[slip_type:end_idx:2]/100) # dip slip values for fault only
    strike_slip_vals = np.array(estSlip[0:end_idx:2]/100)

    # slip is a combination of dip slip and strike slip
    # ** cmi has tensile slip too but for now plotting is in cartesian plane
    # slip = (dip_slip^2 + strike_slip^2) ^ (1/2)
    slip = np.sqrt(np.square(dip_slip_vals) + np.square(strike_slip_vals))

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12, 8))
    rso = ax[0].tripcolor(fault["points"][:,0],
                        fault["points"][:,1], 
                        fault["verts"],
                        facecolors=(slip).flatten(), 
                        vmin=-6, vmax=6)
    cbar1 = fig.colorbar(rso, ax=ax[0], orientation='vertical')
    ax[0].plot(coast.lon+360*(1-lon_corr), coast.lat, color="k", linewidth=0.5)
    cbar1.set_label("Slip (m)")
    ax[0].set(xlim=(xmin-2, xmax), ylim=(ymin, ymax), aspect='equal')
    ax[0].title.set_text("Fault Slip") #graph 1
    ax[0].set_ylabel("Latitude")
    ax[0].set_xlabel("Longitude")

    custom_colors = ['lightsteelblue','lightsteelblue', 'royalblue', 'blue']
    cmap = colors.ListedColormap(custom_colors)
    maxSlip = np.max(slip)
    bounds = [0, maxSlip-1, maxSlip-1, maxSlip-0.1, maxSlip+0.1] 
    norm = colors.BoundaryNorm(bounds, cmap.N)
    rso2 = ax[1].tripcolor(fault["points"][:,0], fault["points"][:,1], fault["verts"], facecolors=(slip).flatten(), cmap=cmap,norm=norm)
    cbar2 = fig.colorbar(rso2, ax=ax[1], boundaries=bounds,orientation='vertical')
    cbar2.set_label('slip magnitude (m)')
    ax[1].plot(coast.lon+360*(1-lon_corr), coast.lat, color="k", linewidth=0.5)
    ax[1].set(xlim=(xmin-2, xmax), ylim=(ymin, ymax), aspect='equal')
    ax[1].title.set_text("Max Slip Location") #graph 2

    plt.savefig("afterslip.pdf")
    plt.close('all')

    return


def displacements(dispMat, allElemBegin, estSlip, predDisp, gps, vecScale, saveFigures=False, allDisp =False, dispSep=False, ratioFig=False) :
    '''plot displacements, including observed, predicted, and displacements separated
    by the component (i.e., displacement from CMI and displacement from fault)'''
    # calculate displacements by which components
    ### VECTOR STYLING ###
    miniVecScale = 500 # mini vec scale used for shorter vectors
    arrowWidth = 0.005
    arrowWidth2 = 0.003
    headLength = 2
    headAxisLength=2
    headWidth = 3

    coast = pd.read_csv("coastline.csv") 

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
    fault_e = np.square(fault_disp[0::3]).reshape(1,-1) # square components
    fault_n = np.square(fault_disp[1::3]).reshape(1,-1)
    totalFaultDisp = np.sqrt(np.sum(np.vstack((fault_e, fault_n)), axis=0))

    plotRatio(gps, totalCmiDisp, totalDisp, saveFigures, ratioFig)
    plotPercent(gps, predDisp, cmi_disp, saveFigures)

    fig, ax = plt.subplots(1, 2, figsize=(10,5))
    ax[0].plot(coast.lon, coast.lat, color="k", linewidth=0.5)
    ax[1].plot(coast.lon, coast.lat, color="k", linewidth=0.5)
    Q= ax[0].quiver(gps.lon, gps.lat, gps.east_vel, gps.north_vel, scale=vecScale, color='k', label='observed', headlength=headLength, headaxislength=headAxisLength, headwidth=headWidth, width=arrowWidth, 
                 lw=0.2, ec='k',)
    ax[0].quiverkey(Q, X = 0.3, Y=0.7, U=100, label='100 cm',labelpos='N', color='k')
    ax[0].set_title("Observed and Predicted Horizontal (cm)")
    ax[0].set_ylim([34, 42])
    ax[0].set_xlim([136, 144])
    ax[0].quiver(gps.lon, gps.lat, predDisp[0::3], predDisp[1::3], scale=vecScale, color='r', label='predicted', headlength=headLength, headaxislength=headAxisLength, headwidth=headWidth, width=arrowWidth, 
                 lw=0.2, ec='k',)
    

    Q0 = ax[1].quiver(gps.lon, gps.lat, np.zeros_like(gps.north_vel), gps.up_vel, scale=vecScale/10, color='k', label='observed', headlength=headLength, headaxislength=headAxisLength, headwidth=headWidth, width=arrowWidth, 
                 lw=0.2, ec='k',)
    ax[1].quiverkey(Q0, X = 0.3, Y=0.7, U=10, label='10 cm',labelpos='N', color='k')
    ax[1].set_title("Observed and Predicted Vertical (cm)")
    ax[1].set_ylim([34, 42])
    ax[1].set_xlim([136, 144])

    ax[1].quiver(gps.lon, gps.lat, np.zeros_like(predDisp[0::3]), predDisp[2::3], scale=vecScale/10, color='r', label='predicted', headlength=headLength, headaxislength=headAxisLength, headwidth=headWidth, width=arrowWidth, 
                 lw=0.2, ec='k',)
    # ax[0].legend()
    ax[1].legend()

    if saveFigures and allDisp:
        plt.savefig('totalDisp.pdf')

    miniVecScale = 500

    fig, ax = plt.subplots(1, 2, figsize=(10,5))
    # plotting of coastline (deleted during ax.clear() so it needs to be updated with each new plotting)
    ax[0].plot(coast.lon, coast.lat, color="k", linewidth=0.5)
    ax[1].plot(coast.lon, coast.lat, color="k", linewidth=0.5)

    # plot the component of horizontal displacement from afterslip
    ax[0].quiver(gps.lon, gps.lat, fault_disp[0::3], fault_disp[1::3], 
                 scale=miniVecScale, color='b', label="fault contribution", 
                 headlength=headLength, headaxislength=headAxisLength, headwidth=headWidth, width=arrowWidth, 
                 lw=0.2, ec='k',)
    
    # plot the component of horizontal displacement from VER / cmi slip
    Q2 = ax[0].quiver(gps.lon, gps.lat, cmi_disp[0::3], cmi_disp[1::3], 
                      scale=miniVecScale, color='y', label="cmi contribution", 
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
                 scale=miniVecScale/5, color='b', label="fault contribution", 
                 headlength=headLength, headaxislength=headAxisLength, headwidth=headWidth, width=arrowWidth, 
                 lw=0.2, ec='k',)
    
    # plot the component of vertical displacement from VER / cmi
    Q3 = ax[1].quiver(gps.lon, gps.lat, ogVec, cmi_disp[2::3], 
                      scale=miniVecScale/5, color='y', label="cmi contribution", 
                      headlength=headLength, headaxislength=headAxisLength, headwidth=headWidth, width=arrowWidth2, 
                      lw=0.2, ec='k',)
    ax[1].quiverkey(Q3, X=0.3, Y=0.8, U=10, label="10 cm", labelpos='N', color='b')
    ax[1].set_ylim([34, 42])
    ax[1].set_xlim([136, 144])
    ax[1].set_title("Vertical contribution")
    ax[1].legend()

    if saveFigures and dispSep:
        plt.savefig('dispByComponent.pdf')
    
    plt.close('all')
    return


# plot the ratio of displacement due to cmi vs total displacements
def plotRatio(gps, totalCmiDisp, totalDisp, saveFigures, ratioFig) :
    ratio = totalCmiDisp / totalDisp
    lon_corr = 1 # longitude correction hardcoded to 1

    coast = pd.read_csv("coastline.csv")

    plt.close('all')
    fig, ax = plt.subplots()
    dots = ax.scatter(gps.lon, gps.lat, c=ratio) # color by ratio of cmi_disp / total disp
    ax.plot(coast.lon+360*(1-lon_corr), coast.lat, color="k", linewidth=0.5)
    ax.set_xlim(130, 145)
    ax.set_ylim(33, 45)
    # Add a customized colorbar
    cbar = fig.colorbar(dots, ax=ax, orientation='vertical', shrink=0.7,
                        label='Ratio of cmi:total', extend='both')

    if saveFigures and ratioFig:
        plt.savefig("ratioCmiToTotal.pdf")
    
    plt.close('all')

    return

def plotPercent(gps, predDisp, cmiDisp, saveFigures):
    percent = 100*(cmiDisp/predDisp) #element wise division
    percent_up = percent[2::3]
    percent_east = percent[0::3]
    percent_north = percent[1::3]

    avgPercent = (percent_up + percent_north + percent_east)/3
    lon_corr = 1 # longitude correction hardcoded to 1

    coast = pd.read_csv("coastline.csv")

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
    ax[0].plot(coast.lon+360*(1-lon_corr), coast.lat, color="k", linewidth=0.5)
    ax[0].set_xlim(130, 145)
    ax[0].set_ylim(33, 45)
    ax[0].set_title("average percent")
    # Add a customized colorbar
    fig.colorbar(dots, ax=ax[0], orientation='vertical', shrink=0.7, extend='both')
    
    dots1 = ax[1].scatter(gpsLon[upMask], gpsLat[upMask], c=percent_up[upMask], vmin=minVal, vmax=maxVal, cmap="PiYG") # color by ratio of cmi_disp / total disp
    ax[1].scatter(gpsLon[~upMask], gpsLat[~upMask], c=percent_up[~upMask], vmin=minVal, vmax=maxVal, cmap="PiYG", alpha=0.3)
    ax[1].plot(coast.lon+360*(1-lon_corr), coast.lat, color="k", linewidth=0.5)
    ax[1].set_xlim(130, 145)
    ax[1].set_ylim(33, 45)
    ax[1].set_title("vertical percent")
    fig.colorbar(dots1, ax=ax[1], orientation='vertical', shrink=0.7, extend='both')
    
    dots2 = ax[2].scatter(gpsLon[eastMask], gpsLat[eastMask], c=percent_east[eastMask], vmin=minVal, vmax=maxVal, cmap="PiYG") # color by ratio of cmi_disp / total disp
    ax[2].scatter(gpsLon[~eastMask], gpsLat[~eastMask], c=percent_east[~eastMask], vmin=minVal, vmax=maxVal, cmap="PiYG", alpha=0.3) 
    ax[2].plot(coast.lon+360*(1-lon_corr), coast.lat, color="k", linewidth=0.5)
    ax[2].set_xlim(130, 145)
    ax[2].set_ylim(33, 45)
    ax[2].set_title("east percent")
    fig.colorbar(dots2, ax=ax[2], orientation='vertical', shrink=0.7, extend='both')
    
    dots3 = ax[3].scatter(gpsLon[northMask], gpsLat[northMask], c=percent_north[northMask], vmin=minVal, vmax=maxVal, cmap="PiYG") # color by ratio of cmi_disp / total disp
    ax[3].scatter(gpsLon[~northMask], gpsLat[~northMask], c=percent_north[~northMask], vmin=minVal, vmax=maxVal, cmap="PiYG", alpha=0.3)
    ax[3].plot(coast.lon+360*(1-lon_corr), coast.lat, color="k", linewidth=0.5)
    ax[3].set_xlim(130, 145)
    ax[3].set_ylim(33, 45)
    ax[3].set_title("north percent")
    fig.colorbar(dots3, ax=ax[3], orientation='vertical', shrink=0.7, 
                        label='percent disp from cmi', extend='both')

    if saveFigures:
        plt.savefig("percentCmi.pdf")
    
    plt.close('all')

    return


def residualPlot(gps, predDisp, vecScale, saveFigures=False, residFig=False) :
    # calculate gps displacement residuals
    # residual = actual - predicted
     ### VECTOR STYLING ###
    arrowWidth = 0.005
    headLength = 2
    headAxisLength=2
    headWidth = 3

    coast = pd.read_csv("coastline.csv")

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

    Q1 = ax[1].quiver(gps.lon, gps.lat, 0, residuals[:,2], scale=vecScale/5, color='r', label="residuals", 
                      headlength=headLength, headaxislength=headAxisLength, headwidth=headWidth, width=arrowWidth, 
                      lw=0.1, ec='k',)
    ax[1].quiverkey(Q1, X=0.3, Y=0.8, U=10, label='10 cm', labelpos='N', color='r')
    ax[1].plot(coast.lon, coast.lat, color="k", linewidth=0.5)
    ax[1].set_ylim([34, 42])
    ax[1].set_xlim([136, 144])
    plt.title("residual displacements (vertical)")

    if saveFigures and residFig:
        plt.savefig("residuals.pdf")
    
    plt.close('all')

    return


# method to output plots similar in style to Diao et al.
# observed displacement, cal_afterslip, cal_cmi, residual (horiz)
def plotLikeDiao(gps, predDisp, vectorLength, vecScale, dispMat, estSlip, allElemBegin, numDays, saveFigures=False, ratioFig=False) :

    residuals = np.empty((len(gps.lon), 3))
    actual = np.hstack((np.array(gps.east_vel).reshape(-1,1), np.array(gps.north_vel).reshape(-1,1), np.array(gps.up_vel).reshape(-1,1)))

    # predicted / calculated values
    calc_east = predDisp[0::3].flatten()
    calc_north = predDisp[1::3].flatten()
    calc_up = predDisp[2::3].flatten()

    residuals[:,0] = actual[:,0] - calc_east
    residuals[:,1] = actual[:,1] - calc_north
    residuals[:,2] = actual[:,2] - calc_up

    coast = pd.read_csv("coastline.csv")
    lon_corr = 1

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
    ax[0].plot(coast.lon+360*(1-lon_corr), coast.lat, color="k", linewidth=0.5) # coastline
    ax[0].set_title("Observed displacements (cm)")
    ax[0].set_ylim([35, 41])
    ax[0].set_xlim([135, 143])
    ax[0].set_aspect('equal', adjustable='box')

    # AFTERSLIP CONTRIBUTION # 

    Q2 = ax[1].quiver(gps.lon, gps.lat, fault_disp[0::3]/numDays, fault_disp[1::3]/numDays, scale_units="xy", scale=vecScale, color='g', label="displacements from afterslip", width=arrowWidth)
    ax[1].quiverkey(Q2, X=0.3, Y=0.8, U=vectorLength, label="2 mm d-1", labelpos='N', color='g')
    ax[1].plot(coast.lon+360*(1-lon_corr), coast.lat, color="k", linewidth=0.5) # coastline
    ax[1].set_ylim([35, 41])
    ax[1].set_xlim([135, 143])
    ax[1].set_aspect('equal', adjustable='box')
    ax[1].set_title("Calculated Afterslip")

    # CMI CONTRIBUTION # 

    Q3 = ax[2].quiver(gps.lon, gps.lat, cmi_disp[0::3]/numDays, cmi_disp[1::3]/numDays, scale_units="xy", scale=vecScale, color='b', label="displacements from cmi slip", width=arrowWidth)
    ax[2].quiverkey(Q3, X=0.3, Y=0.8, U=vectorLength, label="2 mm d-1", labelpos='N', color='b')
    ax[2].plot(coast.lon+360*(1-lon_corr), coast.lat, color="k", linewidth=0.5) # coastline
    ax[2].set_ylim([35, 41])
    ax[2].set_xlim([135, 143])
    ax[2].set_aspect('equal', adjustable='box')
    ax[2].set_title("Calculated CMI slip")

    
    # RESIDUALS #

    Q4 = ax[3].quiver(gps.lon, gps.lat, residuals[:,0]/numDays, residuals[:,1]/numDays, scale_units="xy", scale=vecScale, color='r', label="residuals", width=arrowWidth)
    ax[3].quiverkey(Q4, X=0.3, Y=0.8, U=vectorLength, label="2 mm d-1", labelpos='N', color='r')
    ax[3].plot(coast.lon+360*(1-lon_corr), coast.lat, color="k", linewidth=0.5) # coastline
    ax[3].set_ylim([35, 41])
    ax[3].set_xlim([135, 143])
    ax[3].set_aspect('equal', adjustable='box')
    ax[3].set_title("Residual Displacements (horizontal)")
    plt.xticks([136, 137, 138, 139, 140, 141, 142])  # Set label locations.

    fig.subplots_adjust(wspace=0, hspace=0)

    plt.savefig("diaoFormattedDisplacements.pdf")
    plt.close('all')

    return

def maxSlipMag(estSlip, allElemBegin):
    slip_vals = [estSlip[0:allElemBegin[1]]/100, estSlip[allElemBegin[1]::]/100] # slip values for fault and CMI, converted from cm to m

    maxFaultMag = np.max(np.abs(slip_vals[0]))
    maxCmiMag = np.max(np.abs(slip_vals[1]))

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

def calcMoment(estSlip, allElemBegin, fault, cmi):
    '''calculate the moment for slip on fault vs slip on cmi
        moment = rigidity x area x slip'''

    rigidity = 3e10 #N / m^2 (is this the same though for the CMI at depth?)

    # allElemBegin[1] corresponds to the end of the fault elements and beginning of the cmi elements
    # calc slip mag
    fault_d_slip = np.square(estSlip[1:allElemBegin[1]:2]) # one for dip slip
    fault_s_slip = np.square(estSlip[0:allElemBegin[1]:2]) # zero for strike slip
    faultSlip = np.sqrt(np.sum(np.vstack((fault_d_slip, fault_s_slip)), axis=0))/100 # convert cm -> m

    faultMoment = rigidity * (fault["area"]) * faultSlip
    faultTotalMoment = np.sum(faultMoment) # in Nm

    cmi_d_slip = np.square(estSlip[1+allElemBegin[1]::3])
    cmi_s_slip = np.square(estSlip[0+allElemBegin[1]::3])
    cmiSlip = np.sqrt(np.sum(np.vstack((cmi_d_slip, cmi_s_slip)), axis=0))/100

    cmiMoment = rigidity * (cmi["area"]) * cmiSlip
    cmiTotalMoment = np.sum(cmiMoment)

    return faultTotalMoment, cmiTotalMoment


def numericalData(estSlip, predDisp, dispMat, gps, allElemBegin, fault, cmi, saveData):
    maxFaultMag, maxCmiMag = maxSlipMag(estSlip, allElemBegin)
    rmse = calcRMSE(predDisp, gps)
    faultMoment, cmiMoment = calcMoment(estSlip, allElemBegin, fault, cmi)
    avgFaultRake, avgCmiRake = rakeCalc(estSlip, allElemBegin)
    percentFault, percentCmi, weightedFault, weightedCmi = contributionsPercent(predDisp=predDisp, dispMat=dispMat, estSlip=estSlip, allElemBegin=allElemBegin)

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

    if (saveData) :
        with open("numericalResults.txt", "w") as file:
            file.write(json.dumps(results, indent=4, sort_keys=True))

    return


# save current config settings for later reference
def saveConfig(config):
    with open('configSettings.yaml', 'w') as file:
        yaml.dump(config, file, default_flow_style=False)
    return


def rakeCalc(estSlip, allElemBegin) :
    '''calculates the average rake value for the cmi and fault
    using the slip magnitude to weight the rake value of the given element'''

    # Slip magnitude weighted average of rake # (lots of slip = we trust it more)
    faultDipSlip = estSlip[1:allElemBegin[1]:2] # one for dip slip
    faultStrikeSlip = estSlip[0:allElemBegin[1]:2] # zero for strike slip
    faultSlipMag = np.sqrt(np.square(faultDipSlip) + np.square(faultStrikeSlip))

    cmiDipSlip = estSlip[1+allElemBegin[1]::3]
    cmiStrikeSlip = estSlip[0+allElemBegin[1]::3]
    cmiSlipMag = np.sqrt(np.sum(np.vstack((np.square(cmiDipSlip).reshape(1, -1), np.square(cmiStrikeSlip).reshape(1,-1))), axis=0))

    cmiRake = np.rad2deg(np.arctan2(cmiDipSlip.flatten(), cmiStrikeSlip.flatten())) # (left lateral = 0 rake)
    faultRake = np.rad2deg(np.arctan2(faultDipSlip.flatten(), faultStrikeSlip.flatten()))

    weightedCmiRake = np.average(cmiRake, weights=cmiSlipMag.flatten())
    weightedFaultRake = np.average(faultRake, weights=faultSlipMag.flatten())

    return weightedFaultRake, weightedCmiRake


def contributionsPercent(predDisp, dispMat, estSlip, allElemBegin) :
    ''' function for returning the percentage of total displacement that the cmi is responsible for
    compared to total percentage of displacements that the afterslip is contributing
    the percentages are calculated component (e/n/u) wise because trying to calc percentages 
    using the total displacement vector magnitude compounds rounding differences from np.float64'''

    # calc disp from cmi, beginning from the cmi elements to the end of the cmi elements
    cmi_disp = dispMat[:, allElemBegin[1]:allElemBegin[2]].dot(estSlip[allElemBegin[1]:allElemBegin[2]]) 
   
    # fault disp
    fault_disp = dispMat[:, allElemBegin[0]:allElemBegin[1]].dot(estSlip[allElemBegin[0]:allElemBegin[1]])

    with open('cmiDisp.npy', 'wb') as f:
        np.save(f, cmi_disp)
    with open('faultDisp.npy', 'wb') as f:
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
