import numpy as np
import pandas as pd


def assembleWeights(numMeshes, assembledMat, dispArray, smoothingWeights, allElemBegin, allElemEnd) :
    '''Assembled the weights matrix composed of
    smoothing weight for fault mesh
    smoothing weight for cmi mesh
    constraint matrix weights (fault edge and cmi edge)
    gps data weights '''

    # Start with unit uncertainties 
    # this puts the smoothing weight for the fault mesh and cmi mesh, and leaves gps stations as 1
    weights = np.ones((np.shape(assembledMat)[0], 1)) # matches the size of the assembled meshes matrix

    for mesh_idx in range(numMeshes+2): # numMeshes is faultMesh and cmiMesh, we add in the fault constraint matrix and cmi constraint matrix
        weights[np.size(dispArray) + allElemBegin[mesh_idx]:np.size(dispArray)+allElemEnd[mesh_idx]] = smoothingWeights[mesh_idx]

    return weights

def runInversion(assembledMat, dispMat, weights, dataVector) :
    # Calculate model covariance
    cov = np.linalg.inv(assembledMat.T * weights.T @ assembledMat) 

    # Estimate slip using pre-calculated covariance
    estSlip = cov @ assembledMat.T * weights.T @ dataVector 
    # Predict displacement at stations
    # note that .dot() will behave like @ or matmul because both dispMat and estSlip are 2D arrays
    predDisp = dispMat.dot(estSlip) 
    # run to check sign convention of dip slip (neg = east pos = west on CMI)
    # pred_disp = disp_mat[:, 1+all_elem_beg[1]::3].dot(est_slip[1+all_elem_beg[1]::3])

    return estSlip, predDisp