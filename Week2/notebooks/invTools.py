import scipy.linalg as LA
#import numpy.linalg as LA2

def invertLS(K,y,Seps,doError=False,Sa=None):
    # assume diagonal Seps
    dim = K.shape
    class fits:
        pass
    # rewrite Rodgers inversion as simple least squares (KKx=yy)
    
    if Sa is not None:
        # Compute Sa^-1/2 (should be precomputed next time)
        invSQRTSa = LA.inv(LA.cholesky(Sa).T)
        
        KK = np.vstack((K.copy(),invSQRTSa))
        yy = np.concatenate((y.copy(),np.zeros(dim[1],)))
    else:
        KK = K.copy()
        yy = y.copy()
    
    # Precondition K and y (normalize by noise):
    for i in range(dim[0]):
        KK[i,:] = KK[i,:]/np.sqrt(Seps[i])
        yy[i] = yy[i]/np.sqrt(Seps[i])
    
    # Solve least squares by QR decomposition with pivoting 
    # (turned off pivoting for now because it drove me mad when computing Shat and errors)
    Q,R = LA.qr(KK,mode='economic')
    pp  = np.dot(Q.T, yy)
    
    # Solve for max a posteriori solution vector 
    fits.xhat = np.dot(LA.inv(R),pp)
    #x = LA.lstsq(KK,yy)
    if doError:
        # Error analysis (can be later switched on/off)
        invShat = R.T.dot(R)
        # Compute posterior covariance matrix
        Shat = LA.inv(invShat)
        # Compute Averaging Kernel Matrix
        fits.A = Shat.dot(invShat-LA.inv(Sa))
        fits.DOF = np.trace(A)
        fits.Shat = Shat
        # Compute model
        fits.ymod = K.dot(fits.xhat)
        # Compute residuals
        res = fits.ymod-y
        # Compute reduced chi2 (using the true DOF for the # of measurements)
        fits.chi2r = np.sum((res/np.sqrt(Seps))**2)/(dim[0]-fits.DOF)
        fits.xnorm = LA.norm(fits.xhat)
      
    return fits 