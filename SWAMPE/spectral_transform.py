import numpy as np
import scipy.special as sp
import math


def PmnHmn(J,M,N,mus):
    """
    Calculates the values of associated Legendre polynomials and their 
     derivatives evaluated at Gaussian latitudes (mus) up to wavenumber M 
        
    :param J: number of latitudes
    :type J: int
    :param M: highest wavenumber for associated Legendre polynomials
    :type M: int
    :param N: highest degree of the Legendre functions for m=0
    :type N: int
    :param mus: Gaussian latitudes
    :type mus: array of float
    
    :return: Pmn array (associated legendre polynomials), Hmn array (derivatives of Pmn*(1-x^2)), both evaluated at the Gaussian latitudes mus
    :rtype: array of float
    """
    
    Pmn=np.zeros((J,M+1,N+1))
    Hmn=np.zeros((J,M+1,N+1))
    Pmntemp=np.zeros((J,M+1,N+1))
    Hmntemp=np.zeros((J,M+1,N+1))
    for j in range (0,J):
        Pmntemp[j],Hmntemp[j] = sp.lpmn(M,N,mus[j])
        Hmntemp[j,:,:] = (1-mus[j]**2)*Hmntemp[j,:,:]
        
    for m in range (0,M+1):
        for n in range (m,N+1):
            scaling_term = np.sqrt((((2*n)+1)*math.factorial(n-m))/(2*math.factorial(n+m)))
            Pmn[:,m,n] = scaling_term*Pmntemp[:,m,n]
            Hmn[:,m,n] = scaling_term*Hmntemp[:,m,n]

            if m % 2 == 1:
                if n>0:
                    Pmn[:,m,n]=-Pmn[:,m,n]
                    Hmn[:,m,n]=-Hmn[:,m,n]              
                
    return Pmn, Hmn


def fwd_leg(data,J,M,N,Pmn,w):
    """Calculates the forward legendre transform

        :param data: input to be transformed (usually output of fft)
        :type data: array of float or array of complex
        
        :param J: number of latitudes
        :type J: int
        
        :param M: highest wavenumber for associated Legendre polynomials
        :type M: int
        
        :param N: highest degree of the Legendre functions for m=0
        :type N: int
        
        :param Pmn: associated legendre functions evaluated at the Gaussian latitudes mus  up to wavenumber M
        :type Pmn: array of float
        
        :param w: Gauss Legendre weights
        :type w: array of float

        :return legcoeff: Legendre coefficients (if data was output from FFT, then legcoeff are the spectral coefficients)
        :rtype legcoeff: array of complex
        """
    legterm=np.zeros((J,M+1,N+1),dtype=complex) 
    
    for m in range (0, M+1):
        for j in range (0,J):
            legterm[j,m,:]=w[j]*(data[j,m])*Pmn[j,m,:] 
    legcoeff=np.sum(legterm,0)
    return legcoeff

def fwd_fft_trunc(data,I,M):
    """Calculates and truncates the fast forward Fourier transform of the input

        :param data: array of dimension IxJ (usually the values of state variables at lat-long coordinates)
        :type data: array of float
        
        :param I: number of longitudes
        :type I: int
        
        :param M: highest wavenumber for associated Legendre polynomials
        :type M: int

        :return datam: Fourier coefficients 0 through M
        :rtype datam: array of complex
        """
    datam=np.fft.fft(data/I,I,1)[:,0:M+1]
    
    return datam

def invrs_leg(legcoeff,I,J,M,N,Pmn):
    """Calculates the inverse Legendre transform function
    
    :param legcoeff: Legendre coefficients
    :type legcoeff: array of complex
    
    :param J: number of latitudes
    :type J: int
    
    :param M: highest wavenumber for associated Legendre polynomials
    :type M: int
    
    :param Pmn: associated legendre functions evaluated at the Gaussian latitudes mus  up to wavenumber M
    :type Pmn: array of float
    
    :return: transformed spectral coefficients
    :rtype: array of complex
    """
    
    approxXim=np.zeros((J,I),dtype=complex)
    approxXimPos=np.zeros((J,M+1),dtype=complex) 
    approxXimNeg=np.zeros((J,M),dtype=complex) 
    for m in range (0, M+1):

        approxXimPos[:,m]=np.matmul(Pmn[:,m,m:N+1],(legcoeff[m,m:N+1]))
        
    #use symmetry of the associated Legendre polynomials to compute negative coefficients
        if m !=0:
            negm=-m
            negPmn=((-1)**m)*Pmn[:,m,m:N+1]
            negXileg=((-1)**m)*np.conj(legcoeff[m,m:N+1])
            approxXimNeg[:,negm]=np.matmul(negPmn,negXileg)
    approxXim[:,0:M+1]=approxXimPos

    approxXim[:,I-M:I]=approxXimNeg
    return approxXim

def invrs_fft(approxXim,I):
    """Calculates the inverse Fourier transform
    
    :param approxXim: Fourier coefficients
    :type approxXim: array of complex
    
    :param I: number of longitudes
    :type I: integer

    :return: long-lat coefficients
    :rtype: array of complex
    """
    approxXinew=np.fft.ifft(I*approxXim,I,1);
    return approxXinew

def invrsUV(deltamn,etamn,fmn,I,J,M,N,Pmn,Hmn,tstepcoeffmn,marray):
    
    """
    Computes the wind velocity from the values of vorticity and divergence. This is a diagnostic relationship.
    For details, see Hack and Jakob (1992) equations (5.24)-(5.25).

    Parameters
    ----------
    :param deltamn: Fourier coefficients of divergence
    :type deltamn: array of complex
    
    :param etamn: Fourier coefficients of vorticity
    :type etamn: array of complex

    :param fmn: spectral coefficients of the Coriolis force
    :type etamn: array of float
    
    :param I: number of longitudes 
    :type I: int
    
    :param J: number of latitudes
    :type J: int
    
    :param M:  highest wavenumber for associated Legendre polynomials
    :type M: int
    
    :param N:  highest degree of associated Legendre polynomials
    :type N: int
    
    :param Pmn: values of the associated Legendre polynomials at Gaussian 
    latitudes mus up to wavenumber M
    :type Pmn: array of float
    
    :param Hmn: values of the associated Legendre polynomial derivatives at Gaussian 
    latitudes up to wavenumber M
    :type Hmn: array of float
        
    :param tstepcoeffmn: coefficient to scale spectral components
    :type tstepcoeffmn: array of float
    
    :param marray: array to multiply a quantity by a factor of m ranging from 0 through M.
    :type marray: array of float

    Returns
    -------
    
     :return:
        - Unew 
                Zonal velocity component
        - Vnew 
                Meridional velocity component

    :rtype: array of float

    """
    
    #do not sum over n=0 (see Hack and Jakob 1992 equations 5.24-5.25)
    deltamn[:,0]=0
    etamn[:,0]=0
        
    newUm1=invrs_leg((1j)*np.multiply(np.multiply(marray,deltamn),tstepcoeffmn), I,J, M, N, Pmn)
    newUm2=invrs_leg(np.multiply(etamn-fmn,tstepcoeffmn), I,J, M, N, Hmn)
    
    newVm1=invrs_leg((1j)*np.multiply(np.multiply(marray,etamn-fmn),tstepcoeffmn), I,J, M, N, Pmn)
    newVm2=invrs_leg(np.multiply(deltamn,tstepcoeffmn), I,J, M, N, Hmn)
    

    Unew=-invrs_fft(newUm1-newUm2, I)
    Vnew=-invrs_fft(newVm1+newVm2, I)
    return Unew, Vnew

def diagnostic_eta_delta(Um,Vm, fmn,I,J,M,N,Pmn,Hmn,w,tstepcoeff,mJarray,dt):
    """

    Computes vorticity and divergence from zonal and meridional wind fields. This is a diagnostic relationship.
    For details, see Hack and Jakob (1992) equations (5.26)-(5.27).

   
    
    :param Um: Fourier coefficient of zonal winds
    :type Um: array of float
    :param Vm: Fourier coefficient of meridional winds
    :type Vm: array of float
    :param fmn: spectal coefficients of the Coriolic force
    :type fmn: array of float
    :param I: number of longitudes
    :type I: int
    :param J: number of Gaussian latitudes
    :type J: int
    :param M:  highest wavenumber for associated Legendre polynomials
    :type M: int
    
    :param N:  highest degree of associated Legendre polynomials
    :type N: int
    
    :param Pmn: values of the associated Legendre polynomials at Gaussian 
    latitudes mus up to wavenumber M
    :type Pmn: array of float64
    
    :param Hmn: values of the associated Legendre polynomial derivatives at Gaussian 
    latitudes up to wavenumber M

    :type Hmn: array of float
    :param w: Gauss Legendre weights
    
    :type w: array of float
    :param tstepcoeff: a coefficient for time-stepping of the form 2dt/(a(1-mus^2))
    from Hack and Jakob (1992)
    :type tstepcoeff: array of float
    :param mJarray: coefficients equal to m=0,1,...,M
    :type mJarray: array of float
    :param dt: time step, in seconds
    :type dt: float
     Returns
        -------
        
         :return:
            - neweta 
                    Absolute vorticity
            - newdelta 
                    Divergence
            - etamn 
                    Spectral coefficients of absolute vorticity
            - deltamn 
                    Spectral coefficients of divergence
    
        :rtype: array of float

    """

    coeff=tstepcoeff/(2*dt)
    etacomp1prep=np.multiply(np.multiply(coeff,(1j)*mJarray),Vm)
    etacomp2prep=np.multiply(coeff,Um)
    
    zetamn=fwd_leg(etacomp1prep, J, M, N, Pmn, w)+fwd_leg(etacomp2prep, J, M, N, Hmn, w)
    
    etamn=zetamn+fmn
    
    deltacomp1prep=np.multiply(np.multiply(coeff,(1j)*mJarray),Um)
    deltacomp2prep=np.multiply(coeff,Vm)
    
    deltacomp1=fwd_leg(deltacomp1prep, J, M, N, Pmn, w)
    deltacomp2=fwd_leg(deltacomp2prep, J, M, N, Hmn, w)
    
    deltamn=deltacomp1-deltacomp2
    
    newdeltam=invrs_leg(deltamn, I,J, M, N, Pmn)
    newdelta=invrs_fft(newdeltam, I)
    
    newetam=invrs_leg(etamn, I,J, M, N, Pmn)
    neweta=invrs_fft(newetam, I)
   
    return neweta,newdelta,etamn,deltamn
    
    