#################################################################
#                                                               #
#          SPECTRUM ANALYST                                     #
#                        version: 1.0 - Feb - 2020              #
# @author: Sergio Lins               sergio.lins@roma3.infn.it  #
#################################################################

import sys,os
import numpy as np
import math
from math import factorial
import matplotlib.pyplot as plt
import EnergyLib

"""
1st; calculate background with peakstrip()
2nd; calibrate the energy axis
3rd; verify configuration parameters - configdict and lookup variables -
4th; run getpeakarea() 
"""
global __FANO__
global __NOISE__

__NOISE__, __FANO__ = 80, 0.114

def updatefano(x,y):
    global __FANO__, __NOISE__
    __FANO__, __NOISE__ = x, y
    return

def function(ydata,x):
    
    """ Returns x index value of ydata array """
    
    return ydata[x]

def dif2(ydata,x,gain):
    
    """ Complementary function to getdif2 """
    
    value = (function(ydata, x + 2*gain) - 2*function(ydata, x + gain)\
            + function(ydata, x)) / (gain * gain)
    return value

def getdif2(ydata,xdata,gain):
    
    """ Returns the second differential of ydata """
    
    xinterval = np.pad(xdata,2,'edge')
    yinterval = np.pad(ydata,2,'edge')
    dif2curve = np.zeros([len(yinterval)])
    for x in range(len(xinterval)-2):
        difvalue = dif2(yinterval,x,1)
        dif2curve[x] = difvalue
    dif2curve = dif2curve[1:-3]
    return dif2curve

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
   
    """ This function was taken from https://gist.github.com/krvajal 
    and compared to scipy.signal package, an affidable and renown package
    for signal analysis. """
    """ References
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688 """

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError as msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

def strip(an_array,cycles,width):

    """
    Strips the peaks contained in the input array following an
    iteractive peak clipping method (van Espen, Chapter 4 in Van Grieken, 2002)
    INPUT:
        an_array; np.array
        cycles; int
        width; int
    OUTPUT:
        an_array; np.array
    """

    ##################################################################
    # W IS THE WIDTH OF THE FILTER. THE WINDOW WILL BE (1*W)+1       #
    # W VALUE MUST BE LARGER THAN 2 AND ODD, SINCE 3 IS THE MINIMUM  #
    # SATISFACTORY POLYNOMIAL DEGREE TO SMOOTHEN THE DATA            #
    ##################################################################

    size = an_array.shape[0]
    for k in range(cycles):
        if k >= cycles-8:
            width = int(width/np.sqrt(2))
        for l in range(0, size):
            if l-width < 0: low = 0
            else: low = l-width
            if l+width >= size: high = size-1
            else: high = l+width
            m = (an_array[low] + an_array[high])/2
            if an_array[l] > m: an_array[l] = m
    return an_array

def peakstrip(an_array,cycles,width,*args):
    
    """
    Calculates the continuum/background contribution of the input spectrum following
    the SNIPBG method described in van Espen, Chapter 4 in Van Grieken, 2002.
    INPUT:
        an_array; np.array
        cycles; int
        width; int
        args; list (cycles, width, sav_gol order and sav_gol window can be fully    
        given as a list instead)
    OUTPUT:
        snip_bg; np.array
    """

    TEST = False

    #initialize snip_bg array
    snip_bg = np.zeros(an_array.shape[0])

    #square root the data
    sqr_data = an_array**0.5

    #apply savgol filter to the square root data
    if len(args) > 0:
        savgol_window,order = args[0],args[1]
        try: 
            #smooth_sqr = scipy.signal.savgol_filter(sqr_data,savgol_window,order)
            smooth_sqr = savitzky_golay(sqr_data,savgol_window,order) 

        except:
            raise ValueError
    else: 
        smooth_sqr = savitzky_golay(sqr_data,width,3)
    
    for i in range(smooth_sqr.shape[0]): 
        if smooth_sqr[i] < 0: smooth_sqr[i] = 0
    
    #strip peaks
    snip_bg = strip(smooth_sqr,cycles,width)

    #transform back
    snip_bg **= 2
   
    return snip_bg

def sigma(energy):
    global __NOISE__, __FANO__
    return math.sqrt(((__NOISE__/2.3548)**2)+(3.85*__FANO__*energy))

def shift_center(xarray,yarray):
    
    """ Returns the highest value and its index in yarray and its corresponding value
    in xarray """

    ymax = yarray.max()
    y_list = yarray.tolist()
    idx = y_list.index(ymax)
    return xarray[idx],ymax,idx

def setROI(lookup,xarray,yarray,localconfig):
    
    """
    INPUT: 
        lookup; eV energy (int)
        xarray; np.array
        yarray; np.array
        localconfig; dict
    OUTPUT: 
        low_idx; int
        high_idx; int
        peak_center; int
        isapeak; bool
        - indexes corresponding to 2*FWHM of a gaussian centered
        at eV energy position (int, int, int, bool) 
    """
    
    lookup = int(lookup)
    peak_corr = 0
    isapeak = True
    
    if localconfig.get('bgstrip') == "SNIPBG":
        yarray  = savitzky_golay(yarray,5,3)
    
    for peak_corr in range(2):
        FWHM = 2.3548 * sigma(lookup)
        lowx = (lookup - (FWHM))/1000
        highx = (lookup + (FWHM))/1000
        
        idx = 0
        while xarray[idx] <= lowx:
            idx+=1
        lowx_idx = idx-3
        while xarray[idx] <= highx:
            idx+=1
        highx_idx = idx+3
        
        ROIaxis = xarray[lowx_idx:highx_idx]
        ROIdata = yarray[lowx_idx:highx_idx]
        shift = shift_center(ROIaxis,ROIdata)
        
        if 1.10*(-FWHM/2) < (shift[0]*1000)-lookup < 1.10*(FWHM/2):
            lookup = shift[0]*1000
            peak_corr = 0
        else: 
            lookupcenter = int(len(ROIaxis)/2)
            shift = (0,0,lookupcenter)
            isapeak = False
        
    lowx_idx = lowx_idx + 2
    highx_idx = highx_idx - 3
    peak_center = shift[2]-3
    return lowx_idx,highx_idx,peak_center,isapeak

def getpeakarea(lookup,data,energyaxis,continuum,localconfig,RAW,usedif2,dif2):

    """
    Calculates the netpeak area of a given lookup energy.
    INPUT:
        lookup; int (theoretical peak center eV energy)
        data; np.array (counts)
        energyaxis; np.array (calibrated energy axis KeV)
        continuum; np.array (background counts)
        localconfig; dict (condiguration)
        RAW; np.array (counts, same as data)
        usedif2; bool
        dif2; np.array (second differential of counts)
    OUTPUT:
        area; float
        idx; list (containing setROI function outputs)
    """

    Area = 0
   
    idx = setROI(lookup,energyaxis,data,localconfig)
    isapeak = idx[3]
    xdata = energyaxis[idx[0]:idx[1]]
    ydata = data[idx[0]:idx[1]]
    original_data = RAW[idx[0]:idx[1]] 
    ROIbg = continuum[idx[0]:idx[1]]
    
    if usedif2 == True: 
            dif2 = getdif2(data,energyaxis,1)
            for i in range(len(dif2)):
                if dif2[i] < -1: dif2[i] = dif2[i]
                elif dif2[i] > -1: dif2[i] = 0

    ######################################
    # SIGNAL TO NOISE PEAK TEST CRITERIA #  
    ######################################
    # After Bright, D.S. & Newbury, D.E. 2004

    if isapeak == True: 
        if original_data.sum() - ROIbg.sum() < 3*math.sqrt(abs(ROIbg.sum())): isapeak = False
    
    ##########################
    # 2ND DIFFERENTIAL CHECK #
    ##########################
    if usedif2 == True and isapeak == True:
        sliced_dif2 = dif2[idx[0]:idx[1]]
    
        # checks second differential and calculates the net area - background
        if sliced_dif2[idx[2]] < 0 or sliced_dif2[idx[2]+1] < 0 or sliced_dif2[idx[2]-1] < 0:
            if ROIbg.sum() < ydata.sum(): Area += ydata.sum() - ROIbg.sum()
    
    ##########################
    return Area,idx

def getdata(mca):
    
    """ Extract the data contained in spectrum files 
    INPUT:
        mca; path
    OUTPUT:
        Data; 1D-array """

    name = str(mca)
    name = name.split("\\")[-1]
    name = name.replace('_',' ')
    
    # custom MC generated files
    if 'test' in name or 'obj' in name or 'newtest' in name:
        Data = []
        datafile = open(mca)
        lines = datafile.readlines()
        for line in lines:
            line = line.split()
            try: counts = float(line[1])
            except: counts = float(line[0])
            counts = counts * 10e3
            Data.append(counts)
        Data = np.asarray(Data)

    # this works for mca extension files
    else:
        ObjectData=[]
        datafile = open(mca)
        line = datafile.readline()
        line = line.replace("\r","")
        line = line.replace("\n","")

        # AMPTEK files start with this tag
        if "<<PMCA SPECTRUM>>" in line:
            while "<<DATA>>" not in line:
                line = datafile.readline()
                if line == "": break
            line = datafile.readline()
            while "<<END>>" not in line:
                try: ObjectData.append(int(line))
                except ValueError as exception:
                    datafile.close()
                    raise exception.__class__.__name__
                line = datafile.readline()
                if line == "": break

        # Works if file is just counts per line
        elif line.isdigit():
            while "<<END>>" not in line:
                ObjectData.append(int(line))
                line = datafile.readline()
                if line == "": break
        
        # if file has two columns separated by space or tab
        elif "\t" in line or " " in line: 
            while "<<END>>" not in line:
                counts = line.split("\t")[-1]
                if counts.isdigit(): ObjectData.append(int(counts))
                else:
                    counts = line.split(" ")[-1]
                    ObjectData.append(int(counts))
                line = datafile.readline()
                if line == "": break

        Data = np.asarray(ObjectData)
    datafile.close()
    return Data

def linregress(x, y, sigmay=None, full_output=False):

    # DISCLAIMER #
    # function extracted from PyMca5.PyMcaMath.fitting.RateLaw script

    """
    Linear fit to a straight line following P.R. Bevington:
    "Data Reduction and Error Analysis for the Physical Sciences"
    Parameters
    ----------
    x, y : array_like
        two sets of measurements.  Both arrays should have the same length.
    sigmay : The uncertainty on the y values
    Returns
    -------
    slope : float
        slope of the regression line
    intercept : float
        intercept of the regression line
    r_value : float
        correlation coefficient
    if full_output is true, an additional dictionary is returned with the keys
    sigma_slope: uncertainty on the slope
    sigma_intercept: uncertainty on the intercept
    stderr: float
        square root of the variance
    
    """

    x = np.asarray(x, dtype=np.float).flatten()
    y = np.asarray(y, dtype=np.float).flatten()
    N = y.size
    if sigmay is None:
        sigmay = np.ones((N,), dtype=y.dtype)
    else:
        sigmay = np.asarray(sigmay, dtype=np.float).flatten()
    w = 1.0 / (sigmay * sigmay + (sigmay == 0))

    n = S = w.sum()
    Sx = (w * x).sum()
    Sy = (w * y).sum()    
    Sxx = (w * x * x).sum()
    Sxy = ((w * x * y)).sum()
    Syy = ((w * y * y)).sum()
    # SSxx is identical to delta in Bevington book
    delta = SSxx = (S * Sxx - Sx * Sx)

    tmpValue = Sxx * Sy - Sx * Sxy
    intercept = tmpValue / delta
    SSxy = (S * Sxy - Sx * Sy)
    slope = SSxy / delta
    sigma_slope = np.sqrt(S /delta)
    sigma_intercept = np.sqrt(Sxx / delta)

    SSyy = (n * Syy - Sy * Sy)
    r_value = SSxy / np.sqrt(SSxx * SSyy)
    if r_value > 1.0:
        r_value = 1.0
    if r_value < -1.0:
        r_value = -1.0

    if not full_output:
        return slope, intercept, r_value

    ddict = {}
    # calculate the variance
    if N < 3:
        variance = 0.0
    else:
        variance = ((y - intercept - slope * x) ** 2).sum() / (N - 2)
    ddict["variance"] = variance
    ddict["stderr"] = np.sqrt(variance)
    ddict["slope"] = slope
    ddict["intercept"] = intercept
    ddict["r_value"] = r_value
    ddict["sigma_intercept"] = np.sqrt(Sxx / SSxx)
    ddict["sigma_slope"] = np.sqrt(S / SSxx)
    return slope, intercept, r_value, ddict

def getcalibration(cfile=None):
    param = []
    
    if configdict['calibration'] == 'manual' and cfile==None:

        """ These parameters are to calibrate specific Monte Carlo simulated spectrum,
        change accordingly to your spectrum files """

        param = [[41,1.6],[152,6.4],[174,6.92],[202, 8.04],[280, 11.16],[555,22.16]] 

    elif configdict['calibration'] == 'manual' and cfile!=None:
        calib_file = open(cfile,"r")
        line = calib_file.readline()
        while line != "":
            line = line.replace("\r","")
            line = line.replace("\n","")
            line = line.replace("\t"," ")
            anchor = line.split(" ")
            for number in anchor:
                if not number.isdigit():
                    line = calib_file.readline()
                    break
            param.append(anchor)
            line = calib_file.readline()
        calib_file.close()

    elif configdict['calibration'] == 'from_source':

        """ Gets the calibration parameters from the spectrum file header. Tested with 
        amptek mca files """

        param = []
        mca_file = open(spectrum_path,'r')
        line = mca_file.readline()
        while line != "":
            while "<<CALIBRATION>>" not in line:
                line = mca_file.readline()
                if line == "": break
            while "<<DATA>>" not in line:
                line = mca_file.readline()
                if line == "": break
                line=line.replace('\r','')
                line=line.replace('\n','')
                line=line.replace('\t',' ')
                aux = line.split()
                try: param.append([int(aux[0]),float(aux[1])])
                except: pass 
            line = mca_file.readline()
        mca_file.close()
        if param == []: 
            raise ValueError("Couldn't fetch calibration from source!")
        for parameter in param:
            if len(param) <= 1: 
                raise ValueError("Need at least two calibration points!")
            elif parameter[0] < 0 or parameter[1] < 0: 
                raise ValueError("Cant receive negative values!")
            elif parameter[0] == 0 or parameter[1] == 0:
                raise ValueError("Cant receive zero as parameter!")
        else: pass
    else: 
        raise ValueError
    return param

def calibrate(nchan,cfile=None):
    
    """
    Returns the energy axis and gain of the calibrated axis
    The parameter are taken from config.cfg is calibration is set to manual
    or from the mca files is calibration is set to from_source
    """

    param = getcalibration(cfile)
    x=[]
    y=[]
    for i in range(len(param)):
        for k in range(len(param[i])):
            x.append(param[i][0])
            y.append(param[i][1])
    x=np.asarray([x])
    y=np.asarray([y])
    coefficients=list(linregress(x,y))
    GAIN=coefficients[0]
    B=coefficients[1]
    R=coefficients[2]
    curve = []
    for i in range(nchan):
        curve.append((GAIN*i)+B)
    curve = np.asarray(curve)
    return curve,GAIN

if __name__.endswith('__main__'):         

    """ Initiate variables """

    spectrum_name = "newtest_1.txt"
    configdict = {}
    configdict["calibration"] = "manual"
    configdict["bgstrip"] = "SNIPBG"
    configdict["fill"] = True
    configdict["write"] = False
    configdict["element"] = None
    
    """ Get user input data """

    for arg in sys.argv:
        if arg.startswith("-"): 
            if arg == "-s": 
                try: spectrum_name = sys.argv[sys.argv.index(arg)+1]
                except: raise ValueError("No file name given")
            elif arg == "-cal":
                found = False
                try: cfile = sys.argv[sys.argv.index(arg)+1]
                except: raise ValueError("No calibration file given")
                for parent, child, files in os.walk(os.getcwd()):
                    if cfile in files: 
                        configdict["calibration"] = "manual"
                        cfile = os.path.join(parent, cfile)
                        found = True
                        break
                if found == False: raise FileNotFoundError("Could not locate file {}".format(cfile))
            elif arg == "-e":
                try: configdict["element"] = sys.argv[sys.argv.index(arg)+1]
                except: raise ValueError("No element given")
            elif arg == "-nobg": configdict["bgstrip"] = None
            elif arg == "-nofill": configdict["fill"] = False
            elif arg == "-fs": configdict["calibration"] = "from_source"
            elif arg == "-nosi": updatefano(0.129, 80)
            elif arg == "-w": 
                configdict["write"] = True
                try: output_name = sys.argv[sys.argv.index(arg)+1]
                except: raise ValueError("No file name given")
                if not output_name.endswith(".txt"):
                    raise IOError("Can't create non-txt file {}".format(output_name))
            else: raise ValueError("Command {} not recognized".format(arg))
    try: lookup = EnergyLib.Energies[EnergyLib.ElementList.index(configdict["element"])]*1000
    except: lookup = 100

    """ Build spectrum file path """

    spectrum_path = os.path.join(os.getcwd(),spectrum_name)
    
    """ Create data """
    data = getdata(spectrum_path) 
    energyaxis, GAIN = calibrate(len(data),cfile)
    if configdict["bgstrip"] == "SNIPBG":

        """ peakstrip parameters: 
        cycles, window, sav_gol_width, sav_gol_order
        if -nofill is input by the user, bg variable
        is set as a zero array """

        bg = peakstrip(data,24,5,5,3) 
    elif configdict["bgstrip"] == None:
        bg = np.zeros([len(data)])
    
    """ Calculate peak area and get indexes """

    if configdict["element"] != None:
        net_peak, peak_indexes = getpeakarea(
                lookup, 
                data, 
                energyaxis, 
                bg, 
                configdict, 
                data, 
                True,
                0)
       
        """ Creates an array containing the peak area with the same length 
        of the data array, so they can be plot together """

        lookup_data, lookup_bg = np.zeros([len(data)]), np.zeros([len(data)])
        for i in range(len(data)):
            if peak_indexes[1] >= i and peak_indexes[0] <= i:
                lookup_data[i] += data[i]
                lookup_bg[i] += bg[i]

    """ If -w is input by the user, the cleaned spectrum is written to dosk with
    input name """

    if configdict["write"] == True:
        clean = data-bg
        try: output_file = open(os.path.join(os.getcwd(),output_name), "+w")
        except: raise OSError("Cannot create {}".format(os.path.join(os.getcwd(),output_name)))
        output_file.write("<<START>>\r")
        for i in range(len(data)):
            output_file.write("{}\r".format(data[i]-bg[i]))
        output_file.write("<<END>>")
        output_file.close()

    """ Plots data to screen """

    fig, ax = plt.subplots()
    plt.semilogy(energyaxis, data, label=spectrum_name, color="blue")
    plt.semilogy(energyaxis, bg, label="Continuum", color="yellow")
    if configdict["write"] == True:
        plt.clf()
        plt.semilogy(energyaxis, clean, label="Cleaned Spec", color="orange")
        plt.semilogy(energyaxis, data, label=spectrum_name, color="blue")
        plt.semilogy(energyaxis, bg, label="Continuum", color="yellow")
    elif configdict["element"] != None: 
        if peak_indexes[3] == True: 
            plt.semilogy(
                    energyaxis, 
                    lookup_data, 
                    label="Peak {}".format(net_peak), 
                    color="red")
            if configdict["fill"] == True \
            and configdict["write"] == False \
            and peak_indexes[3] == True:
                plt.fill_between(energyaxis, lookup_data, lookup_bg)
            print("Peak found, net area = {}".format(net_peak))
        else:
            plt.plot(
                    (energyaxis[peak_indexes[0]],energyaxis[peak_indexes[0]]),
                    (0,data.max()), 
                    "k--")
            plt.plot(
                    (energyaxis[peak_indexes[1]],energyaxis[peak_indexes[1]]),
                    (0,data.max()), 
                    "k--")
            print("Peak is off-centered. Cannot detected.")
    
    plt.legend()
    ax.set_ylabel("Counts")
    ax.set_xlabel("Energy (KeV)")
    ax.set_title("Spectrum Plot")
    plt.show()


