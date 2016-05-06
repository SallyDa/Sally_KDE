
import copy
import numpy as np

def bfl_loglog(myarray, ax=None, color='g'):
    '''
    input array of x (myarray[0, :]) and y (myarray[1, :]) values
    find and plot best fit line: yl = m*xl + c
    '''
    logarray = np.log10(myarray)
    
    sumy = sum(logarray[1, :])
    sumx = sum(logarray[0, :])
    N = logarray.shape[1]
    
    xy = list()
    xsquared = list()
    for a in range(0, N):
        xypoint = logarray[0, a] * logarray[1, a]
        xsquaredpoint = logarray[0, a] **2
        xy.append(xypoint)
        xsquared.append(xsquaredpoint)
    sumxy = sum(xy)
    sumxsquared = sum(xsquared)
    
    c = (sumy*sumxsquared - sumx*sumxy)/(N*sumxsquared - sumx*sumx)
    m = (N*sumxy - sumx*sumy)/(N*sumxsquared - sumx*sumx)
    
    if ax == None:
        return m, c
    
    fig_b = ax.plot()

    logx = logarray[0, :]
    xl = 10**logx
    yl = 10**(m*logx + c)

    ax.plot(xl, yl, color=color, linestyle='--')
    ax.set_xscale('log')
    ax.set_yscale('log')

    return fig_b, m, c

def split_kde_data_array(myarray):
    '''
    as input array (myarray) use fig_data output from kde_plot function
    split into 2 arrays, one for positive B and one for negative
    '''
    mya = copy.deepcopy(myarray)
    mid = mya.shape[1]/2
    
    parray = mya[:, mid:]
    narray = mya[:, :mid]
    # multiply the negative B values by -1 to get the absolute B
    narray[0, :] *= -1
    # reverse order of B values (and their corresponding slope values)
    # so that they go up in magnitude (similar to parray)
    narray = narray[:, ::-1]
    
    return parray, narray

def slope_plot(myarray, ax=None):
    '''
    make a plot of the slope of the KDE log-log plot
    as input array (myarray) use fig_data output from kde_plot function
    returns also 2 arrays containing the coordinates that make the 2 lines
    '''
    parray, narray = split_kde_data_array(myarray)
    
    Bg_pos = B_slope_array(parray)
    Bg_neg = B_slope_array(narray)

    if ax == None:
        return 'no figure', Bg_pos, Bg_neg

    fig_g = ax.plot()

    B_pos = Bg_pos[0, :]
    g_pos = Bg_pos[1, :]
    ax.plot(B_pos, g_pos, c='r', label='Positive')

    B_neg = Bg_neg[0, :]
    g_neg = Bg_neg[1, :]
    ax.plot(B_neg, g_neg, c='b', label='Negative')

    ax.set_xlim(5, 3000)
    ax.set_xscale('log')
    ax.set_xlabel('Magnetic flux density [Gauss]')
    ax.set_ylabel('Slope of the KDE')
    ax.set_ylim(-6, 2)

    return fig_g, Bg_pos, Bg_neg

def B_slope_array(myarray):
    '''
    input array consists of a set of x and y values
    make a new array containing the log of the values of the input array
    find the slope of these values g: the change in log(y) with change in log(x)
    convert the log(x) values back to x values
    output an array of x and g values
    '''
    logarray = np.log10(myarray)
    xg_array = slope_withx(logarray) # want to find the slope of the log-log plot

    Bg_array = copy.deepcopy(xg_array)
    Bg_array[0, :] = 10**xg_array[0, :]

    return Bg_array

def cut_slope_array(myarray, cutoff=-1):
    '''
    input an array of slope values (myarray[1, :])
    with their corresponding B values (myarray[0, :])
    an array such as Bg_pos from the slope_plot function  
    cut the array removing any coordinates after the slope value drops below -1
    '''
    mytarray = np.transpose(myarray)
    mylst = mytarray.tolist()
    newlst = list()
    for Bval, slope in mylst:
        if slope >= cutoff:
            newlst.append([Bval, slope])
        else:
            break
    return np.transpose(np.asarray(newlst))

def B_at_second_knee(myarray, cutoff=-2):
    '''
    input an array of slope values (myarray[1, :])
    with their corresponding B values (myarray[0, :])
    such as Bg_pos from the slope_plot function
    go along the slope plot from right to left
    until the slope goes above the cutoff (-2)
    return the B value at this point
    '''
    # reverse order of array    
    mya = myarray[:, ::-1]
    mytarray = np.transpose(mya)
    mylst = mytarray.tolist() # mylst goes from large B to small B
    for Bval, slope in mylst:
        if slope >= cutoff:
            B_2ndknee = Bval
            break
    return B_2ndknee

def cut_array(myarray, Bmin, Bmax):
    '''
    input an array of slope values (myarray[1, :])
    with their corresponding B values (myarray[0, :])
    cut an array (such as Bg_pos) to keep slope values between Bmin and Bmax
    return an array of the same form
    '''
    mytarray = np.transpose(myarray)
    mylst = mytarray.tolist()
    newlst = list()
    for Bval, slope in mylst:
        if Bval >= Bmin and Bval <= Bmax:
            newlst.append([Bval, slope])
    return np.transpose(np.asarray(newlst))
    
def mean_slope(myarray):
    '''
    input an array of slope values (myarray[1, :])
    output the mean of these values
    '''
    return np.mean(myarray[1, :])
    
def mean_slope_between_knees(myarray, Bmin, Bmax):
    '''
    input array eg. Bg_pos from the slope_plot function
    find the mean slope between Bmin and Bmax
    '''
    c_array = cut_array(myarray, Bmin, Bmax)
    meanslope = mean_slope(c_array)
    return meanslope

def B_at_steepest_slope(myarray):
    '''
    input an array of slope values (myarray[1, :])
    with their corresponding B values (myarray[0, :])
    B values must be positive
    input array eg. Bg_pos from the slope_plot function
    
    return the B value corresponding to the steepest slope
    '''
    mya = copy.deepcopy(myarray)
    mya[0, :] = np.log10(mya[0, :])
    logB_gg = slope_withx(mya, logx=True)

    Bggl = logB_gg.tolist()
    pos_ggmin = np.argmin(Bggl[1])
    logB_atsg = Bggl[0][pos_ggmin]
    
    return 10**(logB_atsg)

def slope_withx(myarray, logx=True):
    '''
    input array consists of a set of x (myarray[0, :]) and y (myarray[1, :]) values
    the change in y with change in x is calculated (for each step in x)
    output an array of x and g values

    '''
    numx = myarray.shape[1]
    a = 0

    xg_list = list()
    while a < numx - 2:
        array22 = myarray[:, a:a+2]
        xmean_g = getslope(array22, logx=logx)
        xmean = float(xmean_g[0])
        slope2 = float(xmean_g[1])
        xg_list.append((xmean, slope2))
        a += 1
    xg_array = np.transpose(np.asarray(xg_list))
    
    return xg_array

def getslope(array22, logx=False):
    '''
    takes an array of 2 x values (array22[0,:]) and 2 y values (array22[1,:])
    returns an array containing the mean x value
    and the difference of the y values over the difference between the x values
    '''
    if logx == True:
        x_cen = np.log10((10**array22[0, 0] + 10**array22[0, 1])/2.)
    else:
        x_cen = (array22[0, 0] + array22[0, 1])/2.
    
    slope = (array22[1, 1]-array22[1, 0])/(array22[0, 1]-array22[0, 0])
    g_array = np.zeros((2, 1))
    g_array[0, 0] = x_cen
    g_array[1, 0] = slope
    return g_array
    