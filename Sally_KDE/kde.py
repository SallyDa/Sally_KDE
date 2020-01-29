
import numpy as np
import copy

def get_x_array():
    '''
    makes an array of x values that are evenly spaced when plotted on a log scale
    (from 1 to ~ 3000)
    '''
    logx100 = range(0, 700) # this will give me 1400 different x positions
    x_lst_p = [10**(a/200.) for a in logx100] # 700/200 ~ log10(3000)
    x_lst_n = [b*-1 for b in x_lst_p[::-1]]
    x_lst = x_lst_n + x_lst_p
    x_array = np.asarray(x_lst)
    return x_array

def kde_plot(map3, logx=True, kern='gaussian', normed=True, ax=None, extrasmooth=1, alpha=1, addlegend=True):
    '''
    makes a kernel density estimation plot using a Gaussian or box kernel
    also returns an array of the coordinates that are plotted
    '''
    bmax = 3000 # maximum flux bin
    bmin = 5 # plot graph stating at 5 G on the x axis

    # flux values for all pixels in the selected region
    all_pix = map3.data[np.where(map3.mask == False)]
    if len(all_pix) == 0: # no selected region, take all pixels
        map3.mask = np.zeros(map3.data.shape)
        all_pix = map3.data[np.where(map3.mask == False)]
    try:
        pixl = all_pix.tolist()
    except:
        all_float_pix = [float(pix) for pix in all_pix]
        all_float_pix = np.asarray(all_float_pix)
        pixl = all_float_pix.tolist()

    if logx:
        x = get_x_array()
    else:
        x = np.arange(-bmax, bmax, 5)

    if kern == 'gaussian':
        y_coords = my_kde(x, pixl, extrasmooth=extrasmooth)

    elif kern == 'box':
        y_coords = my_box_kde(x, pixl)

    y_coordsn = copy.deepcopy(y_coords)
    if normed:
        num_pix = len(all_pix)
        y_coords /= float(num_pix)
        y_coordsn /= float(num_pix)

    fig_data = np.array((x, y_coords))

    if ax == None: # return info without image
        return fig_data, all_pix
    
    fig_k = ax.plot()
    hist_pos_label, hist_neg_label = 'Positive', 'Negative'

    ax.plot(x, y_coords, color='r', label=hist_pos_label, alpha=alpha)
    ax.plot(-x, y_coordsn, color='b', label=hist_neg_label, alpha=alpha)

    ax.set_yscale('log')
    if normed:
        ax.set_ylim(10**-6, 10**-1)
        ax.set_ylabel('Probability density')
    else:
        ax.set_ylim(10**-2, 10**3)
        ax.set_ylabel('Pixel frequency')

    if logx == True:
        ax.set_xscale('log')
    ax.set_xlim(bmin, bmax)
    ax.set_xlabel('Magnetic flux density [Gauss]')
    
    if addlegend:
        ax.legend(loc=1)
    
    return fig_k, fig_data, all_pix


def my_kde(x, mylist, extrasmooth=1):
    '''
    applies a Gaussain with standard deviation sd to each of a list of data points (mylist)
    and sums these Gaussian distributions
    returns an array of y values corresponding to an array of x values
    '''
    y = 0
    for dat_p in mylist:
        # use a variable standard deviation std
        if abs(dat_p) > 25:
            std = abs(dat_p)/10. * extrasmooth
        else:
            std = 2.5 * extrasmooth
        yadd = 1 / (std*np.sqrt(2*np.pi)) * np.exp(-(x-dat_p)**2 / (2*std**2))
        y += yadd
    return y

def my_box_kde(x, mylist):
    '''
    applies a box of a certain half width (hbw) to each data point in the list (mylist)
    and sums these box distributions
    returns an array of y values corresponding to an array of x values
    '''
    y = 0
    for datapoint in mylist:
        if abs(datapoint) > 25:
            hbw = abs(datapoint)/5.
        else:
            hbw = 5
        yadd = box((x-datapoint), halfboxwidth=hbw) / (2 * hbw)
        y += yadd
    return y

def box(x, halfboxwidth=1):
    '''
    the output array is made up of 1s, 0s and 0.5s
    each point of the output array corresponds to a point in the input array
    1 if the input array value is inside the box (the value is less than halfboxwidth)
    and 0 if it is outside of the box
    '''
    newarray = np.zeros(x.shape)
    newarray[np.where(abs(x) < halfboxwidth)] = 1
    newarray[np.where(abs(x) == halfboxwidth)] = 0.5
    return newarray
