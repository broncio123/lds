import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import ipywidgets as widgets

def draw_ld(fig, axis, LD, grid_parameters, subplot_title, interactive, cmap='viridis'):
    """
    Draws a Lagrangian descriptor contour plot and a contour plot showing the magnitude of its gradient field.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure where contour plot will be drawn.

    axis : matplotlib.axes._subplots.AxesSubplot
        Axis handle of subplot in figure.

    LD : ndarray, shape(n, )
        Array to be plotted.

    grid_parameters : list of 3-tuples of floats
        Limits and size of mesh per axis.

    subplot_title : string
        Subplot title.

    interactive : bool
        True allows interactively adjusting the gradient contouor plot minimum and maximum.

    cmap: string
        Name of matplotlib colormap for plot.
    """

    #axes
    if type(grid_parameters) == dict:
        #n-DoF systems
        slice_parameters = grid_parameters['slice_parameters'] # 2n-D grid
        dims_slice = np.array(grid_parameters['dims_slice'])
        if len(dims_slice) == 2:
            labels = ['$x$', '$p_x$']
        elif len(dims_slice) == 4:
            labels = np.array(['$x$','$y$','$p_x$','$p_y$'])
        elif len(dims_slice) == 6:
            labels = np.array(['$x$','$y$','$z$','$p_x$','$p_y$','$p_z$'])
        else:
            labels = np.array(['']*len(dims_slice))
        slice_axes_labels = labels[dims_slice==1]
    else:
        #1-DoF systems
        slice_parameters = grid_parameters # 2-D grid
        slice_axes_labels = ['$x$', '$p_x$']

    ax1_min, ax1_max, N1 = slice_parameters[0]
    ax2_min, ax2_max, N2 = slice_parameters[1]

    points_ax1 = np.linspace(ax1_min, ax1_max, N1)
    points_ax2 = np.linspace(ax2_min, ax2_max, N2)

    # plot
    n_levels = 100
    vmin = np.nanmin(LD)
    vmax = np.nanmax(LD)
    interact_step = (vmax-vmin)/n_levels

    con1 = axis.contourf(points_ax1, points_ax2, LD, cmap=cmap, levels=n_levels)
    if interactive:
        @widgets.interact(clim_min=(vmin, vmax-interact_step, interact_step),clim_max=(vmin+interact_step, vmax, interact_step))
        def update(clim_min=vmin,clim_max=vmax):
            clim_max = max(clim_min+interact_step, clim_max)
            con1.set_clim(clim_min,clim_max)

    axins = inset_axes(axis,
               width="5%",
               height="100%",
               loc='lower left',
               bbox_to_anchor=(1.025, 0., 1, 1),
               bbox_transform=axis.transAxes,
               borderpad=0
               )
    ticks_gradient = np.linspace(vmin,vmax, 11)
    fig.colorbar(con1, cax=axins, ticks=ticks_gradient, format='%.1f', orientation='vertical')

    axis.set_title(subplot_title)
    axis.set_xlabel(slice_axes_labels[0])
    axis.set_ylabel(slice_axes_labels[1])

def draw_ld_pair(LD, LD_gradient, grid_parameters, plot_title, interactive, cmap_gradient):
    """
    Lagrangian descriptor plot wrapper.

    Parameters
    ----------
    LD : ndarray, shape(n, )
        Array of Lagrangian Descriptor values.

    LD_gradient : ndarray, shape(n, )
        Array of Lagrangian Descriptor gradient values.

    grid_parameters : list of 3-tuples of floats
        Limits and size of mesh per axis.

    plot_title : string
        Plot title.

    interactive : bool
        True allows interactively adjusting the gradient contouor plot minimum and maximum.

    cmap_gradient : string
        Name of matplotlib colormap for gradient contour plot.

    Returns
    -------
       fig : `~.figure.Figure`

       ax : `.axes.Axes` or array of Axes
    """
    fig, ax = plt.subplots(1, 2, figsize=(7.5,3), dpi=130, sharex=True, sharey=True)
    plt.subplots_adjust(top=0.85, bottom=0.13, wspace=0.34)  #margins to accommodate boundary of interactive figure environment
    plt.suptitle(plot_title)

    draw_ld(fig, ax[0], normalise(LD), grid_parameters, 'LD values', interactive=False)
    draw_ld(fig, ax[1], LD_gradient, grid_parameters, 'LD gradient magnitude', interactive, cmap=cmap_gradient)

    return fig, ax

def normalise(A):
    """
    Normalises an array.

    Parameters
    ----------
    A : ndarray, shape(n, )
        Array of input values.

    Returns
    -------
    Normalised array : ndarray, shape(n, ).
    """
    return (A - np.nanmin(A)) / (np.nanmax(A) - np.nanmin(A))

def get_gradient_magnitude(LD):
    """
    Calculates magnitude of the gradient of input array LD.

    Parameters
    ----------
    LD : ndarray, shape(n, )
        Array of input values.

    Returns
    -------
    gradient_magnitude : ndarray, shape(n, )
        Magnitude of the gradient of input array LD.
    """
    gradient_x, gradient_y = np.gradient(LD)
    gradient_magnitude = np.sqrt(gradient_x**2 + gradient_y**2)
    return normalise(gradient_magnitude)

def draw_all_lds(LD_forward, LD_backward, grid_parameters, tau=np.nan, p_value=np.nan, interactive=False):
    """
    Draws the forward, backward and total Lagrangian descriptor contour plots and a contour plots showing the magnitude of its gradient field.

    Parameters
    ----------
    LD_forward : ndarray, shape(n, )
        Array of Lagrangian Descriptor values in forward time.

    LD_backward : ndarray, shape(n, )
        Array of Lagrangian Descriptor values in backward time.

    grid_parameters : list of 3-tuples of floats
        Limits and size of mesh per axis.

    tau : float, optional
        Time of integration.
        Default is np.nan.

    p_value : float, optional
        Exponent in Lagrangian descriptor definition.
        Default is np.nan.

    interactive : bool, optional
        True allows interactively adjusting the gradient plot minimum and maximum.
        Default is False.

    Returns
    -------
        List of tuples of the form (fig, ax).
    """

    # Prepare method name
    if np.isnan(p_value):
        t_final = np.nan
    else:
        if p_value == 2:
            str_method = 'Arclength LD'
        elif p_value >= 1:
            str_method = r'p-norm LD$(p={'+str(p_value)+r'})$'
        elif p_value == 0:
            str_method = 'Action-based LD'
        elif p_value < 1:
            str_method = r'LD$_{'+str(p_value)+r'}$'
        else:
            str_method = ''
        t_final=np.abs(tau)

    # Plot LDs
    plot_handles=[]

    if len(LD_forward)>0:
        if np.isnan(t_final):
            plot_title=''
        else:
            plot_title = r'Forward {}, $\tau={}$'.format(str_method,t_final)
        LD_forward_gradient = get_gradient_magnitude(LD_forward)
        plot_tuple = draw_ld_pair(LD_forward, LD_forward_gradient, grid_parameters, plot_title, interactive, 'Blues')
        plot_handles.append(plot_tuple)

    if len(LD_backward)>0:
        if np.isnan(t_final):
            plot_title=''
        else:
            plot_title = r'Backward {}, $\tau={}$'.format(str_method,t_final)
        LD_backward_gradient = -get_gradient_magnitude(LD_backward)
        plot_tuple = draw_ld_pair(LD_backward, LD_backward_gradient, grid_parameters, plot_title, interactive, 'Reds_r')
        plot_handles.append(plot_tuple)

    if len(LD_forward)>0 and len(LD_backward)>0:
        if np.isnan(t_final):
            plot_title=''
        else:
            plot_title = r'Total {}, $\tau={}$'.format(str_method,t_final)
        plot_tuple = draw_ld_pair(LD_backward+LD_forward, LD_forward_gradient+LD_backward_gradient, grid_parameters, plot_title, interactive, 'RdBu')
        plot_handles.append(plot_tuple)

    plt.show()

    return plot_handles

__author__ = 'Broncio Aguilar-Sanjuan, Victor-Jose Garcia-Garrido, Vladimir Krajnak, Shibabrat Naik'
__status__ = 'Development'
