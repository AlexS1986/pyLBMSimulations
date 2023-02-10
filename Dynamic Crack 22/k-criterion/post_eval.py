"""Module to compute, compare and plot fields around the crack tip.
Reads in the .csv file of the crack data, e.g. position of crack tip or SIF,
and the .vtk files with the simulated fields.

plot_sif reads and plots the SIF with the data from the .csv-file.
Additional metrics can be added via the add_* methods, such as the analytical
value according to Mandal (add_mandal_orig) or mean (add_mean) and median
(add_median), which compute the necessary values with the provided parameters
and the data to the plot.

For the testing model:
    plot_sif            -> SIF
    add_mandal_orig     -> get analytical value
    add_median          -> get median for comparison
"""

import os
from collections import namedtuple

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as lines
from mpl_toolkits.axes_grid1 import ImageGrid

import numpy as np
from scipy.signal import savgol_filter
from scipy.signal import butter, filtfilt
import bottleneck as bn


def _read_local_output(name: str):
    data_path = './output/' + name + '_w.outloc'
    local_data = namedtuple('Fielddata', 'ID x y z w')

    # read point IDs and positions
    with open(data_path, 'r') as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if not line.startswith('#'):
            break
        if 'Point ID' in line:
            idx_idx = i+1
        if 'Point position' in line:
            pos_idx = i+1
    # idx_lines = lines[idx_idx:pos_idx-1]
    idx_lines = [int(l.lstrip('# ').rstrip('\n').strip()) for l in lines[idx_idx:pos_idx-2]]
    pos_lines = [[float(x.strip()) for x in l.lstrip('# ').rstrip('\n').split()] for i, l in enumerate(lines)
                 if i >= pos_idx and l.startswith('#')]

    # idx_lines = np.array(idx_lines)
    # pos_lines = np.array(pos_lines)

    pt_data = np.genfromtxt(data_path)

    time = pt_data[:,0]
    points = []
    for idx, pos, data in zip(idx_lines, pos_lines, pt_data[:,1:].T):
        pt = local_data(idx, *pos, data)
        points.append(pt)

    return time, points


def _read_model_data(name: str):
    """reads data of lattice, crack, etc. from the file header
    i.e. lattice spacing h, crack velocity a_dot, wave speed c_s
    """
    with open('./output/' + name + '_crack_tips.csv', 'r') as f:
        elems = f.readline().rstrip().split(',')
    a_dot = elems[-3].split('=')[-1].strip()
    a_dot = float(a_dot)
    h = elems[-2].split('=')[-1].strip()
    h = float(h)
    c_s = elems[-1].split('=')[-1].strip()
    c_s = float(c_s)

    return h, a_dot, c_s


def _read_crack_data(name: str):
    """reads the columns of data from the file"""
    crack_data = np.genfromtxt(
        './output/' + name + '_crack_tips.csv',
        comments='#',
        skip_header=1,
        # delimiter='  ',   # cannot handle multiple spaces well, let numpy handle it
        names=True,
    )
    return crack_data


def _from_local(name: str,
                mu: float = 1.0, t0: float = 2.0, traction: float = 0.0025):
    time, pts = _read_local_output(name)

    pts.sort(key=lambda pt: pt.x**2 + pt.y**2)
    points = [(pts[i], pts[i+1]) for i in range(0,len(pts), 2)]

    # idx_t = np.where(time == t0)
    # assert idx_t

    max_an = 1.2732     # max. analyt. K, from DataAN.mat (<- scipy.io)
    k_s = traction * np.sqrt(np.pi)

    sif = {'time': time}
    for p1, p2 in points:

        assert abs(p1.x) == abs(p2.x)
        assert abs(p1.y) == abs(p2.y)
        assert p1.x == -p2.x or p1.y == -p2.y

        # line is coaxial to x or y and at h/2 distance: h/2 <= d
        # d = max(abs(p1.x), abs(p1.y))
        d = np.sqrt(p1.x**2 + p1.y**2)

        delta = abs(p1.w - p2.w)

        k = delta * mu / 4 * np.sqrt(2 * np.pi / d)
        sif[d] = np.asarray(k)

    return sif, max_an * k_s


def _from_mat():
    from scipy.io import loadmat
    data = loadmat('DataAN.mat')

    return data['timeAN'][::200], data['KfactAN'][::200]

def _from_grid(name: str, env: tuple, t_idx: int):
    """read data from vtk, return list of points within environment env around tip with data
    for a given time instance
    set w_dot at points from data
    """
    x_min, y_min, x_max, y_max = env

    points, data = _read_grid_data(name, t_idx)

    env_pts = []
    for j, coord in enumerate(points):
        x, y, z = coord
        if (x_min <= x <= x_max) and (y_min <= y <= y_max):       # eliminate AuxPoints from grid
            pt = FieldPoint(x, y, z)
            pt.w = data['dispw'][j]
            pt.w_dot = data['w_dot'][j]
            env_pts.append(pt)

    return env_pts


def _build_env_fields(name: str, ts: list[int],
                      n_max: int = 8, m_max: int = 8, sum_error: bool = False):
    """get arrays of field data for grid and replicate at time instances (indices) ts
    parameters: name <- string, for file names
                n_max <- int,
    """
    def get_points_grid(t_idx):
        """finds the points in the environment around the crack tip"""
        # grid points at t=0 with data and (r, phi) related to tip_grid
        x, y, z = ct_xyz[t_idx]
        tip_grid = FieldPoint(x, y, z)
        env_box = (x - n_max * h, y - m_max * h,
                   x + n_max * h, y + m_max * h)

        # points in environment, already with w, w_dot set
        points_grid = _from_grid(name, env_box, t_idx)
        _update_locale(points_grid, tip_grid, cr_dir)

        return points_grid

    # numerical data
    crack_data = _read_crack_data(name)
    time = crack_data['t']
    sif = crack_data['SIF_1']
    ct_xyz = np.stack((crack_data['x_1'], crack_data['y_1'], crack_data['z_1']), axis=-1)

    # read lattice spacing h, wave speed c_s and crack propagation speed a_dot
    h, a_dot, c_s = _read_model_data(name)

    """replicate points around crack tip with local coordinate system
    assuming tip is initially right in the middle between 4 LatticePoints (e.g. (0, 0, 0))
    and assuming the relative position of tip and points does not change over time,
    i.e. translational invariance w.r.t. lattice spacing
    """
    p_r0 = np.array((0.0, 0.0, 0.0))            # initial coordinates of crack tip
    tip_repli = FieldPoint(*p_r0)               # initial crack tip

    # get the propagation direction of the crack tip
    cr_dir = Vector(*ct_xyz[1]) - Vector(*ct_xyz[0])
    cr_dir = cr_dir / cr_dir.absolute_value()

    # pass on data for plotting in compact form via namedtuple
    ct_plot = [FieldPoint(*ct_xyz[t]) for t in ts]
    num_data = namedtuple('PlotData', 'ct cdir h')
    plot_data = num_data(ct_plot, cr_dir, h)

    # instantiate the point sof the replicate field
    points_repli = [FieldPoint((n + 1 / 2) * h, (m + 1 / 2) * h, 0)
                    for n in range(-n_max, n_max) for m in range(-m_max, m_max)]

    # set initial distance and angle
    _update_locale(points_repli, tip_repli, cr_dir)

    val_repli = np.empty((len(ts), 2 * n_max, 2 * m_max))
    val_grid = np.empty_like(val_repli)
    val_err = np.empty_like(val_repli)
    err_cum = np.empty(len(time))

    t_count = 0
    for idx, k in enumerate(sif):
        # r: replicate field points, g: grid (or lattice) points, e: error at points
        if idx in ts:       # only keep for plot at defined instances
            # get points and sort for grid
            pts_grid = get_points_grid(idx)
            pts_grid.sort(key=lambda p: p.xy_coords())

            # put field values in arrays, get errors
            tmp_val_r = np.array([pt.compute_w_dot(k, a_dot, c_s) for pt in points_repli])
            tmp_val_g = np.array([pt.w_dot for pt in pts_grid])
            tmp_val_e = np.abs(tmp_val_r - tmp_val_g)
            err_cum += np.sum(tmp_val_e)

            #
            val_err[t_count] = tmp_val_e.reshape(2 * n_max, 2 * m_max)
            val_grid[t_count] = tmp_val_g.reshape(2 * n_max, 2 * m_max)
            val_repli[t_count] = tmp_val_r.reshape(2 * n_max, 2 * m_max)
            t_count += 1
        elif sum_error:      # sum of point-wise errors for all times
            tmp_val_r = np.array([pt.w_dot for pt in points_repli])
            tmp_val_g = np.array([pt.w_dot for pt in get_points_grid(idx)])
            tmp_val_e = np.abs(tmp_val_r - tmp_val_g)
            err_cum += np.sum(tmp_val_e)

    err_cum / (n_max * m_max * 4)

    return val_grid, val_repli, val_err, err_cum, plot_data


def _split_time(t_in, t_len):
    """get time interval of length t_len from the end"""
    t0 = t_in[-1] - t_len
    n = np.where(t_in > t0)
    t = t_in[n]

    return t, n


def _line_plot(fig: plt.Figure, y: float, color: str = 'gray'):
    ax = fig.axes[0]
    x = ax.get_xlim()
    ax.plot(
        x,
        [y, y],
        linestyle='--',
        c='gray',
    )


def _shift_axis_right(fig: plt.Figure):
    ax2 = fig.axes[-1]
    ax2.spines["right"].set_position(("axes", 1.15))


def get_k_from_g(g: np.array, mu: float = 1.0):
    # return np.sqrt(g * 2 * mu)
    return np.sqrt(np.abs(g) * 2 * mu)


def smooth(k: list[float], w: int = 125, p: int = 2):
    """use filter to smooth over values"""
    k_hat = savgol_filter(
        k,
        w,    # window size used for filtering
        p,      # order of fitted polynomial
    )
    return k_hat


def rollavg(k: list[float], n: int = 25):
    """moving mean to smooth over values"""
    return bn.move_mean(k, window=n, min_count=None)


def butter_lowpass_filter(data, cutoff, order):
    """Butterworth lowpass filter, twice apllied"""
    n = len(data)
    nyq = n // 2

    normal_cutoff = cutoff / nyq
    # Get the filter coefficients
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data)
    return y


def sif_matczynski(traction: float, v_crk: float, c_s: float, height: float):
    """compute the analytical value for a dynamical crack in a long, narrow strip
    (Matczynski 1972)
    """
    beta = np.sqrt(1 - (v_crk / c_s) ** 2)
    n = traction * np.sqrt(beta * height / np.pi)
    return n


def sif_matczynski_34(w_0: float, v: float, height: float, mu: float = 1.0):
    """compute the analytical value for a dynamical crack in a long, narrow strip
    (Matczynski 1972)
    """
    beta = np.sqrt(1 - v**2)
    n = mu * w_0 * np.sqrt(beta / (np.pi * height))
    return n


def sif_mandal(traction: float, v_crk: float, c_s: float, height: float):
    """compute the analytical value for a dynamical crack in a long, narrow strip
    (Matczynski 1972)
    """
    beta = np.sqrt(1 - (v_crk / c_s) ** 2)
    sh = 2 * beta * height
    k = traction * np.sqrt(sh / (sh + 1))
    return k


def sif_mandal_orig(w_0: float, v_crk: float, c_s: float, height: float, mu: float = 1.0,
                    verbose: bool = True):
    """compute the analytical value for a dynamical crack in a long, narrow strip
    (Matczynski 1972)
    """
    beta = np.sqrt(1 - (v_crk / c_s) ** 2)
    k = mu * w_0 * np.sqrt(2 * beta / (height * (2 * height * beta + 1)))

    if verbose:
        print("mandal original:", k)
    return k


def least_squares(xdata: np.array, ydata: np.array, guess: float = 1, model: str = 'constant'):
    from scipy.optimize import curve_fit

    def fc(x, a):
        return a

    def fl(x, a, b):
        return a * x + b

    if model in ['l', 'linear', '1']:
        func = fl
    else:
        func = fc

    popt, pcov = curve_fit(func, xdata, ydata)

    return func(xdata, *popt), popt, pcov


def add_legend(fig: plt.Figure):
    handels, labels = [], []
    for ax in fig.axes:
        h, l = ax.get_legend_handles_labels()
        handels.extend(h)
        labels.extend(l)
    ax0 = fig.axes[0]

    ax.legend(
        handels,
        labels,
        bbox_to_anchor=(0., 1.02, 1., .102),
        loc='lower left',
        ncol=2,
        mode="expand",
        borderaxespad=0.,
    )

    return fig


def add_regression(fig: plt.Figure, name: str,
                   guess: float = 1, model: str = 'c', color: str = 'tab:purple'):
    crack_data = _read_crack_data(name)
    time = crack_data['t']
    sif = crack_data['SIF_1']

    n = len(time) // 2

    k_data = sif[-n:]
    t = time[-n:]

    k_fit, popt, pcov = least_squares(t, k_data, guess=guess, model=model)
    k_fit = k_fit * np.ones_like(t)

    ax = fig.axes[0]
    ax.plot(t, k_fit, '-', c=color, label='least squares')

    print('popt:', popt)
    return fig


def add_stat_mean(fig: plt.Figure, name: str,
                   color: str = 'tab:brown', tf: float = 15):
    """statistical evaluation of SIF: arithmetic mean & standard deviation"""

    crack_data = _read_crack_data(name)
    time = crack_data['t']
    sif = crack_data['SIF_1']

    t, n = _split_time(time, tf)
    k_data = sif[n]

    k_m = np.mean(k_data)
    sigma = np.std(k_data)

    ax = fig.axes[0]
    ax.plot(t, k_m * np.ones_like(t), linestyle='-', c=color)
    ax.plot(t, (k_m + sigma) * np.ones_like(t), linestyle=':', c=color)
    ax.plot(t, (k_m - sigma) * np.ones_like(t), linestyle=':', c=color)

    print("mean:", k_m, "std derivation:", sigma)
    return [k_m, fig]


def add_stat_median(fig: plt.Figure, name: str,
                     color: str = 'tab:olive', tf: float = 15):
    """statistical evaluation of SIF: median and 25% percentile"""

    crack_data = _read_crack_data(name)
    time = crack_data['t']
    sif = crack_data['SIF_1']

    t, n = _split_time(time, tf)
    k_data = sif[n]

    k_m = np.median(k_data)
    per25, per75 = np.percentile(k_data, [25, 75])

    ax = fig.axes[0]
    ax.plot(t, k_m * np.ones_like(t), linestyle='-', c=color)
    ax.plot(t, per75 * np.ones_like(t), linestyle=':', c=color)
    ax.plot(t, per25 * np.ones_like(t), linestyle=':', c=color)

    print("median:", k_m, "25%:", k_m - per25, "75%:", per75 - k_m)
    return [k_m, fig]


def add_lowpass(fig: plt.Figure, name: str,
                cutoff: float = 1, order: int = 1, color='tab:green'):
    """add Butterwoth lowpass filter to figure"""
    crack_data = _read_crack_data(name)
    time = crack_data['t']
    sif = crack_data['SIF_1']

    n = len(time) // 2

    k_data = sif[-n:]
    t = time[-n:]
    k = butter_lowpass_filter(k_data, cutoff, order)

    ax = fig.axes[0]
    ax.plot(t, k, linestyle='-', c=color)

    return fig


def add_k_crit(fig: plt.Figure, k_crit : float,
               color: str = 'tab:gray'):
    """add line with critical k-value to figure
    (Mandal 2020)
    """
    _line_plot(fig, k_crit)
    return fig


def add_length(fig: plt.Figure, name: str, color: str = 'tab:purple'):
    cd = _read_crack_data(name)
    a = cd['a']
    t = cd['t']

    ax = fig.axes[0]

    ax_a = ax.twinx()     # length a on right
    ax_a.set_ylabel('a')
    ax_a.plot(t, a, '-', c=color, label="crack length")

    ax_a.yaxis.label.set_color(color)
    ax_a.tick_params(axis='y', colors=color)

    if len(fig.axes) > 2:
        _shift_axis_right(fig)

    return fig


def add_velocity(fig: plt.Figure, name: str, color: str = 'tab:cyan', tip: int = 1):
    v_tip = 'v_' + str(tip)

    cd = _read_crack_data(name)
    v = cd[v_tip]
    t = cd['t']

    _, a_dot, _ = _read_model_data(name)

    ax = fig.axes[0]
    ax_v = ax.twinx()               # speed on right

    ax_v.set_ylabel('v')
    ax_v.plot(t, v, '.', markersize=4, c=color, label="crack velocity "+v_tip)

    ax_v.yaxis.label.set_color(color)
    ax_v.tick_params(axis='y', colors=color)
    cur_y, _ = ax_v.get_ylim()
    ax_v.set_ylim(cur_y, a_dot)
    ax_v.grid(axis='y')

    if len(fig.axes) > 2:
        _shift_axis_right(fig)

    return fig


def add_velocity_steps(fig: plt.Figure, name: str, color: str = 'tab:cyan', tip: int = 1):
    v_tip = 'v_' + str(tip)

    cd = _read_crack_data(name)
    v_crk = cd[v_tip]
    time = cd['t']

    ax = fig.axes[0]
    ax_v = ax.twinx()               # speed on right

    i = 1
    while i < len(time):
        t = time[i-1], time[i]
        v = v_crk[i]

        ax_v.plot(t, (v, v), '-', markersize=4, c=color, label="crack velocity "+v_tip)
        i += 1

    ax_v.set_ylabel('v')
    # ax_v.plot(t, v, '.', markersize=4, c=color, label="crack velocity "+v_tip)

    ax_v.spines['right'].set_color(color)
    ax_v.yaxis.label.set_color(color)
    ax_v.tick_params(axis='y', colors=color)

    if len(fig.axes) > 2:
        _shift_axis_right(ax_v, color)

    return fig


def add_mandal_orig(fig: plt.Figure, name: str, w_0: float, height: float,
                    mu: float = 1.0, color: str = 'tab:gray', verbose: bool = True):
    """add line with analytical value to figure - original problem
    (Mandal 2020)
    """
    _, v_crk, c_s = _read_model_data(name)

    sif = sif_mandal_orig(w_0, v_crk, c_s, height, mu=mu, verbose=verbose)
    _line_plot(fig, sif, color=color)
    return [sif, fig]


def add_mandal(fig: plt.Figure, traction: float, height: float,
                    color: str = 'tab:gray'):
    """add line with analytical value to figure - transformed problem
    (Mandal 2020)
    """
    _, v_crk, c_s = _read_model_data(name)

    sif = sif_mandal(traction, v_crk, c_s, height)
    _line_plot(fig, sif, color=color)
    return fig


def add_matczynski_34(fig: plt.Figure, name: str, w_0: float, height: float,
                     mu: float = 1.0, color: str = 'tab:gray'):
    """add line with analytical value to figure - transformed problem
    (Mandal 2020)
    """
    _, v_crk, c_s = _read_model_data(name)
    v = v_crk / c_s

    sif = sif_matczynski_34(w_0, v, height, mu=mu)
    _line_plot(fig, sif, color=color)
    return fig


def plot_tip_fields(name: str, ts: list[int],
                    n_max: int = 8, m_max: int = 8, sum_error: bool = False):
    """plot the fields (grid and replicate, point-wise error) and the crack tip + crack
    get values; build ImageGrid;
    """
    # first, get values at mesh sites
    val_grid, val_repli, val_err, err_cum, plot_data = _build_env_fields(
        name,
        ts,
        n_max=n_max,
        m_max=m_max,
        sum_error=sum_error,
    )

    # build an ImageGrid for plot output
    num_plot = len(ts)
    fig = plt.figure()
    grid = ImageGrid(
        fig, 111,                   # similar to subplot(111)
        nrows_ncols=(3, num_plot),  # creates grid of axes with rows for grid, repli, err
        axes_pad=0.1,               # pad between axes in inch
        label_mode='1',             # left and lower axes labeled
        cbar_mode='edge',           # colorbar per row
        cbar_location='right',      # where to put colorbars
    )

    # append arrays to loop over
    val_arr = np.append(val_grid, val_repli, axis=0)
    vmin, vmax = val_arr.min(), val_arr.max()       # min, max of all fields
    xx, yy = np.mgrid[-n_max:n_max + 1, -m_max:m_max + 1]
    for i, val in enumerate(val_arr):
        grid[i].pcolormesh(         # get axes from grid
            xx, yy, val,            # actual plotting of values
            cmap='jet',
            edgecolors='white',     # white edges to seperate points
            linewidth=0.5,
            shading='flat',         # no interpolation, shows values at points
            vmin=vmin,
            vmax=vmax,
        )

    emin, emax = val_err.min(), val_err.max()       # min, max of errors
    j = 2 * num_plot
    for val in val_err:
        ax = grid[j]                # get axes from grid
        ax.pcolormesh(
            xx, yy, val,            # actual plotting of values
            cmap='gray_r',
            edgecolors='white',     # white edges to seperate points
            linewidth=0.5,
            shading='flat',         # no interpolation, shows values at points
            vmin=emin,
            vmax=emax,
        )
        j += 1

    tip = FieldPoint(0, 0, 0)
    cr_dir = plot_data.cdir
    cr_dir.x *= n_max
    cr_dir.y *= m_max

    x1, y1, _ = tip.coords()
    x2, y2, _ = cr_dir.coords()

    prop_crk = '--k'
    prop_pt = 'ok'
    for i in range(num_plot):
        # crack tip
        grid[i].plot((x1, x2), (y1, y2), prop_crk)
        grid[i].plot(x1, y1, prop_pt)
        # crack as line
        grid[i + num_plot].plot((x1, x2), (y1, y2), prop_crk)
        grid[i + num_plot].plot(x1, y1, prop_pt)

    fig.set_frameon(False)

    return fig


def add_conf_k(fig: plt.Figure, name: str,
               mu: float = 1.0, norm: float = 1.0,
               color: str = 'tab:orange', tip: int = 1):
    """Add the SIF as computed from the configuration forces"""
    v_tip = 'v_' + str(tip)

    crack_data = _read_crack_data(name)
    t = crack_data['t']

    ax = fig.axes[0]
    colors = iter(('tab:orange', 'tab:red', 'tab:yellow'))

    if not hasattr(tip, '__len__'):
        tip = tip,
    for i in tip:
        g_tip = 'ConF_' + str(i)
        g = crack_data[g_tip]
        k = get_k_from_g(g) / norm
        color = next(colors)

        ax.plot(
            t, k,
            'x',
            markersize=4,
            c=color,
            label='K(G)_' + str(i),
        )
    return fig


def add_an_data(fig: plt.Figure,
                color: str = 'tab:orange'):
    t, k = _from_mat()
    ax = fig.axes[0]
    ax.plot(t.T, k.T, '-', c=color)

    return fig


def plot_sif(name: str, norm: float = 1.0,
             tip: [int, tuple[int]] = 1):
    """plot stress intensity factor (SIF) and crack length over time
    may also apply filter to SIF for smoothing
    parameters: name <- string, for file names
                plot_a <- boolean, plot length over time or not
    """

    crack_data = _read_crack_data(name)
    t = crack_data['t']

    fig, ax = plt.subplots()
    colors = iter(('tab:blue', 'tab:green', 'tab:olive'))

    if not hasattr(tip, '__len__'):
        tip = tip,
    for i in tip:
        k_tip = 'SIF_' + str(i)
        k = crack_data[k_tip] / norm
        color = next(colors)

        ax.plot(
            t, k,
            '+',
            markersize=4,
            c=color,
            label=k_tip,
        )

    ax.set_xlabel('t')
    if norm == 1.0:
        ax.set_ylabel('K_III')
    else:
        ax.set_ylabel('K_III / K_c')
    fig.set_tight_layout(True)

    return fig


def plot_k_max(name: str,
                  mu: float = 1.0, t0: float = 2.0, h: int = 256):

    sif, k_an = _from_local(name, mu=mu, t0=t0)

    time = sif.pop('time')
    k_max = np.array([max(k_t) for k_t in sif.values()]) / np.sqrt(2)
    d = np.array(list(sif.keys()))

    fig, ax = plt.subplots()

    ax.plot(d, k_max / k_an, '.-', c='tab:blue', markersize=4, label="SIF h=1/{}".format(h))
    ax.plot(d, np.ones_like(d), '--', c='grey', markersize=4, label="SIF analytical")

    ax.set_xlabel('d')
    ax.set_ylabel('K_III')

    ax.grid()
    fig.set_tight_layout(True)

    return fig


def plot_k_err(name: str, comp: float,
               tf: float = 15, mu: float = 1.0, h: int = 256):

    sif, _ = _from_local(name, mu=mu)

    time = sif.pop('time')
    t, n = _split_time(time, tf)

    k_med = np.array([np.median(k_t[n]) for k_t in sif.values()])
    d = np.array(list(sif.keys()))

    err1 = 1 - k_med / comp

    fig, ax_k = plt.subplots()
    ax_e = ax_k.twinx()

    ax_k.plot(d, k_med / comp, 'o-', c='tab:blue', markersize=4, label="SIF h=1/{}".format(h))
    ax_e.plot(d, err1, '+-', c='tab:red', markersize=4, label="absolute")

    ax_k.set_xlabel('d')
    ax_k.set_ylabel('K_III')
    ax_e.set_ylabel('deviation')

    ax_e.set_ylim((-0.5, 0.5))
    ax_e.grid()
    fig.set_tight_layout(True)

    return fig


def plot_k_at_d(name: str, d: float,
                  mu: float = 1.0, t0: float = 2.0, h: int = 256, comp: float = None):

    sif, k_an = _from_local(name, mu=mu, t0=t0)

    time = sif.pop('time')

    if not d in sif:
        dist = [r for r in sif if r >= d]
        d = dist[0]

    k = sif[d]

    print(d)

    if comp:
        k_an = comp

    fig, ax = plt.subplots()

    ax.plot(time, k, '.-', c='tab:blue', markersize=4, label="SIF h=1/{}".format(h))
    ax.plot(time, np.ones_like(time) * k_an, '--', c='grey', markersize=4, label="SIF analytical")

    ax.set_xlabel('d')
    ax.set_ylabel('K_III')

    ax.grid()
    fig.set_tight_layout(True)

    return fig


def plot_k_at_d_inf(name: str, d: float,
                  mu: float = 1.0, t0: float = 2.0, h: int = 256, comp: float = None):

    sif, k_an = _from_local(name, mu=mu, t0=t0)

    time = sif.pop('time')

    if not d in sif:
        dist = [r for r in sif if r >= d]
        d = dist[0]

    k = sif[d] * np.sqrt(2)

    print(d)

    if comp:
        k_an = comp

    fig, ax = plt.subplots()

    ax.plot(time / 2, k, '.-', c='tab:blue', markersize=4, label="SIF h=1/{}".format(h))
    ax.plot(time / 2, np.ones_like(time), '--', c='grey', markersize=4, label="SIF analytical")

    ax.set_xlabel('d')
    ax.set_ylabel('K_III')

    ax.grid()
    fig.set_tight_layout(True)

    return fig


def save_plot(fig: plt.Figure, path: str, tex: bool = False, dpi: int = 300):
    """Saves a plot to disk
    The default is a .png file.
    If tex is set to True a .pgf is generated. This requires XeLatex to be installed.
    """
    if not tex:
        path += '.png'
    else:               # save as .pgf for LaTeX documents
        mpl.use("pgf")
        mpl.rcParams.update({
            # 'pgf.texsystem': 'pdflatex',  # not working with pdflatex due to bug in matplotlib
            'pgf.texsystem': 'xelatex',
            'font.family': 'serif',
            'font.size': 10,
            'text.usetex': True,
            'pgf.rcfonts': False,
        })
        path += '.pgf'
    fig.savefig(path, dpi=dpi)


def set_size(fig: plt.Figure, name:str, size: tuple[float] = (12, 8)):
    a, b = size     # in cm
    a *= 0.3937     # to inch
    b *= 0.3937

    crack_data = _read_crack_data(name)
    t = crack_data['t']
    t0 = t[0]
    tf = t[-1]

    fig.set_size_inches((a, b))
    ax = fig.axes[0]
    ax.set_xlim((t0, tf))

    return fig
