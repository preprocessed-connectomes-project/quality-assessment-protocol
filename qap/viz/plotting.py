#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
import math
import os.path as op
import numpy as np
import nibabel as nb
import pandas as pd
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.backends.backend_pdf import FigureCanvasPdf as FigureCanvas
import seaborn as sns


def calculate_gray_plot(func_file, mask_file):
    import numpy as np
    from matplotlib import pyplot as plt
    import nibabel as nb
    from sklearn import preprocessing

    # read in and mask the data 
    func = nb.load(func_file)
    mask = nb.load(mask_file)
    mask_affine = mask.affine
    mask = mask.get_data()

    new_shape = func.shape[0]*func.shape[1]*func.shape[2]
    func = func.get_data()
    func = np.resize(func, (new_shape,func.shape[-1] ))
    func = func[np.where( mask.ravel() == 1 )[0], :]

    # z-score the rows so that they can be viewed on the same image
    func = preprocessing.scale(func, axis=1)

    #now try some basic reordering
    from sklearn.cluster.bicluster import SpectralBiclustering
    model = SpectralBiclustering(n_clusters=(5,1),n_components=5, method='log',random_state=0)
    model.fit(func)
    reordered_func=func[np.argsort(model.row_labels_)]

    #create a matrix with cluster ids, 
    #make sure cluster ids > func values so we can display both easily. 
    out_clusters=np.zeros(np.shape(func[:,0]))
    func_max = func.max()+1
    for i in range(0,5):
        out_clusters[model.get_indices(i)[0]]=int(i*2+func_max)
    reordered_clusters=out_clusters[np.argsort(model.row_labels_)]

    #calculte each cluster global signal
    cluster_gs = []
    for i in np.unique(out_clusters):
        c = func[np.where( out_clusters.ravel() == i )[0], :]
        time = func.shape[-1]
        output = [0]*time
    
        for i in range(time):
            output[i] = c[:,i].mean()
        output = (output - min(output))/(max(output) - min(output))
        cluster_gs.append(output)


    return reordered_func, reordered_clusters, cluster_gs

def plot_measures(df, measures, ncols=4, title='Group level report',
                  subject=None, figsize=(8.27, 11.69)):
    import matplotlib.gridspec as gridspec
    print "plot measures"
    nmeasures = len(measures)
    nrows = nmeasures // ncols
    if nmeasures % ncols > 0:
        nrows += 1

    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(nrows, ncols)

    axes = []

    for i, mname in enumerate(measures):
        axes.append(plt.subplot(gs[i]))
        axes[-1].set_xlabel(mname)
        sns.distplot(
            df[[mname]], ax=axes[-1], color="b", rug=True,  norm_hist=True)

        # labels = np.array(axes[-1].get_xticklabels())
        # labels[2:-2] = ''
        axes[-1].set_xticklabels([])
        plt.ticklabel_format(style='sci', axis='y', scilimits=(-1, 1))

        if subject is not None:
            subid = subject
            subdf = df.loc[df['Participant'] == subid]
            sessions = np.atleast_1d(subdf[['Session']]).reshape(-1).tolist()

            for ss in sessions:
                sesdf = subdf.loc[subdf['Session'] == ss]
                scans = np.atleast_1d(sesdf[['Series']]).reshape(-1).tolist()

                for sc in scans:
                    scndf = subdf.loc[sesdf['Series'] == sc]
                    plot_vline(
                        scndf.iloc[0][mname], '%s_%s' % (ss, sc), axes[-1])

    fig.suptitle(title)
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.subplots_adjust(top=0.85)
    plt.close()
    return fig


def plot_all(df, groups, subject=None, figsize=(11.69, 5),
             strip_nsubj=10, title='Summary report'):
    import matplotlib.gridspec as gridspec
    print "plot_all"
    colnames = [v for gnames in groups for v in gnames]
    lengs = [len(el) for el in groups]
    ncols = np.sum(lengs)

    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(1, len(groups), width_ratios=lengs)

    subjects = sorted(pd.unique(df.Participant.ravel()))
    nsubj = len(subjects)

    if subject:
        if subject not in subjects:
            return None

    subid = subject

    df["Participant"] = df["Participant"].astype(str)

    axes = []
    for i, snames in enumerate(groups):
        axes.append(plt.subplot(gs[i]))

        if nsubj > strip_nsubj:
            sns.violinplot(data=df[snames], ax=axes[-1])
        else:
            stdf = df.copy()
            if subid is not None:
                stdf = stdf.loc[stdf['Participant'] != subid]
            try:
                sns.stripplot(data=stdf[snames], ax=axes[-1], jitter=0.25)
            except KeyError:
                # handle the possibility of one or multiple phase-encoding
                # directions for GSR measure
                if "Ghost_" in snames[0]:
                    for sname in snames:
                        try:
                            sns.stripplot(data=stdf[[sname]], ax=axes[-1],
                                          jitter=0.25)
                        except KeyError:
                            pass

        axes[-1].set_xticklabels(
            [el.get_text() for el in axes[-1].get_xticklabels()],
            rotation='vertical')
        plt.ticklabel_format(style='sci', axis='y', scilimits=(-1, 1))
        # df[snames].plot(kind='box', ax=axes[-1])

        # If we know the subject, place a star for each scan
        if subject is not None:
            subdf = df.loc[df['Participant'] == str(subid)]
            scans = sorted(pd.unique(subdf.Series.ravel()))
            nstars = len(scans)
            for j, s in enumerate(snames):
                vals = []
                for k, scid in enumerate(scans):
                    val = subdf.loc[df.Series == scid, [s]].iloc[0, 0]
                    vals.append(val)
                if len(vals) != nstars:
                    continue

                pos = [j]
                if nstars > 1:
                    pos = np.linspace(j-0.3, j+0.3, num=nstars)
                axes[-1].plot(
                    pos, vals, ms=9, mew=.8, linestyle='None',
                    color='w', marker='*', markeredgecolor='k',
                    zorder=10)

    fig.suptitle(title)
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.subplots_adjust(top=0.85)
    plt.close()
    return fig


def plot_mosaic(nifti_file, title=None, overlay_mask=None,
                figsize=(11.7, 8.3)):
    from six import string_types
    from pylab import cm
    print "plot_mosaic"
    if isinstance(nifti_file, string_types):
        nii = nb.load(nifti_file)
        mean_data = nii.get_data()
    else:
        mean_data = nifti_file

    z_vals = np.array(range(0, mean_data.shape[2]))
    # Reduce the number of slices shown
    if mean_data.shape[2] > 70:
        rem = 15
        # Crop inferior and posterior
        mean_data = mean_data[..., rem:-rem]
        z_vals = z_vals[rem:-rem]
        # Discard one every two slices
        mean_data = mean_data[..., ::2]
        z_vals = z_vals[::2]

    n_images = mean_data.shape[2]
    row, col = _calc_rows_columns(figsize[0] / figsize[1], n_images)

    if overlay_mask:
        overlay_data = nb.load(overlay_mask).get_data()

    # create figures
    fig = plt.Figure(figsize=figsize)
    FigureCanvas(fig)

    fig.subplots_adjust(top=0.85)
    for image, z_val in enumerate(z_vals):
        ax = fig.add_subplot(row, col, image + 1)
        data_mask = np.logical_not(np.isnan(mean_data))
        if overlay_mask:
            ax.set_rasterized(True)

        ax.imshow(np.fliplr(mean_data[:, :, image].T), vmin=np.percentile(
            mean_data[data_mask], 0.5),
            vmax=np.percentile(mean_data[data_mask], 99.5),
            cmap=cm.Greys_r, interpolation='nearest', origin='lower')

        if overlay_mask:
            cmap = cm.Reds  # @UndefinedVariable
            cmap._init()
            alphas = np.linspace(0, 0.75, cmap.N + 3)
            cmap._lut[:, -1] = alphas
            ax.imshow(np.fliplr(overlay_data[:, :, image].T), vmin=0, vmax=1,
                      cmap=cmap, interpolation='nearest', origin='lower')

        ax.annotate(
            str(z_val), xy=(.95, .015), xycoords='axes fraction',
            fontsize=10, color='white', horizontalalignment='right',
            verticalalignment='bottom')

        ax.axis('off')

    fig.subplots_adjust(
        left=0.05, right=0.95, bottom=0.05, top=0.95, wspace=0.01, hspace=0.1)

    if not title:
        _, title = op.split(nifti_file)
        title += " (last modified: %s)" % time.ctime(
            op.getmtime(nifti_file))
    fig.suptitle(title, fontsize='10')
    return fig


def plot_fd(meanfd_file, dvars, global_signal, metadata, figsize=(11.7, 8.3), mean_fd_dist=None, title='Mean FD, DVARS, Global Signal'):
    fd_power = _calc_fd(meanfd_file)
    global_signal = (global_signal - min(global_signal))/(max(global_signal) - min(global_signal))
    x = metadata[0]
    a = dvars[x]
    dvars = a['Standardized DVARS']
    fig = plt.Figure(figsize=figsize)
    FigureCanvas(fig)

    if mean_fd_dist:
        grid = GridSpec(2, 4)
    else:
        grid = GridSpec(1, 2, width_ratios=[3, 1])
        grid.update(hspace=1.0, right=0.95, left=0.1, bottom=0.2)

    ax = fig.add_subplot(grid[0, :-1])
    fd = ax.plot(fd_power, label='Mean FD')
    d, = ax.plot(dvars, label='DVARS')
    gs, = ax.plot(global_signal, label='Global Signal')
    ax.set_xlim((0, len(fd_power)))
    ax.set_ylabel("Frame Displacement [mm], DVARS and Global Signal")
    ax.set_xlabel("Frame number")
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels)
    ylim = ax.get_ylim()

    ax1 = fig.add_subplot(grid[0, -1])
    sns.distplot(fd_power, vertical=True, ax=ax1)
    ax1.set_ylim((min(fd_power), max(fd_power)))

    if mean_fd_dist:
        ax = fig.add_subplot(grid[1, :])
        sns.distplot(mean_fd_dist, ax=ax)
        ax.set_xlabel("Mean Frame Displacement (over all subjects) [mm]")
        mean_fd = fd_power.mean()
        label = r'$\overline{\text{FD}}$ = %g' % mean_fd
        plot_vline(mean_fd, label, ax=ax)

    fig.suptitle(title)
    return fig


def plot_dist(
        main_file, mask_file, xlabel, distribution=None, xlabel2=None,
        figsize=(11.7, 8.3)):
    print "plot_dist"
    data = _get_values_inside_a_mask(main_file, mask_file)

    fig = plt.Figure(figsize=figsize)
    FigureCanvas(fig)

    gs = GridSpec(2, 1)
    ax = fig.add_subplot(gs[0, 0])
    sns.distplot(data.astype(np.double), kde=False, bins=100, ax=ax)
    ax.set_xlabel(xlabel)

    ax = fig.add_subplot(gs[1, 0])
    sns.distplot(np.array(distribution).astype(np.double), ax=ax)
    cur_val = np.median(data)
    label = "%g" % cur_val
    plot_vline(cur_val, label, ax=ax)
    ax.set_xlabel(xlabel2)

    return fig


def plot_vline(cur_val, label, ax):
    print "plot_vline"
    ax.axvline(cur_val)
    ylim = ax.get_ylim()
    vloc = (ylim[0] + ylim[1]) / 2.0
    xlim = ax.get_xlim()
    pad = (xlim[0] + xlim[1]) / 100.0
    ax.text(cur_val - pad, vloc, label, color="blue", rotation=90,
            verticalalignment='center', horizontalalignment='right')


def _calc_rows_columns(ratio, n_images):
    print "_calc_rows_columns"
    rows = 1
    for _ in range(100):
        columns = math.floor(ratio * rows)
        total = rows * columns
        if total > n_images:
            break

        columns = math.ceil(ratio * rows)
        total = rows * columns
        if total > n_images:
            break
        rows += 1
    return rows, columns


def _calc_fd(fd_file):
    """Calculate the frame-wise displacement (FD) given the FD vector.

    :type fd_file: str
    :param fd_file: The filepath to the frame-wise displacement 1D vector
    file.
    :rtype: NumPy array
    :return: The frame-wise displacement (FD) array.
    """
    print "_calc_fd"
    lines = open(fd_file, 'r').readlines()
    rows = [[float(x) for x in line.split()] for line in lines]
    cols = np.array([list(col) for col in zip(*rows)])

    translations = np.transpose(np.abs(np.diff(cols[0:3, :])))
    rotations = np.transpose(np.abs(np.diff(cols[3:6, :])))

    FD_power = np.sum(translations, axis=1) + \
        (50 * 3.141 / 180) * np.sum(rotations, axis=1)

    # FD is zero for the first time point
    FD_power = np.insert(FD_power, 0, 0)

    return FD_power


def _get_mean_fd_distribution(fd_files):
    print "_get_mean_fd_distribution"
    mean_FDs = []
    max_FDs = []
    for fd_file in fd_files:
        FD_power = _calc_fd(fd_file)
        mean_FDs.append(FD_power.mean())
        max_FDs.append(FD_power.max())

    return mean_FDs, max_FDs


def _get_values_inside_a_mask(main_file, mask_file):
    print "_get_values_inside_a_mask"
    main_nii = nb.load(main_file)
    main_data = main_nii.get_data()
    nan_mask = np.logical_not(np.isnan(main_data))
    mask = nb.load(mask_file).get_data() > 0

    data = main_data[np.logical_and(nan_mask, mask)]
    return data

def grayplot(func_file, mask_file, meanfd_file, dvars, global_signal, metadata, figsize=(11.7, 8.3), title='Timeseries Plot'):
    from sklearn import preprocessing
    fd_power = _calc_fd(meanfd_file)
    gray_matrix, color_matrix, cluster_gs = calculate_gray_plot(func_file, mask_file)

    #zscore data
    fd_power = preprocessing.scale(fd_power)
    dvars = preprocessing.scale(dvars)
    global_signal = preprocessing.scale(global_signal)


    #create grid with 2 rows
    figsize=(11, 8)
    fig = plt.Figure(figsize=figsize)
    #fig,ax1 = plt.subplots()
    FigureCanvas(fig)

    grid = GridSpec(2, 1, wspace=0.01)
    ax1 = fig.add_subplot(grid[0, 0])
    ax2 = fig.add_subplot(grid[1, 0])

    #gray plot part
    #add brain labels to grayplot
    gray_plot = np.insert(gray_matrix, gray_matrix.shape[1], color_matrix, axis=1)
    gray_plot = np.insert(gray_plot, gray_matrix.shape[1], color_matrix, axis=1)

    #plot grayplot
    from numpy.ma import masked_array
    func_max = gray_matrix.max()+1
    gray = masked_array(gray_plot,gray_plot>=func_max)
    colors = masked_array(gray_plot,gray_plot<func_max)

    aspect = (float(gray_matrix.shape[1])/float(gray_matrix.shape[0]))/3.0
    pa = ax1.imshow(gray,interpolation='None',cmap='gray',aspect=aspect)
    pb = ax1.imshow(colors,interpolation='None',cmap='rainbow',aspect=aspect)
    ax1.set_xlabel('Frame number')
    ax1.set_ylabel('Voxels')
    ax1.yaxis.set_ticklabels([])
    ax1.grid(False)

    #dvars plot
    fd = ax2.plot(fd_power, label='Mean FD', color='#7bc71b')
    d, = ax2.plot(dvars, label='DVARS',color='#FFC830')
    gs, = ax2.plot(global_signal, label='Global Signal',color='#FF5624')
    cluster_colors = ['#ffffff','#7f00ff', '#2adcdc', '#d4dc7f', '#ff0000']
    for i in range(5):
        c, = ax2.plot(cluster_gs[i], color=cluster_colors[i], linewidth=0.7, label='cluster gs')

    #adjust fig params
    ax2.set_xlim((0, len(fd_power)))
    ax2.set_xlabel("Frame number")
    ylim = ax2.get_ylim()
    fig.suptitle(title)
    ax2.grid(True)

    #remove cluster gs labls
    handles, labels = ax2.get_legend_handles_labels()
    newLabels, newHandles = [], []
    for handle, label in zip(handles, labels):
      if label != 'cluster gs':
        newLabels.append(label)
        newHandles.append(handle)
    ax2.legend(newHandles, newLabels) 
    
    plt.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0)

    #create nii file with cluster values
    mask = nb.load(mask_file)
    mask_affine = mask.affine
    mask = mask.get_data()
    out_clusters_img = np.zeros(np.prod(mask.shape)) 
    out_clusters_img[mask.ravel()==1] = cluster_gs 
    out_clusters_img = out_clusters_img.reshape(mask.shape)
    clust_img = nb.Nifti1Image(out_clusters_img, mask_affine)
    
    #grid.update(hspace=1.0, right=0.95, left=0.1, bottom=0.2)
    #ax = fig.add_subplot(grid[0, :-1])

    return fig, clust_img