# -*- coding: utf-8 -*-
# coding=utf8

import matplotlib
# note: some matplotlib backends are broken and cause segfaults in fig.save.
# try commenting out and in the use of Cairo in case of problems ...
#try:
matplotlib.use('Cairo')
#except ImportError: pass

import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import matplotlib.text
import matplotlib.lines
import matplotlib.patches

from scipy import interpolate
import numpy as np

import math

def add_xlabel(axes, text, *args, **kwargs):
    label = axes.set_xlabel(text, size='large', ha='right', *args, **kwargs)
    label.set_position((1.0, 0.03))
    return label

def add_ylabel(axes, text, *args, **kwargs):
    label = axes.set_ylabel(text, size='large', va='top', *args, **kwargs)
    label.set_position((-0.03, 1.0))
    return label

# plotdata represents the data of a single curve in a plot, including drawing options, legend, etc.
class plotdata:
    def __init__(self, color = '#000000', legend = None, as_function = False, lw = 2, legend_order = 0):
        self.x = []
        self.y = []
        self.xmax = None # required only for one-binned histograms ...
        self.legend = legend
        self.legend_order = legend_order
        self.yerrors = None
        # capsize for the error bars:
        self.capsize = 1.5
        # how to draw the yerrors; valid are "bars" for error bars (using color, capsize and lw), 'bars0' to also show the marker for y values of 0,
        # or 'area' to draw only a shaded area.
        self.yerrors_mode = 'bars'
        self.yerrors_fill_alpha = 0.5
        self.xerrors = None
        # filling will be done only is fill_color is not None
        self.fill_color = None
        # in case filling is done, the x value range to fill. The default "None"
        # is to fill everything (i.e., whole x-range). Otherwise, specify tuple (xmin, xmax) here.
        self.fill_xrange = None
        self.fill_to_y = 0.0
        self.color = color
        self.marker = 'None'
        self.markersize = 1.0
        self.lw = lw
        self.fmt = '-'
        # an array of bands; a band is a three-tuple (y1, y2, color). y1 and y2 are 
        # arrays of y values.
        # bands to draw first should come first 
        self.bands = None
        self.band_lw = 0
        self.bands_fill = True
        self.as_function = as_function
        self.draw_histo = True
        self.draw_line = True
        
    # make a histogram of the given values
    def histogram(self, values, xmin, xmax, nbins, errors = False):
        xmin, xmax = float(xmin), float(xmax)
        self.xmax = xmax
        self.x = [xmin + (xmax - xmin) / nbins * i for i in range(nbins)]
        self.y = [0.0] * nbins
        if errors: self.yerrors = [0.0] * nbins
        for v in values:
            ibin = int((v - xmin) / (xmax - xmin) * nbins)
            if ibin < 0 or ibin >= nbins: continue
            self.y[ibin] += 1
            if errors: self.yerrors[ibin] += 1
        if errors: self.yerrors = map(math.sqrt, self.yerrors)

    # scale all y values by factor
    def scale_y(self, factor):
        self.y = [y*factor for y in self.y]
        if self.bands is None: return
        for band in self.bands:
            band[0][:] = [y * factor for y in band[0]]
            band[1][:] = [y * factor for y in band[1]]

    def scale_x(self, factor):
        self.x = [x*factor for x in self.x]
        self.xmax *= factor

    # set data according to the "histo triple" h = (xmin, xmax, data) or Histogram instance
    def histo_triple(self, h):
        binwidth = (h[1] - h[0]) / len(h[2])
        self.x = [h[0] + i * binwidth for i in range(len(h[2]))]
        self.y = h[2][:]
        self.xmax = h[1]
        
    # set values from Histogram object:
    def set_histogram(self, histo):
        self.xmax = histo.get_xmax()
        self.x = [histo.get_x_low(i) for i in range(histo.get_nbins())]
        self.y = histo.get_values()
        self.yerrors = histo.get_uncertainties()
        
    
    # replace x, y and bands by a smoothed version, obtained by cubic interpolation
    # evaluated n times more points than original
    #
    # s is the amount of smoothing (default None is the scipy default)
    #
    # relunc is the relative uncertainty for each y value to assume for s
    #
    # in case more control is needed, make the smoothing outside and re-set the x,y values.
    #
    # should only be used in cases self.as_function = True
    def smooth(self, n = 3, s = None, relunc = 0.05, miny_factor = 0.0):
        oldx = self.x[:]
        # assume a 5% uncertainty for the smoothing
        y_average = sum(self.y) / len(self.y)
        tck = interpolate.splrep(oldx, self.y, w = [1 / (relunc * max(self.y[i], y_average * miny_factor)) for i in range(len(oldx))], s = s)
        self.x = list(np.linspace(min(self.x), max(self.x), n * len(self.x)))
        self.y = interpolate.splev(self.x, tck)
        if self.bands is None: return
        for band in self.bands:
            tck = interpolate.splrep(oldx, band[0], w = [1 / (relunc * band[0][i]) for i in range(len(oldx))], s = s)
            band[0][:] = interpolate.splev(self.x, tck)
            tck = interpolate.splrep(oldx, band[1], w = [1 / (relunc * band[1][i]) for i in range(len(oldx))], s = s)
            band[1][:] = interpolate.splev(self.x, tck)
        
    # ofile is a string (filename) or a handle to an open file
    def write_txt(self, ofile):
        if type(ofile)==str: ofile = open(ofile, 'w')
        ofile.write('# x; y')
        if self.yerrors is not None: ofile.write('; yerror')
        if self.bands is not None:
            for k in range(len(self.bands)):
                ofile.write('; band %d low; band %d high' % (k, k))
        ofile.write("\n")
        for i in range(len(self.x)):
            ofile.write("%10.5g %10.5g " % (self.x[i], self.y[i]))
            if self.yerrors is not None:
                ofile.write("%10.5g " % self.yerrors[i])
            if self.bands is not None:
                for k in range(len(self.bands)):
                    ofile.write("%10.5g %10.5g" % (self.bands[k][0][i], self.bands[k][1][i]))
            ofile.write("\n")
    
    # infile is the filename
    def read_txt(self, infile):
        # values is a list of lines in the file; each line is a list of floats
        values = []
        for line in file(infile):
            if len(line) == 0: continue
            if line[0] == '#': continue
            line_values = map(lambda s: float(s), line.split())
            # check that the number of values in the lines agree:
            if len(values) > 0:
                if len(values[0]) != len(line_values): raise RuntimeError, "number of values given is inconsistent!"
            values.append(line_values)
        n_values = len(values[0])
        have_yerrors = (n_values % 2 == 1)
        # read x, y values:
        self.x = [row[0] for row in values]
        self.y = [row[1] for row in values]
        if have_yerrors:
            self.yerrors = [row[2] for row in values]
        # read bands:
        n_bands = (n_values - 2) / 2
        self.bands = []
        colors = ['#ffff00', '#00ff00']
        yerror_offset = 0
        if have_yerrors: yerror_offset = 1
        for i in range(n_bands):
            band = ([row[2+2*i +  yerror_offset] for row in values], [row[3+2*i  + yerror_offset] for row in values], colors[i % len(colors)])
            self.bands.append(band)
        
        

## \brief Make a plot and write it to an output file
#
#
# histos is a list / tuple of plotdata instances, xlabel and ylabel are lables for the axes, outname is the output filename.
#
# logx and logy control whether the x/y-scale should be logarithmic.
#
# If set, ax_modifier should be a function / callable. It will be called with the matplotlib Axes instance as only argument.
# This can be used for arbitrary mainulation, such as adding other objects to the plot, etc.
#
# title_ul and title_ur are strings shown on the upper left and upper right of the plot, resp.
#
# extra_legend_items is a list of extra items to add to the legend, as tuple (handle, lable) where handle is a matplotlib Artist
# and label the legend string. As special case, we also allow handle to be a string in which case it is assumed to encode a color.
# For other cases, see matplotlib legend guide for details.
#
# xmin, xmax, ymin, ymax control the region shown in the plot. the default is to determine the region automatically (by matplotblib)
def plot(histos, xlabel, ylabel, outname = None, logy = False, logx = False, ax_modifier=None, title_ul=None, title_ur = None,
 extra_legend_items = [], xmin = None, xmax=None, ymin=None, ymax=None, legend_args = {}, fig = None, figsize_cm = (15, 12), fontsize = 10, axes_creation_args = {}):
    #extra_legend_items = extra_legend_items[:]
    legend_items = []
    cm = 1.0/2.54
    fsize = figsize_cm[0]*cm, figsize_cm[1]*cm
    fp = fm.FontProperties(size = fontsize)
    if fig is None:
        fig = plt.figure(figsize = fsize)
    rect = axes_creation_args.pop('rect', (0.15, 0.15, 0.8, 0.75))
    ax = fig.add_axes(rect, **axes_creation_args)
    if logy: ax.set_yscale('log')
    if logx: ax.set_xscale('log')
    add_xlabel(ax, xlabel, fontproperties=fp)
    add_ylabel(ax, ylabel, fontproperties=fp)
    if title_ul is not None:
        if type(title_ul) == type([]):
            yoffset = 1.02
            for s in title_ul:
                ax.text(0.0 if yoffset>1 else 0.02, yoffset, s, transform = ax.transAxes, ha='left', va = 'bottom' if yoffset > 1 else 'top')
                yoffset -= fontsize * 1.0 / (72 * fsize[1]) * 1.5
        else:
            ax.text(0.0, 1.02, title_ul, transform = ax.transAxes, ha='left', va='bottom')
    if title_ur is not None: ax.text(1.0, 1.02, title_ur, transform = ax.transAxes, ha='right', va='bottom')
    draw_legend = False
    for histo in histos:
        legend_added = False
        assert len(histo.x)==len(histo.y), "number of x,y coordinates not the same for '%s'" % histo.legend
        if histo.legend: draw_legend = True
        # allow empty "dummy" plots which have legend but no content:
        if len(histo.x)==0:
            if histo.fill_color is not None:
                legend_items.append((histo.legend_order, matplotlib.patches.Rectangle((0, 0), 1, 1, fc=histo.fill_color, ec=histo.color, lw=histo.lw), histo.legend))
            else:
                legend_items.append((histo.legend_order, matplotlib.lines.Line2D((0, 1, 2), (0, 0, 0), color=histo.color, lw=histo.lw), histo.legend))
            continue
        if histo.bands is not None:
            for band in histo.bands:
                if histo.bands_fill:
                    ax.fill_between(histo.x, band[0], band[1], lw=histo.band_lw, facecolor=band[2], color=band[2])
                else:
                    xs = histo.x + [x for x in reversed(histo.x)]
                    ys = band[0] + [y for y in reversed(band[1])]
                    xs.append(xs[0])
                    ys.append(ys[0])
                    ax.plot(xs, ys, lw=histo.band_lw, color=band[2])
        if not histo.as_function:
            if histo.yerrors is None or histo.draw_histo:
               new_x = [histo.x[0]]
               for x in histo.x[1:]: new_x += [x]*2
               new_x += [histo.xmax]
               new_y = []
               for y in histo.y: new_y += [y]*2
               if logy and ymin is not None:
                    for i in range(len(new_y)): new_y[i] = max(new_y[i], ymin)
               if histo.fill_color is not None:
                    if histo.fill_xrange is None:
                        ax.fill_between(new_x, new_y, [0] * len(new_y), lw=histo.lw, color=histo.color, facecolor = histo.fill_color)
                        legend_items.append((histo.legend_order, matplotlib.patches.Rectangle((0, 0), 1, 1, fc=histo.fill_color, ec=histo.color, lw=histo.lw), histo.legend))
                    else:
                        x_clipped = [histo.fill_xrange[0]]
                        y_clipped = []
                        for x,y in zip(new_x, new_y):
                            if x >= histo.fill_xrange[0]:
                                if len(y_clipped)==0: y_clipped.append(y)
                                if x <= histo.fill_xrange[1]:
                                    x_clipped.append(x)
                                    y_clipped.append(y)
                                else:
                                    x_clipped.append(histo.fill_xrange[1])
                                    y_clipped.append(y_clipped[-1])
                                    break
                        ax.fill_between(x_clipped, y_clipped, [histo.fill_to_y] * len(y_clipped), lw=0, color=None, facecolor = histo.fill_color)
                        ax.plot(new_x, new_y, histo.fmt, lw=histo.lw, color=histo.color)
                        legend_items.append((histo.legend_order, matplotlib.patches.Rectangle((0, 0), 1, 1, fc=histo.fill_color, ec=histo.color, lw=histo.lw), histo.legend))
                        #legend_items.append((histo.legend_order, matplotlib.lines.Line2D((0, 1, 2), (0, 0, 0), color=histo.color, lw=histo.lw), histo.legend))
               else:
                   ax.plot(new_x, new_y, histo.fmt, lw=histo.lw, color=histo.color)
                   legend_items.append((histo.legend_order, matplotlib.lines.Line2D((0, 1, 2), (0, 0, 0), color=histo.color, lw=histo.lw, ls = histo.fmt), histo.legend))
               legend_added = True
            # if histo.yerrors is set, draw with errorbars, shifted by 1/2 binwidth ...
            if histo.yerrors is not None:
                if histo.yerrors_mode.startswith('bars'):
                    low_x = histo.x + [histo.xmax]
                    x_centers = [0.5 * (low_x[i] + low_x[i+1]) for i in range(len(histo.x))]
                    ys = histo.y
                    yerrors = histo.yerrors
                    if histo.yerrors_mode != 'bars0':
                        x_new, y_new, ye_new = [], [], []
                        for x,y,ye in zip(x_centers, ys, yerrors):
                            if y != 0:
                                x_new.append(x)
                                y_new.append(y)
                                ye_new.append(ye)
                        x_centers, ys, yerrors = x_new, y_new, ye_new
                    ax.plot(x_centers, ys, marker = histo.marker, markersize=histo.markersize, ls='None', mew = 0.0, mfc = histo.color)
                    ax.errorbar(x_centers, ys, yerrors, ecolor = histo.color, capsize = histo.capsize, lw = histo.lw, fmt = None)
                    if not legend_added:
                        legend_items.append((histo.legend_order, matplotlib.lines.Line2D((0, 1, 2), (0, 0, 0), color=histo.color, marker=histo.marker, markevery=(1,10), mew=0.0, markersize=histo.markersize, lw=histo.lw), histo.legend))
                        legend_added = True
                else:
                    new_x = [histo.x[0]]
                    for x in histo.x[1:]: new_x += [x]*2
                    new_x += [histo.xmax]
                    new_y_low, new_y_high = [], []
                    for y, ye in zip(histo.y, histo.yerrors):
                        new_y_low += [y - ye]*2
                        new_y_high += [y + ye]*2
                    ax.fill_between(new_x, new_y_high, new_y_low, lw=histo.lw, color = histo.color, facecolor = histo.fill_color, alpha = histo.yerrors_fill_alpha)
                    if not legend_added:
                        legend_items.append((histo.legend_order, matplotlib.patches.Rectangle((0, 0), 1, 1, fc=histo.fill_color, ec=histo.color, lw=histo.lw, alpha = histo.yerrors_fill_alpha), histo.legend))
                        legend_added = True
        else:
            if histo.yerrors is not None:
                if histo.yerrors_mode == 'bars':
                    lw = histo.lw
                    if histo.draw_line is False: lw = 0
                    ax.errorbar(histo.x, histo.y, histo.yerrors, elinewidth = histo.lw, lw=lw, label=histo.legend, color=histo.color, marker=histo.marker, markersize=histo.markersize)                    
                else:
                    raise RuntimeError, "yerrors_mode!='bars' for as_function=True not supported!"
            else:
                if histo.fill_color is not None:
                    ax.fill_between(histo.x, histo.y, [histo.fill_to_y] * len(histo.y), lw=histo.lw, label=histo.legend, color=histo.color, facecolor = histo.fill_color)
                else:
                    ax.plot(histo.x, histo.y, histo.fmt, lw=histo.lw, color=histo.color, marker=histo.marker, markersize=histo.markersize)
                    legend_items.append((histo.legend_order, matplotlib.lines.Line2D((0, 1, 2), (0, 0, 0), color=histo.color, ls = histo.fmt, lw=histo.lw), histo.legend))
    legend_items = [(h,l) for (o,h,l) in sorted(legend_items, cmp = lambda x,y: cmp(x[0], y[0]))]
    if draw_legend:
        handles, labels = ax.get_legend_handles_labels()
        handles = handles[:]
        labels = labels[:]
        for h, l in  legend_items + extra_legend_items:
            labels.append(l)
            if type(h)==str:
                h = matplotlib.patches.Rectangle((0, 0), 1, 1, fc=h)
            handles.append(h)
        ax.legend(handles, labels, prop = fp, numpoints = 1, **legend_args)
    #if ax.get_legend() is not None:
    #    map(lambda line: line.set_lw(1.5), ax.get_legend().get_lines())

    if ymin!=None:
        ax.set_ylim(ymin=ymin)
    if ymax!=None:
        ax.set_ylim(ymax=ymax)
    if xmin!=None:
        ax.set_xlim(xmin=xmin)
    if xmax is None:
        if len(histos) > 0:
            xmax = max([max(pd.x) for pd in histos])
            xmax = max([xmax] + [pd.xmax for pd in histos if pd.xmax is not None])
    if xmax is not None:
        ax.set_xlim(xmax=xmax)
    
    if ax_modifier!=None: ax_modifier(ax)
    if outname is not None: fig.savefig(outname)
    del fig
    
def make_stack(pdatas):
    for i in range(len(pdatas)):
        for j in range(i+1, len(pdatas)):
            pdatas[i].y = map(lambda x: x[0] + x[1], zip(pdatas[i].y, pdatas[j].y))
            if pdatas[i].yerrors is not None and pdatas[j].yerrors is not None:
                pdatas[i].yerrors = map(lambda x: math.sqrt(x[0]**2 + x[1]**2), zip(pdatas[i].yerrors, pdatas[j].yerrors))


def scatter_ax_m(x, y, xycol = None, s = 8):
   if xycol is not None:
       if type(xycol)==str:
           color = xycol
       else:
           color = [xycol(xv, yv) for xv, yv in zip(x,y)]
   else:
       color = '#000000'
   return lambda ax: (ax.scatter(x,y, color = color, s = s))


