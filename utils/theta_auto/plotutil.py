# coding=utf8

import matplotlib
# note: png with cairo is broken on ubuntu 10.10: segfaults sometimes in figure.save ...
# but not using it leads to VERY slow plotting, especially over ssh!
#matplotlib.use('Cairo')

import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import matplotlib.text
import matplotlib.lines
import matplotlib.patches

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
    def __init__(self):
        self.x = []
        self.y = []
        self.legend = None
        self.yerrors = None
        self.xerrors = None
        self.fill_color = None
        self.color = None
        self.marker = 'None'
        self.lw = 2
        self.fmt = '-'
        # an array of bands; a band is a three-tuple (y1, y2, color). y1 and y2 are 
        # arrays of y values.
        # bands to draw first should come first 
        self.bands = None
        self.band_lw = 0
        self.bands_fill = True
        self.as_function = False
        self.draw_line = True
        
    # make a histogram of the given values
    def histogram(self, values, xmin, xmax, nbins):
        self.x = [xmin + (xmax - xmin) / nbins * i for i in range(nbins)]
        self.y = [0.0] * nbins
        for v in values:
            ibin = int((v - xmin) / (xmax - xmin) * nbins)
            if ibin < 0 or ibin >= nbins: continue
            self.y[ibin] += 1
        
    # ofile is a string (filename) or a handle to a open file
    def write_txt(self, ofile):
        if type(ofile)==str: ofile = open(ofile, 'w')
        ofile.write('# x; y')
        if self.bands is not None:
            for k in range(len(self.bands)):
                ofile.write('; band %d low; band %d high' % (k, k))
        ofile.write("\n")
        for i in range(len(self.x)):
            ofile.write("%10.5g %10.5g " % (self.x[i], self.y[i]))
            if self.bands is not None:
                for k in range(len(self.bands)):
                    ofile.write("%10.5g %10.5g" % (self.bands[k][0][i], self.bands[k][1][i]))
            ofile.write("\n")

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
def plot(histos, xlabel, ylabel, outname, logy = False, logx = False, ax_modifier=None, title_ul=None, title_ur = None,
extra_legend_items = [], xmin = None, xmax=None, ymin=None, ymax=None):
    cm = 1.0/2.54
    fsize = 15*cm, 12*cm
    fp = fm.FontProperties(size=10)
    fig = plt.figure(figsize = fsize)
    rect = 0.15, 0.15, 0.8, 0.75
    ax = fig.add_axes(rect)
    if logy:  ax.set_yscale('log')
    if logx: ax.set_xscale('log')
    add_xlabel(ax, xlabel, fontproperties=fp)
    add_ylabel(ax, ylabel, fontproperties=fp)
    if title_ul is not None: ax.text(0.0, 1.02, title, transform = ax.transAxes, ha='left', va='bottom')
    if title_ur is not None: ax.text(1.0, 1.02, title_ur, transform = ax.transAxes, ha='right', va='bottom')
    draw_legend = False
    for histo in histos:
        assert len(histo.x)==len(histo.y), "number of x,y coordinates not the same for '%s'" % histo.legend
        if histo.legend: draw_legend = True
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
            # histo.x is assumed to contain the lower bin edges in this case ...
            if len(histo.x) >= 2:  x_binwidth = histo.x[1] - histo.x[0]
            else: x_binwidth = 1.0
            # if histo.yerrors is set, draw with errorbars, shifted by 1/2 binwidth ...
            if histo.yerrors is not None:
               new_x = [x + 0.5 * x_binwidth for x in histo.x]
               ax.errorbar(new_x, histo.y, histo.yerrors, label=histo.legend, ecolor = histo.color, marker='o', ms = 2, capsize = 0, fmt = None)
            else:
               new_x = [histo.x[0]]
               for x in histo.x[1:]: new_x += [x]*2
               new_x += [histo.x[-1] + x_binwidth]
               new_y = []
               for y in histo.y: new_y += [y]*2
               if histo.fill_color is not None:
                   ax.fill_between(new_x, new_y, [0] * len(new_y), lw=histo.lw, label=histo.legend, color=histo.color, facecolor = histo.fill_color)
               else:
                   ax.plot(new_x, new_y, histo.fmt, lw=histo.lw, label=histo.legend, color=histo.color)
               #ax.bar(map(lambda x: x - 0.5  * x_binwidth, histo.x), histo.y, width=x_binwidth, color=histo.color, label=histo.legend, ec=histo.color)
        else:
            if histo.yerrors is not None:
                lw = histo.lw
                if histo.draw_line is False: lw = 0
                ax.errorbar(histo.x, histo.y, histo.yerrors, elinewidth = histo.lw, lw=lw, label=histo.legend, color=histo.color, marker=histo.marker)
            else:
                ax.plot(histo.x, histo.y, histo.fmt, lw=histo.lw, label=histo.legend, color=histo.color, marker=histo.marker)

    if draw_legend:
        handles, labels = ax.get_legend_handles_labels()
        handles = handles[:]
        labels = labels[:]
        for h, l in extra_legend_items:
            labels.append(l)
            if type(h)==str:
                h = matplotlib.patches.Rectangle((0, 0), 1, 1, fc=h)
            handles.append(h)
        ax.legend(handles, labels, prop = fp)
    if ax.get_legend() is not None:
        map(lambda line: line.set_lw(1.5), ax.get_legend().get_lines())

    if ymin!=None:
        ax.set_ylim(ymin=ymin)
    if ymax!=None:
        ax.set_ylim(ymax=ymax)
    if xmin!=None:
        ax.set_xlim(xmin=xmin)
    if xmax!=None:
        ax.set_xlim(xmax=xmax)
    
    if ax_modifier!=None: ax_modifier(ax)
    fig.savefig(outname)
    del fig
    
def make_stack(pdatas):
    for i in range(len(pdatas)):
        for j in range(i+1, len(pdatas)):
            pdatas[i].y = map(lambda x: x[0] + x[1], zip(pdatas[i].y, pdatas[j].y))

