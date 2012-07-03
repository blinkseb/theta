# -*- coding: utf-8 -*-
import config, sqlite3, os.path, glob, re, array, math, copy
import scipy.stats
import Model
from plotutil import *
import Report

inf = float("inf")
cl_1sigma = 0.68268949213708585
cl_2sigma = 0.95449973610364158

# returns a dictionary float x --> signal process group id
#
# As default, the x value is the first number encountered in the signal process group id string. 
#
# options can contain a dictionary 'spid_to_xvalue' in which case that number is used.
def get_x_to_sp(spgids, **options):
    x_to_sp = {}
    next_x = 0
    for sp in spgids:
        if 'spid_to_xvalue' in options and sp in options['spid_to_xvalue']: x = options['spid_to_xvalue'][sp]
        else: x = extract_number(sp)
        if x is None:
            warning("cannot find x-value for signal process id '%s', using %d (try passing the option spid_to_xvalue = {'%s': XXX})" % (sp, next_x, sp))
            x = next_x
            next_x += 1
        x_to_sp[x] = sp
    assert len(spgids) == len(x_to_sp)
    return x_to_sp


def make_tuple(*args): return tuple(args)


class Copyable:
    def copy(self):
        return copy.deepcopy(self)


# count how often predicate evaluates to true on iterable
def count(pred, iterable):
    result = 0
    for i in iterable:
        if pred(i): result += 1
    return result
    

def reldiff(d1, d2):
    return abs(d1 - d2) / max(abs(d1), abs(d2))

# in theta, names must begin with a letter and consist only of A-Za-z0-9_- and not contain '__'
def transform_name_to_theta(name):
    result = ''
    for c in name:
        if c >= 'a' and c <= 'z' or c >= 'A' and c <= 'Z' or c >='0' and c <='9' or c in ('-', '_'): result += c
        else: result += '_'
    result = result.replace('__', '_')
    if result[0] >= '0' and result[0] <= '9' or result[0]=='-': result = 'tn_' + result
    return result

# get mean and rms estimate from list l:
def get_mean_width(l):
   n = len(l) * 1.0
   assert n > 0
   mean = sum(l) / n
   if n == 1: width = float('inf')
   else: width = math.sqrt(sum([(x - mean)**2 for x in l]) / (n-1))
   return mean, width

# get truncated mean / width similar to MarkovChainMC method from combine:
def get_trunc_mean_width(l):
   n = len(l)
   l_sorted = sorted(l)
   assert n > 0
   median = l_sorted[n/2]
   width = l_sorted[3*n/4] - l_sorted[n/4]
   l2 = [x for x in l if x >= median-width and x <= median+width]
   return get_mean_width(l2)

   
def p_to_Z(p_value):
   return -scipy.stats.norm.ppf(p_value)
   
def Z_to_p(z_value):
    return scipy.stats.norm.sf(z_value)
    
def get_p(n, n0):
    p = n*1.0 / n0
    p_error = max(math.sqrt(p*(1-p) / n0), 1.0 / n0)
    return p, p_error
    
def get_Z(n, n0):
    p, p_error = get_p(n, n0)
    Z0 = p_to_Z(p)
    Zplus = p_to_Z(p + p_error)
    if p_error >= p: return Z0, inf
    else: Zminus = p_to_Z(p - p_error)
    return Z0, max(abs(Z0 - Zplus), abs(Z0 - Zminus))
   
  
def extract_number(s):
    r = re.compile('(\d+)')
    m = r.search(s)
    if m is None: return None
    return float(m.group(1))

## \brief Populate an instance of plotutil.plotdata with results from bayesian_quantiles or cls_limits.
#
# The parameter 'quantiles' is the return value from \ref bayesian_quantiles or cls_limits
#
# include_band is a boolean indicating whether to also include the +-1sigma and +-2sigma bands, or only the median line.
#
# options:
# - spid_to_xvalue: a dictionary mapping signal process name to values to be used as x axis for the band plot. As default, the first integer
#     in the signal process name is used for the x axis value.
#
# returns one plotutil.plotdata instance containing the 'observed' (or median expected) limit and 'expected' bands.
def limit_band_plot(quantiles, include_band, quantile = 0.95, **options):
    #expected results maps (process name) -> (median, band1, band2)
    # where band1 and band2 are tuples for the central 68 and 95%, resp.
    results = {}
    signal_processes = set()
    for sp in quantiles:
        # ignore uncertainties or other special entries:
        if '__' in sp: continue
        signal_processes.add(sp)
    # map process names and x values:
    x_to_sp = get_x_to_sp(signal_processes, **options)
    pd = plotdata()
    pd.color = '#000000'
    if include_band: pd.color = '#aaaaaa'
    pd.as_function = True
    pd.x = sorted(list(x_to_sp.keys()))
    pd.y = []
    if include_band:
        pd.bands = [([], [], '#00ff00'), ([], [], '#00aa00')]
    else: pd.yerrors = []
    for x in pd.x:
        sp = x_to_sp[x]
        data = sorted(quantiles[sp][quantile])
        n = len(data)
        median, band1, band2 = (data[n / 2], (data[int(0.16 * n)], data[int(0.84 * n)]), (data[int(0.025 * n)], data[int(0.975 * n)]))
        if not include_band:
            mean, error = get_trunc_mean_width(data)
            pd.y.append(mean)
            pd.yerrors.append(error / math.sqrt(len(data)))
            continue
        pd.y.append(median)
        pd.bands[1][0].append(band1[0])
        pd.bands[1][1].append(band1[1])
        pd.bands[0][0].append(band2[0])
        pd.bands[0][1].append(band2[1])
    return pd
    

## \brief Make expected / observed limits plots and write them to the report
#
# name and shortname are used to distinguish different calls to this function. name is used in the report
# and should be the 'human-firendly' version, while shortname is used for filenames and should be the 'computer-friendly'
# version.
def report_limit_band_plot(expected_limits, observed_limits, name, shortname, write_table = True):
    plotsdir = os.path.join(config.workdir, 'plots')
    plots = []
    extra_legend_items = []
    if expected_limits is not None:
        expected_limits.legend = 'median expected limit'
        extra_legend_items.append((expected_limits.bands[0][2], '$\\pm 2\\sigma$ expected limit'))
        extra_legend_items.append((expected_limits.bands[1][2], '$\\pm 1\\sigma$ expected limit'))
        plots.append(expected_limits)
    if observed_limits is not None:
        observed_limits.legend = 'observed limit'
        plots.append(observed_limits)
    if len(plots) == 0: return
    config.report.new_section('Limits %s' % name)
    if write_table:
        result_table = Report.table()
        result_table.add_column('process', 'signal process')
        if expected_limits is not None:
            result_table.add_column('exp', 'expected limit')
            result_table.add_column('exp1', 'expected limit (central 1sigma)')
            result_table.add_column('exp2', 'expected limit (central 2sigma)')
        if observed_limits is not None:
            result_table.add_column('obs', 'observed limit')
        x_values = []
        if expected_limits is not None: x_values = expected_limits.x
        else: x_values = observed_limits.x
        for i in range(len(x_values)):
            result_table.set_column('process', '%g' % x_values[i])
            if expected_limits is not None:
                result_table.set_column('exp', '%.3g' % expected_limits.y[i])
                result_table.set_column('exp1', '%.3g--%.3g' % (expected_limits.bands[1][0][i], expected_limits.bands[1][1][i]))
                result_table.set_column('exp2', '%.3g--%.3g' % (expected_limits.bands[0][0][i], expected_limits.bands[0][1][i]))
            if observed_limits is not None:
                result_table.set_column('obs', '%.3g +- %.3g' % (observed_limits.y[i], observed_limits.yerrors[i]))
            result_table.add_row()
        config.report.add_html(result_table.html())
    plot(plots, 'signal process', 'upper limit', os.path.join(plotsdir, 'limit_band_plot-%s.png' % shortname), extra_legend_items=extra_legend_items)
    plot(plots, 'signal process', 'upper limit', os.path.join(plotsdir, 'limit_band_plot-log-%s.png' % shortname), logy = True, extra_legend_items=extra_legend_items)
    config.report.add_html('<p><img src="plots/limit_band_plot-%s.png" /></p>' % shortname)
    config.report.add_html('<p><img src="plots/limit_band_plot-log-%s.png" /></p>' % shortname)
    return plots

def info(s):
    if not config.suppress_info:
        print "[INFO] ", s

def warning(s):
    print "[WARN] ", s

# return a list of result rows for the given query on the .db filename.
def sql_singlefile(filename, query, return_colnames = False):
    if not os.path.exists(filename): raise RuntimeError, "sql: the file %s does not exist!" % filename
    conn = sqlite3.connect(filename)
    c = conn.cursor()
    try:  c.execute(query)
    except Exception, ex:
        print "exception executing %s on file %s: %s" % (query, filename, str(ex))
        raise ex
    result = c.fetchall()
    if return_colnames:
        desc = c.description
    c.close()
    conn.close()
    if return_colnames: return result, desc
    return result

def sql(filename_pattern, query):
    result = []
    for f in glob.glob(filename_pattern):
        result.extend(sql_singlefile(f, query))
    return result
