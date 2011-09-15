import config, sqlite3, os.path, glob, re
import numpy.linalg, array, math

import theta_interface, plotutil

import Model

# returns a Distribution object, given the model, signal process and nuisance_prior specification ('shape:X;rate:Y'...)
def nuisance_prior_distribution(model, spec):
    if spec.__class__ == Model.Distribution: result = Model.Distribution.merge(spec, spec)  # note: merging copies ...
    else:
        result = Model.Distribution()
        min_p0, max_p1 = len(spec), 0
        shape_spec, rate_spec = None, None
        if 'shape:' in spec:
            p0 = spec.find('shape:')
            p1 = spec.find(';', p0)
            if p1 == -1: p1 = len(spec)
            min_p0 = min(min_p0, p0)
            max_p1 = max(max_p1, p1)
            shape_spec = spec[p0 + len('shape:'):p1]
        if 'rate:' in spec:
            p0 = spec.find('rate:')
            p1 = spec.find(';', p0)
            if p1 == -1: p1 = len(spec)
            min_p0 = min(min_p0, p0)
            max_p1 = max(max_p1, p1)
            rate_spec = spec[p0 + len('rate:'):p1]
        if (min_p0, max_p1) != (0, len(spec)): raise RuntimeError, 'nuisance_prior_distribution: could not parse specification %s (%d, %d)' % (spec, min_p0, max_p1)
        rc_pars, sc_pars = model.get_rate_shape_parameters()
        if shape_spec is not None:
            if shape_spec not in ('fix', 'free'): raise RuntimeError, 'unknown shape specification "%s"' % shape_spec
            for par in sc_pars:
                dist_pars = model.distribution.get_distribution(par)
                if shape_spec == 'fix':
                    #print 'fixing %s' % par
                    dist_pars['width'] = 0.0
                    result.set_distribution(par, **dist_pars)
                elif shape_spec == 'free':
                    #print 'freeing %s' % par
                    dist_pars['width'] = float("inf")
                    result.set_distribution(par, **dist_pars)
        if rate_spec is not None:
            if rate_spec not in ('fix', 'free'): raise RuntimeError, 'unknown rate specification "%s"' % rate_spec
            for par in rc_pars:
                dist_pars = model.distribution.get_distribution(par)
                if rate_spec == 'fix':
                    #print 'fixing %s' % par
                    dist_pars['width'] = 0.0
                    result.set_distribution(par, **dist_pars)
                elif rate_spec == 'free':
                    #print 'freeing %s' % par
                    dist_pars['width'] = float("inf")
                    result.set_distribution(par, **dist_pars)
    result = Model.Distribution.merge(model.distribution, result)
    return result

def get_mean_width(l):
   n = len(l)
   assert n > 0
   mean = sum(l) / n
   if n == 1: width = float('inf')
   else: width = math.sqrt(sum([(x - mean)**2 for x in l]) / (n-1))
   return mean, width

# returns the theta config dictionary, given a signal_prior specification ('flat' / 'fix:X')
def signal_prior_dict(spec):
    if type(spec) == str:
        if spec.startswith('flat'):
            if spec.startswith('flat:'):
                res = re.match('flat:\[([^,]+),(.*)\]', spec)
                if res is None: raise RuntimeError, "signal_prior specification '%s' invalid (does not match range specification syntax)" % spec
                xmin, xmax = float(res.group(1)), float(res.group(2))
            else:
                if spec!='flat': raise RuntimeError, "signal_prior specification '%s' invalid" % spec
                xmin, xmax = 1e-12, float("inf")
            value = 0.5 * (xmax - xmin)
            if value==float("inf"): value = 1.0
            signal_prior_dict = {'type': 'flat_distribution', 'beta_signal': {'range': [xmin, xmax], 'fix-sample-value': value}}
        elif spec.startswith('fix:'):
            v = float(spec[4:])
            signal_prior_dict = theta_interface.delta_distribution(beta_signal = v)
        else: raise RuntimeError, "signal_prior specification '%s' unknown" % spec
    else:
        if type(spec) != dict: raise RuntimeError, "signal_prior specification has to be a string ('flat' / 'fix:X') or a dictionary!"
        signal_prior_dict = spec
    return signal_prior_dict

# returns a tuple (data_source dict, beta_signal distribution dict), given the input_spec ('toys:X' / 'toys-scan:N,min,max' / 'data')
# the name is always 'source'.
#
# possible options: 'toydata_seed'
def data_source_dict(model, input_spec, **options):
    if input_spec == 'data':
        source_dict = {'type': 'histo_source', 'name': 'source'}
        for o in model.observables:
            source_dict[o] = theta_interface.get_histo_cfg(model.data_histos[o])
        return (source_dict, theta_interface.delta_distribution(beta_signal = 0.0))
    elif input_spec.startswith('toys:') or input_spec.startswith('toys-asimov:'):
        beta_signal_value = float(input_spec[input_spec.find(':') + 1:])
        seed = 1
        if 'toydata_seed' in options: seed = int(options['toydata_seed'])
        result = {'type': 'model_source', 'name': 'source', 'model': '@model', 'rnd_gen': {'seed': seed}}
        if input_spec.startswith('toys-asimov:'): result['dice_poisson'] = False;
        return (result, theta_interface.delta_distribution(beta_signal = beta_signal_value))
    else: raise RuntimeError, 'data_source dict: not implemented for input = "%s"' % input_spec

# given (x,y) points of a function find the x value where the maximum is taken, suitable
# especially for binned data.
# xs must be monotonically increasing.
# Uses quadratic interpolation with the three points around the maximum to make the result more accurate.
# If the interpolated maximum is found outside the x range, the returned x value will be cut off at "half a binwidth", i.e. at
#         xs[0] - 0.5 * (xs[1] - xs[0])       and        xs[-1] + 0.5 * (xs[-1] - xs[-2])
# Only works well for non-noisy data.
def argmax_xy(xs, ys):
    assert len(xs) == len(ys)
    assert len(xs) >= 3
    maxy = -float("inf")
    maxi = 0
    for i in range(len(xs)):
        if ys[i] > maxy:
            maxy = ys[i]
            maxi = i
    if maxi==0: imin, imax = 0, 3
    elif maxi==len(xs)-1: imin, imax = len(xs)-3, len(xs)
    else: imin, imax = maxi-1, maxi+2
    xleft, x0, xright = xs[imin:imax]
    yleft, y0, yright = ys[imin:imax]
    a,b,c = numpy.linalg.solve([[xleft**2, xleft, 1], [x0**2, x0, 1], [xright**2, xright, 1]], [yleft, y0, yright])
    if a>=0: return xs[maxi]
    return -b/(2*a)
    
# same as argmax_xy, but use an instance of poltutil.plotdata as argument
def argmax(pd): return argmax_xy(pd.x, pd.y)

def extract_number(s):
    r = re.compile('(\d+)')
    m = r.search(s)
    if m is None: return None
    return float(m.group(1))


def info(s):
    if not config.suppress_info:
        print "[INFO] ", s

# return a list of result rows for the given query on the .db filename.
def sql_singlefile(filename, query):
    if not os.path.exists(filename): raise RuntimeError, "sql: the file %s does not exist!" % filename
    conn = sqlite3.connect(filename)
    c = conn.cursor()
    try:  c.execute(query)
    except Exception, ex:
        print "exception executing %s on file %s: %s" % (query, filename, str(ex))
        raise ex
    result = c.fetchall()
    c.close()
    conn.close()
    return result

def sql(filename_pattern, query):
    result = []
    for f in glob.glob(filename_pattern):
        result.extend(sql_singlefile(f, query))
    return result

def plotdata_from_histoColumn(h):
    a = array.array('d')
    a.fromstring(h)
    pdata = plotutil.plotdata()
    xmin = a[0]
    xmax = a[1]
    nbins = len(a) - 4
    binwidth = (xmax - xmin) / nbins
    pdata.x = [xmin + i*binwidth for i in range(nbins)]
    pdata.y = a[3:-1]
    pdata.color = '#000000'
    return pdata

