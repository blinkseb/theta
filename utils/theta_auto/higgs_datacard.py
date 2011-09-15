# make a Model form a datacard as described in 
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideHiggsAnalysisCombinedLimit

import math

from Model import *

# line is either a string or a tuple (string, int)
def get_cmds(line):
    if type(line) == tuple:
        line = line[0]
    cmds = [c.strip() for c in line.split()]
    return [c for c in cmds if c!='']


# in theta, names must begin with a letter and consist only of A-Za-z0-9_-
def transform_name_to_theta(name):
    result = ''
    for c in name:
        if c >= 'a' and c <= 'z' or c >= 'A' and c <= 'Z' or c >='0' and c <='9' or c in ('-', '_'): result += c
        else: result += '_'
    result = result.replace('__', '_')
    if result[0] >= '0' and result[0] <= '9' or result[0]=='-': result = '_' + result
    return result
    
# filter_channel is a function which, for each channel name (as given in the model configuration in fname), returns
# True if this channel should be kept and False otherwise. The default is to keep all channels.
def build_model(fname, filter_channel = lambda chan: True):
    model = Model()
    lines = [l.strip() for l in file(fname)]
    lines = [(lines[i], i+1) for i in range(len(lines)) if not lines[i].startswith('#') and lines[i]!='' and not lines[i].startswith('--')]
    
    cmds = get_cmds(lines[0])
    while cmds[0] != 'imax':
        print 'WARNING: ignoring line %d ("%s") at beginning of file as first token is "%s", not "imax", although not marked as comment' % (lines[0][1], lines[0][0], cmds[0])
        lines = lines[1:]
        cmds = get_cmds(lines[0])
    assert cmds[0]=='imax', "Line %d: Expected imax statement as first statement in the file" % lines[0][1]
    imax = cmds[1]
    if imax !='*': imax = int(imax)
    lines = lines[1:]

    cmds = get_cmds(lines[0])
    assert cmds[0]=='jmax', "Line %d: Expected 'jmax' statement directly after 'imax' statement" % lines[0][1]
    #jmax = int(cmds[1])
    lines = lines[1:]

    cmds = get_cmds(lines[0])
    assert cmds[0]=='kmax', "Line %d: Expected 'kmax' statement directly after 'jmax' statement" % lines[0][1]
    if cmds[1] == '*': kmax = -1
    else: kmax = int(cmds[1])
    lines = lines[1:]

    cmds = get_cmds(lines[0])
    assert cmds[0].lower() in ('bin', 'observation'), "Line %d: Expected 'bin' or 'observation' statement" % lines[0][1]
    if cmds[0].lower() == 'bin':
        # prepend a 'c' so we can use numbers as channel names:
        channel_labels = [ 'c' + c for c in cmds[1:]]
        if imax=='*': imax = len(channel_labels)
        assert len(channel_labels) == imax, "Line %d: Number of processes from 'imax' and number of labels given in 'bin' line (%s) mismatch" % (lines[0][1], str(channel_labels))
        lines = lines[1:]
        cmds = get_cmds(lines[0])
    else:
        channel_labels = [ 'c%d' % i for i in range(1, imax + 1)]
    observables = set(channel_labels)
    assert cmds[0].lower()=='observation', "Line %d: Expected 'observation' statement directly after fist 'bin' statement" % lines[0][1]
    observed_flt = [float(o) for o in cmds[1:]]
    observed_int = map(lambda f: int(f), observed_flt)
    if observed_flt != observed_int: raise RuntimeError, "Line %d: non-integer events given in 'observed' statement!" % lines[0][1]
    if imax=='*': imax = len(observed_int)
    assert len(observed_int) == imax, "Line %d: Number of processes from 'imax' and number of bins given in 'observed' mismatch: imax=%d, given in observed: %d" % (lines[0][1], imax, len(observed))
    for i in range(len(channel_labels)):
        if not filter_channel(channel_labels[i][1:]): continue
        model.set_data_histogram(channel_labels[i], (0.0, 1.0, [observed_flt[i]]))
    lines = lines[1:]

    cmds = get_cmds(lines[0])
    assert cmds[0] == 'bin', "Line %d: Expected 'bin' statement"% lines[0][1]
    # save the channel 'headers', to be used for parsing the next line:
    channels_for_table = cmds[1:]
    for c in channels_for_table:
        if 'c' + c not in channel_labels: raise RuntimeError, "Line % d: unknown channel '%s'" % (lines[0][1], c)
    lines = lines[1:]
    n_cols = len(channels_for_table)

    cmds = get_cmds(lines[0])
    assert cmds[0]=='process'
    processes_for_table = map(transform_name_to_theta, cmds[1:])
    if len(processes_for_table) != n_cols:
        raise RuntimeError, "Line %d: 'bin' statement and 'process' statement have different number of elements" % lines[0][1]
    lines = lines[1:]

    cmds = get_cmds(lines[0])
    assert cmds[0]=='process', "Line %d: Expected second 'process' line directly after first" % lines[0][1]
    process_ids_for_table = [int(s) for s in cmds[1:]]
    if n_cols != len(process_ids_for_table):
        raise RuntimeError, "Line %d: 'process' statements have different number of elements" % lines[0][1]
    lines = lines[1:]

    # check process label / id consistency:
    p_l2i = {}
    p_i2l = {}
    for i in range(n_cols):
        p_l2i[processes_for_table[i]] = process_ids_for_table[i]
        p_i2l[process_ids_for_table[i]] = processes_for_table[i]
    # go through again to make check, also save signal processes:
    signal_processes = set()
    for i in range(n_cols):
        if p_l2i[processes_for_table[i]] != process_ids_for_table[i] or p_i2l[process_ids_for_table[i]] != processes_for_table[i]:
            raise RuntimeError, "Line %d: mapping process id <-> process label (defined via the two 'process' lines) is not one-to-one as expected!" % lines[0][1]
        if p_l2i[processes_for_table[i]] <= 0:
            signal_processes.add(processes_for_table[i])

    cmds = get_cmds(lines[0])
    assert cmds[0]=='rate', "Line %d: Expected 'rate' statement after the two 'process' statements" % lines[0][1]
    if n_cols != len(cmds)-1:
        raise RuntimeError, "Line %d: 'rate' statement does specify the wrong number of elements" % lines[0][1]
    for i in range(n_cols):
        if not filter_channel(channels_for_table[i]): continue
        o, p = 'c%s' % channels_for_table[i], processes_for_table[i]
        n_exp = float(cmds[i+1])
        #print o,p,n_exp
        hf = HistogramFunction()
        hf.set_nominal_histo((0.0, 1.0, [n_exp]))
        model.set_histogram_function(o, p, hf)
    lines = lines[1:]
    
    kmax  = len(lines)

    if kmax != len(lines):
        raise RuntimeError, "Line %d--end: wrong number of lines for systematics (expected kmax=%d, got %d)" % (lines[0][1], kmax, len(lines))
    
    for i in range(kmax):
        cmds = get_cmds(lines[i])
        assert len(cmds) >= len(processes_for_table) + 2, "Line %d: wrong number of entries for uncertainty '%s'" % (lines[0][1], cmds[0])
        uncertainty = transform_name_to_theta('u' + cmds[0])
        if cmds[1] == 'gmN':
            values = cmds[3:]
            n_affected = 0
            k = float(cmds[2])
            for icol in range(n_cols):
                if values[icol]=='-': continue
                val = float(values[icol])
                if val==0.0: continue
                if not filter_channel(channels_for_table[icol]): continue
                obsname ='c%s' % channels_for_table[icol]
                procname = processes_for_table[icol]
                #print cmds[0], k, obsname, procname, values[icol], val
                # add the same parameter (+the factor in the table) as coefficient:
                model.get_coeff(obsname, procname).add_factor('id', parameter = 'delta_%s' % uncertainty)
                n_affected += 1
                n_exp = model.get_histogram_function(obsname, procname).get_nominal_histo()[2][0]
                #print n_exp
                if abs(n_exp - val*k)/max(n_exp, val*k) > 0.03:
                    raise RuntimeError, "gmN uncertainty %s for process %s is inconsistent: the rate expectation should match k*theta but N_exp=%f, k*theta=%f!" % (cmds[0], procname, n_exp, val*k)
            if n_affected > 0:
                k = float(cmds[2])
                hf = HistogramFunction()
                hf.set_nominal_histo((0.0, 1.0, [k]))
                obs_sb = '%s_sideband' % uncertainty
                model.set_histogram_function(obs_sb, 'proc_sb', hf)
                model.set_data_histogram(obs_sb, (0.0, 1.0, [k]))
                model.get_coeff(obs_sb, 'proc_sb').add_factor('id', parameter = 'delta_%s' % uncertainty)
                # the maximum likelihood estimate for the delta parameter is 1.0
                model.distribution.set_distribution('delta_%s' % uncertainty, 'gauss', mean = 1.0, width = float("inf"), range = (0.0, float("inf")))
        elif cmds[1] == 'lnN':
            n_affected = 0
            values = cmds[2:]
            for icol in range(n_cols):
                if values[icol]=='-': continue
                if not filter_channel(channels_for_table[icol]): continue
                if '/' in values[icol]:
                    p = values[icol].find('/')
                    lambda_minus = -math.log(float(values[icol][0:p]))
                    lambda_plus = math.log(float(values[icol][p+1:]))
                else:
                    lambda_minus = math.log(float(values[icol]))
                    lambda_plus = lambda_minus
                obsname ='c%s' % channels_for_table[icol]
                procname = processes_for_table[icol]
                n_affected += 1
                #print cmds[0], obsname, procname, lambda_minus, lambda_plus
                model.get_coeff(obsname, procname).add_factor('exp', parameter = 'delta_%s' % uncertainty, lambda_minus = lambda_minus, lambda_plus = lambda_plus)
            if n_affected > 0:
                model.distribution.set_distribution('delta_%s' % uncertainty, 'gauss', mean = 0.0, width = 1.0, range = (float("-inf"), float("inf")))
        elif cmds[1] == 'gmM':
            values = cmds[2:]
            values_f = set([float(s) for s in values if float(s)!=0.0])
            if len(values_f)>1: raise RunetimeError, "gmM does not support different uncertainties"
            if len(values_f)==0: continue
            n_affected = 0
            for icol in range(n_cols):
                if not filter_channel(channels_for_table[icol]): continue
                obsname ='c%s' % channels_for_table[icol]
                procname = processes_for_table[icol]
                model.get_coeff(obsname, procname).add_factor('id', parameter = 'delta_%s' % uncertainty)
                n_affected += 1
            if n_affected > 0:
                model.distribution.set_distribution('delta_%s' % uncertainty, 'gamma', mean = 1.0, width = float(values[icol]), range = (0.0, float("inf")))
        else: raise RuntimeError, "Line %d: unknown uncertainty type %s" % (lines[0][1], cmds[1])
    #print 'signal_processes', signal_processes
    model.set_signal_processes(list(signal_processes))
    return model

