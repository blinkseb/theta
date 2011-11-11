# all classes here are immutable: once constructed, the contents cannot be changed.
# (instead, the methods return new instances).

class FunctionBase:
    def get_cfg(self): pass
    def get_parameters(self): pass

    def __eq__(self, other): return type(self) is type(other) and self.__dict__ == other.__dict__
    def __neq__(self, other): return not (self.__eq__(other))

    def __hash__(self): return hash(self.get_cfg())

    def __mul__(self, other):
        if type(other) in (float, int): other = MultiplyFunction([float(other)])
        assert isinstance(other, FunctionBase)
        return MultiplyFunction([self, other])

    def __add__(self, other):
        if type(other) in (float, int): other = MultiplyFunction([float(other)])
        assert isinstance(other, FunctionBase)
        return AddFunction([self, other])


class MultiplyFunction(FunctionBase):
    # factors is a list of FunctionBase instances, strings (parameter names), or convertable to float.
    def __init__(self, factors):
        parameters = set()
        ffactors = []
        pfactors = []
        constant = 1.0
        for f in factors:
            if type(f) == str:
                pfactors.append(f)
                parameters.add(f)
            elif isinstance(f, MultiplyFunction):
                ffactors.extend(f.function_factors)
                pfactors.extend(f.parameter_factors)
                parameters.update(f.get_parameters())
                constant *= f.constant
            elif isinstance(f, FunctionBase):
                ffactors.append(f)
                parameters.update(f.get_parameters())
            else:
                constant *= float(f)
        self.parameters = frozenset(parameters)
        self.function_factors = tuple(ffactors)
        self.parameter_factors = tuple(pfactors)
        self.constant = constant

    def get_cfg(self):
        result = {'type': 'multiply', 'factors': [f.get_cfg() for f in self.function_factors] + [p for p in self.parameter_factors]}
        if self.constant != 1.0: result['factors'].append(self.constant)
        return result

    def __str__(self):
        terms = []
        terms.extend(self.parameter_factors)
        terms.extend(map(str, self.function_factors))
        if self.constant != 1.0 or len(terms)==0:
            terms.append('%f' % self.constant)
        return ' * '.join(terms)

    def get_parameters(self):
        return self.parameters
        
    def __mul__(self, other):
        res = FunctionBase.__mul__(self, other)
        return res._get_optimized()

    # optimize the exp_functions (if more than one) by using only one ExpFunction which calculates the same result.
    # returns a new MultiplyFunction instance
    def _get_optimized(self):
        factors = list(self.parameter_factors)
        factors.append(self.constant)
        exp = None
        for f in self.function_factors:
            if isinstance(f, ExpFunction):
                if exp is None: exp = f
                else: exp = exp * f
            else: factors.append(f)
        if exp is not None: factors.append(exp)
        return MultiplyFunction(factors)


# add some other functions
class AddFunction(FunctionBase):
    def __init__(self, funcs):
        parameters = set()
        functions = []
        for f in funcs:
            parameters.update(f.get_parameters())
            if isinstance(f, AddFunction): functions.extend(f.functions)
            else: functions.append(f)
        self.functions = tuple(functions)
        self.parameters = frozenset(parameters)

    def get_cfg(self):
        result = {'type': 'add', 'addends': [f.get_cfg() for f in self.functions]}
        return result

    def get_parameters(self):
        return self.parameters

    def __str__(self):
        return '(' + ' + '.join(map(str, self.functions)) + ')'


# sqrt(w1*p1**2 + w2*p2**2 + ...)
class AddSquared(FunctionBase):
    def __init__(self, pars, weights):
        assert len(pars) == len(weights)
        self.parameters = tuple([str(p) for p in pars])
        self.weights = tuple([float(w) for w in weights])

    def get_cfg(self):
        result = {'type': 'add_squared', 'parameters': self.parameters, 'weights': self.weights}
        return result

    def get_parameters(self):
        return self.parameters

# some methods working with FunctionBase objects and strings:
def get_cfg(f):
    if type(f) == str: return f
    else: return f.get_cfg()
    
def get_parameters(f):
    if type(f) == str: return set([f])
    else: return f.get_parameters()

def close(a, b):
    if a*b==0: return a==0 and b==0
    return abs(a-b) / max(abs(a), abs(b)) < 1e-10

# the negative logarithm of a multivariate Gauss
class NLGauss(FunctionBase):
    # rows is a list / tuple, etc. of functions (or parameters). mu is a list of floats, covariance is a list of list of floats.
    def __init__(self, rows, mu, covariance):
        assert len(rows) == len(mu) and len(rows) == len(covariance)
        n = len(rows)
        self.rows = tuple(rows)
        self.mu = tuple([float(m) for m in mu])
        cov = []
        parameters = set()
        for r in self.rows:
            parameters.update(get_parameters(r))
        self.parameters = frozenset(parameters)
        for irow in range(n):
            row = [float(e) for e in covariance[irow]]
            assert len(row) == n
            cov.append(row)
        #print cov
        for i in range(n):
            for j in range(i, n):
                assert close(cov[i][j], cov[j][i]), "covariance not symmetric (i,j) = (%d, %d); elements: %f, %f" % (i,j,cov[i][j], cov[j][i])
                symm_elem = 0.5 * (cov[i][j] + cov[j][i])
                cov[i][j] = symm_elem
                cov[j][i] = symm_elem
        self.covariance = tuple([tuple(row) for row in cov])

    def get_cfg(self):
        result = {'type': 'nl_gauss', 'rows': [get_cfg(r) for r in self.rows], 'mu': self.mu, 'covariance': self.covariance}
        return result
    """
    def __add__(self, other):
        if isinstance(other, NLGauss):
            overlapping_rows = set(self.rows).intersection(other.rows)
            new_rows = list(overlapping_rows)
            for r in self.rows:
                if r not in overlapping_rows: new_rows.append(r)
            for r in other.rows:
                if r not in overlapping_rows: new_rows.append(r)
            n = len(new_rows)
            for i in range(n):
                r = new_rows[i]
                if r in self.rows and r in other.rows: 
                    
        else: FunctionBase.__add__(self, other)
    """

    def get_parameters(self): return self.parameters

class ExpFunction(FunctionBase):
    def __init__(self, parameters, lambdas_plus, lambdas_minus):
        self.parameters = tuple(parameters)
        self.lambdas_plus = tuple(lambdas_plus)
        self.lambdas_minus = tuple(lambdas_minus)

    def get_lambdas_plus(self): return self.lambdas_plus
    def get_lambdas_minus(self): return self.lambdas_minus
    def get_parameters(self): return self.parameters

    def __mul__(self, other):
        if isinstance(other, ExpFunction):
            factors = []
            # parameter name to list [lambda_minus, lambda_plus]
            parameters_to_lmp = {}
            pars = list(other.parameters) + list(self.parameters)
            lambdas_plus = list(other.lambdas_plus) + list(self.lambdas_plus)
            lambdas_minus = list(other.lambdas_minus) + list(self.lambdas_minus)
            for i in range(len(pars)):
                p = pars[i]
                if p in parameters_to_lmp:
                    parameters_to_lmp[p][0] += lambdas_minus[i]
                    parameters_to_lmp[p][1] += lambdas_plus[i]
                else: parameters_to_lmp[p] = [lambdas_minus[i], lambdas_plus[i]]
            parameters = []
            lambdas_plus = []
            lambdas_minus = []
            for p in parameters_to_lmp:
                parameters.append(p)
                lambdas_minus.append(parameters_to_lmp[p][0])
                lambdas_plus.append(parameters_to_lmp[p][1])
            return ExpFunction(parameters, lambdas_plus, lambdas_minus)
        else: return FunctionBase.__mul__(self, other)

    def __str__(self):
        terms = []
        for i in range(len(self.parameters)):
            if self.lambdas_plus[i] == self.lambdas_minus[i]: terms.append(self.parameters[i] + ' * %f' % self.lambdas_plus[i])
            else:
                terms.append('%s * (%s > 0 ? : %f : %f)' % (self.parameters[i], self.parameters[i], self.lambdas_plus[i], self.lambdas_minus[i]))
        return 'exp(' + ' + '.join(terms) + ')'

    def get_cfg(self):
        return {'type': 'exp_function', 'parameters': self.parameters, 'lambdas_plus': self.lambdas_plus, 'lambdas_minus': self.lambdas_minus}


