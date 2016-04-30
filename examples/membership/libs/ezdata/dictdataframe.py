'''
DictDataFrame, a simplistic column based dataframe

The :class:`DataFrame` container allows easier manipulations of the data but is
basically a wrapper of many existing function around a dictionary object.

.. note::

    * tested with python 2.7, & 3.4
    * tested compatible with pandas (not required)
    * requirements: numpy

:author: Morgan Fouesneau
'''
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np
import sys
import itertools
import operator

PY3 = sys.version_info[0] > 2

if PY3:
    iteritems = operator.methodcaller('items')
    itervalues = operator.methodcaller('values')
    basestring = (str, bytes)
else:
    range = xrange
    from itertools import izip as zip
    iteritems = operator.methodcaller('iteritems')
    itervalues = operator.methodcaller('itervalues')
    basestring = (str, unicode)

try:
    from .plotter import Plotter
except ImportError:
    Plotter = None


__all__ = ['DictDataFrame']


def pretty_size_print(num_bytes):
    """
    Output number of bytes in a human readable format

    Parameters
    ----------
    num_bytes: int
        number of bytes to convert

    returns
    -------
    output: str
        string representation of the size with appropriate unit scale
    """
    if num_bytes is None:
        return

    KiB = 1024
    MiB = KiB * KiB
    GiB = KiB * MiB
    TiB = KiB * GiB
    PiB = KiB * TiB
    EiB = KiB * PiB
    ZiB = KiB * EiB
    YiB = KiB * ZiB

    if num_bytes > YiB:
        output = '%.3g YB' % (num_bytes / YiB)
    elif num_bytes > ZiB:
        output = '%.3g ZB' % (num_bytes / ZiB)
    elif num_bytes > EiB:
        output = '%.3g EB' % (num_bytes / EiB)
    elif num_bytes > PiB:
        output = '%.3g PB' % (num_bytes / PiB)
    elif num_bytes > TiB:
        output = '%.3g TB' % (num_bytes / TiB)
    elif num_bytes > GiB:
        output = '%.3g GB' % (num_bytes / GiB)
    elif num_bytes > MiB:
        output = '%.3g MB' % (num_bytes / MiB)
    elif num_bytes > KiB:
        output = '%.3g KB' % (num_bytes / KiB)
    else:
        output = '%.3g Bytes' % (num_bytes)

    return output


class DictDataFrame(dict):
    """
    A simple-ish dictionary like structure allowing usage as array on non
    constant multi-dimensional column data.

    It initializes like a normal dictionary and can be used as such.
    A few divergence points though: some default methods such as :func:`len`
    may refer to the lines and not the columns as a normal dictionary would.

    This data object implements also array slicing, shape, dtypes and some data
    functions (sortby, groupby, where, select, etc)
    """
    def __init__(self, *args, **kwargs):
        """ A dictionary constructor and attributes declaration """
        dict.__init__(self, *args, **kwargs)
        self.__dict__ = self   # give access to everything directly

    def __len__(self):
        """ Returns the number of rows """
        return len(self[list(self.keys())[0]])

    @property
    def nrows(self):
        """ Number of rows in the dataset """
        return len(self)

    @property
    def ncols(self):
        """ Number of columns in the dataset """
        return dict.__len__(self)

    @classmethod
    def from_lines(cls, it):
        """ Create a DataFrame object from its lines instead of columns

        Parameters
        ----------
        it: iterable
            sequence of lines with the same keys (expecting dict like structure)

        Returns
        -------
        df: DataFrame
            a new object
        """
        d = {}
        n = 0
        for line in it:
            for k in line.keys():
                d.setdefault(k, [np.atleast_1d(np.nan)] * n).append(np.atleast_1d(line[k]))
        for k,v in dict.items(d):
            dict.__setitem__(d, k, np.squeeze(np.vstack(v)))

        return cls(d)

    def __repr__(self):
        txt = 'DataFrame ({0:s})\n'.format(pretty_size_print(self.nbytes))
        txt += '\n'.join([str((k, v.dtype, v.shape)) for (k,v) in self.items()])
        return txt

    @property
    def nbytes(self):
        """ number of bytes of the object """
        n = sum(k.nbytes if hasattr(k, 'nbytes') else sys.getsizeof(k)
                for k in self.__dict__.values())
        return n

    def __getitem__(self, k):
        try:
            return dict.__getitem__(self, k)
        except Exception as e:
            print(e)
            return self.__class__({a:v[k] for a,v in self.items()})

    @property
    def dtype(self):
        """ the dtypes of each column of the dataset """
        return dict((k, v.dtype) for (k,v) in self.items())

    @property
    def shape(self):
        """ dict of shapes """
        return dict((k, v.shape) for (k,v) in self.items())

    def groupby(self, key):
        """ create an iterator which returns (key, DataFrame) grouped by each
        value of key(value) """
        for k, index in self.arg_groupby(key):
            d = {a: b[index] for a,b in self.items()}
            yield k, self.__class__(d)

    def arg_groupby(self, key):
        """ create an iterator which returns (key, index) grouped by each
        value of key(value) """
        val = self.evalexpr(key)
        ind = sorted(zip(val, range(len(val))), key=lambda x:x[0])

        for k, grp in itertools.groupby(ind, lambda x: x[0]):
            index = [k[1] for k in grp]
            yield k, index

    def __iter__(self):
        """ Iterator on the lines of the dataframe """
        return self.iterlines()

    def iterlines(self):
        """ Iterator on the lines of the dataframe """
        return self.lines()

    def lines(self):
        """ Iterator on the lines of the dataframe """
        for k in range(self.nrows):
            yield self[k]

    def columns(self):
        """ Iterator on the columns
        refers to :func:`dict.items`
        """
        return dict.items(self)

    def where(self, condition, condvars=None, **kwargs):
        """ Read table data fulfilling the given `condition`.
        Only the rows fulfilling the `condition` are included in the result.

        Parameters
        ----------
        query: generator
            generator of records from a query

        condition : str
            The evaluation of this condition should return True or False the
            condition can use field names and variables defined in condvars

        condvars: dict
            dictionary of variables necessary to evaluate the condition.

        Returns
        -------
        it: generator
            iterator on the query content. Each iteration contains one selected
            entry.

        ..note:
            there is no prior check on the variables and names
        """
        for line in self.lines():
            if eval(condition, dict(line), condvars):
                yield line

    def sortby(self, key, reverse=False, copy=False):
        """
        Parameters
        ----------
        key: str
            key to sort on. Must be in the data

        reverse: bool
            if set sort by decreasing order

        copy: bool
            if set returns a new dataframe

        Returns
        -------
        it: DataFrame or None
            new dataframe only if copy is False
        """
        val = self.evalexpr(key)
        ind = np.argsort(val)
        if reverse:
            ind = ind[::-1]
        if not copy:
            for k in self.keys():
                dict.__setitem__(self, k, dict.__getitem__(self, k)[ind])
        else:
            d = {}
            for k in self.keys():
                d[k] = dict.__getitem__(self, k)[ind]
            return self.__class__(d)

    def select(self, keys, caseless=False):
        """ Read table data but returns only selected fields.

        Parameters
        ----------
        keys: str, sequence of str
            field names to select.
            Can be a single field names as follow:
            'RA', or ['RA', 'DEC'], or 'RA,DEC', or 'RA DEC'

        caseless: bool
            if set, do not pay attention to case.

        Returns
        -------
        df: DataFrame
            reduced dataframe

        ..note:
            there is no prior check on the variables and names
        """
        if keys == '*':
            return self

        if caseless:
            _keys = ''.join([k.lower() for k in keys])
            df = self.__class__(dict( (k,v) for k,v in self.items() if (k.lower() in _keys)))
        else:
            df = self.__class__(dict( (k,v) for k,v in self.items() if k in keys))
        return df

    def pickle_dump(self, fname):
        """ create a pickle dump of the dataset """
        import pickle
        with open(fname, 'wb') as f:
            pickle.dump(self, f)

    @classmethod
    def unpickle(cls, fname):
        """ restore a previously pickled object """
        import pickle
        with open(fname, 'rb') as f:
            return pickle.load(f)

    def multigroupby(self, key, *args):
        """
        Returns nested grouped DataFrames given the multiple keys

        Parameters
        ----------
        key1, key2, ...: sequence
            keys over which indexing the data

        Returns
        -------
        it: generator
            (key1, (key2, (... keyn, {})))
        """
        return _df_multigroupby(self, key, *args)

    def aggregate(self, func, keys, args=(), kwargs={}):
        """ Apply func on groups within the data

        Parameters
        ----------
        func: callable
            function to apply
        keys: sequence(str)
            sequence of keys defining the groups
        args: tuple
            optional arguments to func (will be given at the end)
        kwargs: dict
            optional keywords to func

        Returns
        -------
        seq: sequence
            flattened sequence of keys and value
            (key1, key2, ... keyn, {})
        """
        pv = [(k, list(v)) for k, v in self.multigroupby(*keys)]
        return _df_multigroupby_aggregate(pv, func=func)

    @property
    def Plotter(self):
        """ Plotter instance related to this dataset.
        Requires plotter add-on to work """
        if Plotter is None:
            raise AttributeError('the add-on was not found, this property is not available')
        else:
            return Plotter(self)

    def evalexpr(self, expr, exprvars=None, dtype=float):
        """ evaluate expression based on the data and external variables
            all np function can be used (log, exp, pi...)

        Parameters
        ----------
        data: dict or dict-like structure
            data frame / dict-like structure containing named columns

        expr: str
            expression to evaluate on the table
            includes mathematical operations and attribute names

        exprvars: dictionary, optional
            A dictionary that replaces the local operands in current frame.

        dtype: dtype definition
            dtype of the output array

        Returns
        -------
        out : np.array
            array of the result
        """
        return evalexpr(self, expr, exprvars=exprvars, dtype=dtype)


def _df_multigroupby(ary, *args):
    """
    Generate nested df based on multiple grouping keys

    Parameters
    ----------
    ary: dataFrame, dict like structure

    args: str or sequence
        column(s) to index the DF
    """
    if len(args) <= 0:
        yield ary
    elif len(args) > 1:
        nested = True
    else:
        nested = False

    val = ary[args[0]]
    ind = sorted(zip(val, range(len(val))), key=lambda x:x[0])

    for k, grp in itertools.groupby(ind, lambda x:x[0]):
        index = [v[1] for v in grp]
        d = ary.__class__({a: np.array([b[i] for i in index])
                           for a, b in ary.items()})
        if nested:
            yield k, _df_multigroupby(d, *args[1:])
        else:
            yield k, d


def _df_multigroupby_aggregate(pv, func=lambda x:x):
    """
    Generate a flattened structure from multigroupby result

    Parameters
    ----------
    pv: dataFrame, dict like structure
    result from :func:`_df_multigroupby`

    func: callable
        reduce the data according to this function (default: identity)

    Returns
    -------
    seq: sequence
        flattened sequence of keys and value
    """
    def aggregate(a, b=()):
        data = []
        for k, v in a:
            if type(v) in (list, tuple,):
                data.append(aggregate(v, b=(k,)))
            else:
                data.append(b + (k, func(v)))
        return data
    return list(itertools.chain(*aggregate(pv)))


def evalexpr(data, expr, exprvars=None, dtype=float):
    """ evaluate expression based on the data and external variables
        all np function can be used (log, exp, pi...)

    Parameters
    ----------
    data: dict or dict-like structure
        data frame / dict-like structure containing named columns

    expr: str
        expression to evaluate on the table
        includes mathematical operations and attribute names

    exprvars: dictionary, optional
        A dictionary that replaces the local operands in current frame.

    dtype: dtype definition
        dtype of the output array

    Returns
    -------
    out : np.array
        array of the result
    """
    _globals = {}
    keys = []
    if hasattr(data, 'keys'):
        keys += list(data.keys())
    if hasattr(getattr(data, 'dtype', None), 'names'):
        keys += list(data.dtype.names)
    if hasattr(data, '_aliases'):
        # SimpleTable specials
        keys += list(data._aliases.keys())
    keys = set(keys)
    if expr in keys:
        return data[expr]
    for k in keys:
        if k in expr:
            _globals[k] = data[k]

    if exprvars is not None:
        if (not (hasattr(exprvars, 'items'))):
            raise AttributeError("Expecting a dictionary-like as condvars with an `items` method")
        for k, v in ( exprvars.items() ):
            _globals[k] = v

    # evaluate expression, to obtain the final filter
    # r = np.empty( self.nrows, dtype=dtype)
    r = eval(expr, _globals, np.__dict__)

    return np.array(r, dtype=dtype)
