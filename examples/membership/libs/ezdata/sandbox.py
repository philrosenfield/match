""" Anything that is not in the mainstream yet """

def crosstab(d, *args, **kwargs):
    """
    Compute a simple cross-tabulation of two (or more) factors. By default
    computes a frequency table of the factors unless an array of values and an
    aggregation function are passed

    Parameters
    ----------
    d: dict like structure
        data container (SimpleTable or DictDataFrame)
    index : array-like, Series, or list of arrays/Series
        Values to group by in the rows
    columns : array-like, Series, or list of arrays/Series
        Values to group by in the columns
    values : array-like, optional
        Array of values to aggregate according to the factors
    aggfunc : function, optional
        If no values array is passed, computes a frequency table
    rownames : sequence, default None
        If passed, must match number of row arrays passed
    colnames : sequence, default None
        If passed, must match number of column arrays passed
    margins : boolean, default False
        Add row/column margins (subtotals)
    dropna : boolean, default True
        Do not include columns whose entries are all NaN

    Returns
    -------
    crosstab : DataFrame
        result
    """


import itertools
import numpy as np

import scipy.sparse as sps


def pivottab(d, rowkey, columnkey, reducefn, duplicates=False):
    """ A simple pivot table generator based on 2 key indexes and a reduce function

    Parameters
    ----------
    d: dict like structure
        data container (SimpleTable or DictDataFrame)
    rowkey : str
        Values to group by in the rows
    columns : str
        Values to group by in the columns
    reducefn: callable or str
        reduce function used in aggregate
        if `str` argument, use itemgetter
    duplicate: bool
        in case of repeated entries for a same row/column combination, values
        will be added. (should not append after :func:`aggregate`)

    Returns
    -------
    piv: ndarray
        2d array with value(row, column)
    row: ndarray
        1d array of the row index values
    column: ndarray
        1d array of the column index values
    """
    if hasattr(d, 'aggregate'):
        agg = d.aggregate(reducefn, (rowkey, columnkey))
    else:
        agg = aggregate(d, reducefn, (rowkey, columnkey))
    data = np.array(list(agg))
    rows, row_pos = np.unique(data[:, 0], return_inverse=True)
    cols, col_pos = np.unique(data[:, 1], return_inverse=True)

    if duplicates:
        pivot_table = sps.coo_matrix((data[:, 2], (row_pos, col_pos)),
                                     shape=(len(rows), len(cols))).A
    else:
        pivot_table = np.zeros((len(rows), len(cols)), dtype=data.dtype)
        pivot_table[row_pos, col_pos] = data[:, 2]

    return pivot_table, rows, cols


def aggregate(d, func, keys, args=(), kwargs={}):
    """ Apply func on groups within the data

    Parameters
    ----------
    d: dict like structure
        data container
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
    pv = [(k, list(v)) for k, v in _df_multigroupby(d, *keys)]

    def _aggregate(a, b=()):
        data = []
        for k, v in a:
            if type(v) in (list, tuple,):
                data.append(_aggregate(v, b=(k,)))
            else:
                data.append(b + (k, func(v)))
        return data

    return list(itertools.chain(*_aggregate(pv)))


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
