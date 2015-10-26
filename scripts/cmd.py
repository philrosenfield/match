"""Class for reading the output.cmd file from calcsfh"""
import os
import numpy as np

from .fileio import read_match_cmd

__all__ = ['CMD']

class CMD(object):
    """
    A quikly made object to read the MATCH CMD file and hold paramters to
    automatically make plots with the same color scale as other MATCH CMD files.
    """
    def __init__(self, filename):
        self.cmd = read_match_cmd(filename)
        self.figname = os.path.split(filename)[1] + '.png'
        labels = ['${\\rm %s}$' % i for i in ('data', 'model', 'diff', 'sig')]
        labels[1] = '${\\rm %s}$' % self.figname.split('.')[0].replace('_', '\ ')
        self.labels = labels
        self.load_match_cmd(filename)

    def load_match_cmd(self, filename):
        """
        pgcmd needs hesses and extent. Optional are max_* which set the vmins
        and vmaxs.
        """
        self.nmagbin = len(np.unique(self.cmd['mag']))
        self.ncolbin = len(np.unique(self.cmd['color']))
        self.data = self.cmd['Nobs'].reshape(self.nmagbin, self.ncolbin)
        self.model = self.cmd['Nsim'].reshape(self.nmagbin, self.ncolbin)
        self.diff = self.cmd['diff'].reshape(self.nmagbin, self.ncolbin)
        self.sig = self.cmd['sig'].reshape(self.nmagbin, self.ncolbin)
        self.hesses = [self.data, self.model, self.diff, self.sig]
        self.extent = [self.cmd['color'][0], self.cmd['color'][-1],
                       self.cmd['mag'][-1], self.cmd['mag'][0]]
        self.max_counts = np.nanmax(np.concatenate([self.data, self.model]))
        self.max_diff = np.nanmax(np.abs(self.diff))
        self.max_sig = np.nanmax(np.abs(self.sig))
