import os

# should manually change it in the mpl config file ... this is stupid.
if 'Linux' in os.uname():
    import matplotlib as mpl
    mpl.use('Agg')
