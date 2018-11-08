
## Make sure plots are create inline
#%matplotlib inline
## The usual packages (numpy, matplotlib, etc)
import matplotlib.pyplot as plt
import numpy as np
# nicer figures using ggg plot style.
plt.style.use('ggplot')
from IPython.core.display import HTML
def css_styling():
    styles = open("../python/styles/custom.css", "r").read()
    return HTML(styles)
css_styling()
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=False)
from IPython.core.pylabtools import figsize
figsize(11,11/1.618)
