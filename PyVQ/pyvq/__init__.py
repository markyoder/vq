import pylab as plt
import matplotlib as mpl
#mpl.use('agg')

from pyvq import *
import betas

__all__=['pyvq', 'pyvq.betas']

if __name__=='__main__':
	# being loaded as root module (like from command line).
	mpl.use('agg')	# not sure if this is the right way/place to do this. in any case, we want to distinguisn when pyvq is being loaded from command-line vs interactive.
else:
	plt.ion()
	
