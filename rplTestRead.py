import matplotlib as mpl
import hyperspy.api as hs
import numpy as np
%matplotlib inline 

s = hs.load("CuAg850_SOI4.rpl",reader="rpl",signal_type = "EDS_TEM").T
s.axes_manager.signal_axes[0].units = 'keV'

s.set_elements(['Ag','Cu'])
s.metadata.Sample


s.get_lines_intensity(['Ag_La'], plot_result=True)
#im = s.get_lines_intensity()
#hs.plot.plot_images(im, tight_layout=True, cmap='RdYlBu_r', axes_decor='off',
#   colorbar='single', vmin='1th', vmax='99th', scalebar='all',
#    scalebar_color='black', suptitle_fontsize=16,
#    padding={'top':0.8, 'bottom':0.10, 'left':0.05,
#             'right':0.85, 'wspace':0.20, 'hspace':0.10})