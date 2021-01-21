#Test av att plotta spektra med HyperSpy, det spektrat som användes var genererat från DSTA-II mha deras exempelscript. 
import matplotlib as mpl
import hyperspy.api as hs


s = hs.load("testspectraAu.msa",reader="msa",signal_type = "EDS_TEM")
s.plot()
