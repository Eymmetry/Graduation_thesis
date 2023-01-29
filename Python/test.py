import matplotlib.pyplot as plt
import numpy as np
import math
import pickle
np.fft.ifft([0,1j,0,1,0,-1j])

real=7.945427933E-5
imag=2.120862628E-4
e=np.sqrt((real)**2+(imag)**2)
loge=20*np.log10(e)
print("e=",e)
print("loge=",loge)