using Plots

using FFTW

x = [0,im,0,1,0,-im]
print(ifft(x))

