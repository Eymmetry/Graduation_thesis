# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 17:30:44 2015
program code for yanagisawa
produce color plot with subtracting reference signal

@author: Sho
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl



from matplotlib import pyplot as plt



params = {
    'font.family': 'Calibri', #使用するフォント名
    'mathtext.default': 'regular', #数式で通常使用するフォント．'regular'にすると数式のフォントが通常使用に設定されたフォントと同じになる
#    'font.weight': 'bold',#boldに出来る
    'font.size': 8, #全体のフォントサイズ．凡例のサイズ調整は何故かうまく動作せず．個々のフォントサイズを変えたい場合は以下を調整
    'axes.labelsize': 30, # 軸ラベルのフォントサイズ
    'legend.fontsize': 8, # 凡例の文字の大きさ
    'xtick.labelsize': 25, # x軸の数値の文字の大きさ
    'ytick.labelsize': 25, # y軸の数値の文字の大きさ
    'lines.linewidth'   : 2.0, 
    'figure.subplot.hspace' : 0.1 ,
    'figure.subplot.left'     : 0.2,
    'figure.subplot.bottom'     : 0.1,
    'figure.autolayout' : False,
    'axes.linewidth' : 1.0,
    'xtick.major.width'    : 1,
    'xtick.major.size'    : 2,
    'ytick.major.width'    : 1,
    'ytick.major.size'    : 2,
    'figure.figsize': [24.0/2.54, 24.0/2.54],
}
 
mpl.rcParams.update(params)
#mpl.rcParams['pdf.fonttype'] = 42

gyoave=1
retuave=1

load = True
if load:
    POWERdata= np.genfromtxt("DeltaPnorm.txt")
#    POWERdata= np.genfromtxt("@traceV5_typeB_normalized_sattelite.txt")
    tracedata= np.genfromtxt("tracedata.txt")
#    tracedata=POWERdata[0,:]
    freq= np.genfromtxt("freq0.txt")
    field= np.genfromtxt("gaussX_0_0.txt")
#    field= np.delete(field, 0, axis=0)

#   gauss=-0.0325223+0.189416*field-0.000494966*field**2+0.000536273*field**3+0.000418497*field**4-7.95163e-5*field**5-5.38487e-5*field**6+3.22656e-6*field**7+2.79955e-6*field**8-5.36458e-8*field**9-4.82089e-8*field**10
#　上の式はSAW用プローバーにおける磁場校正の式、下は知らん
#    gauss=gauss*100-1.63884
#    gauss=-3.1085+13.361*field-1.5152*field**2+0.7780*field**3+11.439*field**4-21.788*field**5+11.375*field**6
    

    freq=freq*1e-9
#    gauss=gauss*1e+2-1.169
    gauss=field-1.3
#定数部分はオフセット分（電磁石への入力電圧を0にしたときに残る磁場)

#gauss = gauss[::-1]
gauss.sort()


#from guiqwt.builder import make
#def imshow(win, x, y, data, yreverse=False, title=None):
#    image = make.xyimage(x, y, data, interpolation="nearest", title=title )
#    plot = win.get_plot()
#    plot.add_item(image)
#    plot.do_autoscale()
#    plot.replot()
    
#    return plot, image



#from guiqwt.plot import ImageDialog
#import guidata
#_app = guidata.qapplication()



plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.direction'] = 'in'
X, Y = np.meshgrid(freq,gauss)
Z = POWERdata*1e4
plt.subplots_adjust(wspace=100,hspace=100)
ax1=plt.subplot(1, 1, 1)
#plt.pcolor(X, Y, Z, cmap='rainbow')
plt.pcolor(X, Y, Z)
plt.colorbar(shrink=1,label='$\Delta \mathit{P}^{norm} \\times 10^{4}$ ')
#plt.xlabel("$\mathit{f}$ /GHz",fontsize=30)
plt.ylabel("$\mathit{\mu}_{0}\mathit{H}$ /mT",fontsize=30)
plt.ylim([-20,20])
plt.xlim([0.5,3])
plt.clim(0, 0.01) 
#plt.text(1.2,7,"750 nm Pt(40)",fontsize=50,color="white")
#750 nmだけなら１．５３５でよい
#plt.yticks([-6,-4,-2,0,2,4,6], ["-6","-4","-2","0","2","4","6"])
#plt.clim(vmin=0,vmax=40.549) 
#plt.setp(ax1.get_xticklabels(), visible=False)
plt.tick_params(labelbottom=False, labelleft=True)

#x = np.arange(-10 ,-4 , 0.05)
#plt.plot(x-x+1.604978125,x,color='cyan',linestyle="dashed")
#x = np.arange(-3 ,3.5 , 0.05)
#plt.plot(x-x+1.604978125,x,color='cyan',linestyle="dashed",label='SAW')
#x = np.arange(4.5 ,10 , 0.05)
#plt.plot(x-x+1.604978125,x,color='cyan',linestyle="dashed")


ramda=1.788411e-6
#ramda=1.3269e-6
M=0.98
gamma=2.203e5
mu = 4*np.pi/10**7
A=0
#A=1.6e-6
d=20e-9
Hk=0 #単位はmT
x = np.arange(0 , 3 , 0.05)


def B(ramda):
      return 2*A/(M*10**4/4/np.pi)*(2*np.pi/ramda)**2*10**-8
#B(ramda)は0.000469[T]
def p(ramda):
      return (1-np.exp(-d/ramda*2*np.pi))/(d/ramda*2*np.pi)  
#p(ramda)は0.9704
def g(x):
      return Hk+(-B(ramda)-M/2*p(ramda)+((M/2*p(ramda))**2+(2*np.pi*mu*x*10**9/gamma)**2)**0.5)*1e3
plt.plot(x,g(x),color='red',linestyle="dashed",label='MSBVW')
plt.plot(x,-g(x),color='red',linestyle="dashed")


#plt.legend(loc="upper right",fontsize=23)
#plt.text(1.563,8,"(a)NiFe/Cu",fontsize=40,color="white")


plt.axes([.200, .21, .56, .25])
plt.plot(freq,10*np.log10(tracedata),"-")
plt.xlim([0.5,3])
plt.ylim([-130,0])
plt.subplots_adjust(bottom=0.47) 
plt.xlabel("$\mathit{f}$ /GHz",fontsize=30)
plt.ylabel("$\mathit{S}_{21}$ /dB",fontsize=30)
#x = np.arange(-70 ,-10 , 0.05)
#plt.plot(x-x+1.604978125,x,color='cyan',linestyle="dashed")
#plt.xticks([0.01,2], ["0.01","2"])
plt.yticks([-20,-40,-60,-80,-100], ["-20","-40","-60","-80","-100"])
#plt.text(1.563,-60,"(b)",fontsize=40,color="black")

plt.savefig("@@@PowerImage.png", dpi=300)
plt.show()

