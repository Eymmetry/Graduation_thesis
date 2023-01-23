import matplotlib.pyplot as plt
#import matplotlib as mpl
import numpy as np
import math
import pickle
import os

rptpypt=0.18603296098519714
rptpy=0.023398251556292975
rmnpy=0.07488608297935191
rptpymn=0.08390254537713687
rpy=0.18270118603207425

lptpypt=0.04386752263479449
lptpy=0.02890588787224217
lmnpy=0.070827733631309
lptpymn=0.09136695321952298
lpy=0.14319174613787836

mptpypt=(rptpypt+lptpypt)/2
mptpy=(rptpy+lptpy)/2
mmnpy=(rmnpy+lmnpy)/2
mptpymn=(rptpymn+lptpymn)/2
mpy=(rpy+lpy)/2

sample=np.array([1,2,3,4,5])
right=np.array([rptpypt,rptpy,rmnpy,rptpymn,rpy])
left=np.array([lptpypt,lptpy,lmnpy,lptpymn,lpy])
mid=np.array([mptpypt,mptpy,mmnpy,mptpymn,mpy])

plt.scatter(sample,right,label="Right Peak")
plt.scatter(sample,left,label="Left Peak")
plt.scatter(sample,mid,label = "meanvalue")
plt.legend()
plt.autoscale()

plt.xticks([1,2,3,4,5])
plt.xlabel("sample No.", fontsize=20)
plt.ylabel('$\Delta \mathit{P}^{norm} \\times 10^{4}$ ',fontsize = 20)

plt.savefig("p-sample")
plt.show()
plt.close()