# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 15:29:51 2022

@author: shinozaki ryo
"""
import matplotlib.pyplot as plt
import numpy as np
import math
import pickle
#-----tracedata---------------------------------------------------------
#データの読み込み
switch_trace = False #事前にメモリしたtracedataとのベクトル差を作りたいときはこっちをTrueにする
#ただし、ifft前に周波数信号のジャンプを取り除くために窓関数を掛けたりしていないし、折り返しのやり方も下と異なる不完全な部分なのでココを使いたくなったら次のセクションを参考に直して下さい。
if switch_trace:
    Real = np.genfromtxt("Real_mem.txt") #重要！このデータは電圧比を表していて無次元です。同時にVNAがはき出すLogMagというデータはもちろんdB単位になっています。
    Imag = np.genfromtxt("Img_mem.txt")
    freq = np.genfromtxt("freq0.txt")
    
    #Sパラメータと絶対値を計算
    Strace = Real + Imag*1j #Sパラメータ
    Magtrace = abs(Strace) #絶対値
    
    #フーリエ変換のために周波数を0GHzからにして折り返す。折り返すのは精度を確保するため。折り返した部分は複素共役にする
    #参考：https://mashiroyuya.hatenablog.com/entry/ifftpython
    df = freq[1]-freq[0]
    offset = np.floor(freq[0]/df)
    F0 = []
    n = int(2*len(freq) + 2*offset + 1)
    for i in range(n): #周波数帯を0GHｚからにして折り返した分と中央の0
        if i < len(freq):
            F0.append(Strace[len(freq)- i - 1].conjugate())
        elif i > int(len(freq) + 2*offset):
            F0.append(Strace[i-(len(freq) + 2*int(offset) + 1)])
        else:
            F0.append(0)
    F = np.array(F0) #Numpy配列化

    #ifftでtimedomainにしたときの横軸容易
    #t = 2*np.arange(2*len(freq) + 2*int(offset) + 1)/(df*2.0*(2*len(freq) + 2*offset + 1))
    t = np.arange(2*len(freq) + 2*int(offset) + 1)

#ifft
    f0 = np.fft.ifft(F)

    #最初に到達した表面弾性波だけで切り取り
    for i in range(n):
        if (i > 323) and (i < 973):
            f0[i] = f0[i]
        else:
            f0[i] = 0

#fft
    fafter = np.fft.fft(f0)
    Strace_reduction = fafter[n-len(freq):n]
    Magtrace_reduction = abs(Strace_reduction)
#-----tracedata------------------------------------------------------------
Real = np.genfromtxt("Real_data0_0.txt")
#-----Spara---------
switch_DeltaP = True
if switch_DeltaP:#switch
    Pin=10**(-5/10)*1e-3 #入射電力は-5dBmに設定。ここでは単位をWに直す。
    #角度依存性で切る
    for i in range(1): #0°～360°まで37測定分の解析
        Real = np.genfromtxt("Real_data0_0.txt")
        Imag = np.genfromtxt("Imag_data0_0.txt")
        field = np.genfromtxt("gaussX_0_0.txt")
        freq = np.genfromtxt("freq0.txt")

        #変数の設定
        Nfreq=len(freq)
        Nfield=len(field)
        N_Sreal=len(Real[:,0])
        N_Simag=len(Imag[:,0])

        freqmin=freq[0]#最小周波数
        freqmax=freq[Nfreq-1]#最大周波数
        Deltaf=(freqmax-freqmin)/(Nfreq-1)

        #説明しよう！mp.linespaceとは等差数列をはき出す関数である
        #ここではifftを施すために、周波数配列の区間を0からにし、後でifft後の精度を確保するための折り返しを行う準備をする
        #折り返すといっても、0GHzのところで折り返すわけではない
        f0=np.linspace(0, freqmin, num = int(freqmin/Deltaf), endpoint = True, retstep = False, dtype = None)
        f=np.append(f0,freq)
        fturn=np.linspace(freqmax, freqmax*2, num = len(f), endpoint = True, retstep = False, dtype = None)
        f2=np.append(f,fturn)

        Dt=1/np.amax(f2)#1/4GHz  4GHzはサンプリング周波数
        ifT=range(0,len(f2))*Dt

#秒数でゲーティング範囲指定したいときはこっち
        gst=120e-9#gating_start
        ged=215e-9#gating_end
        X=int(gst/Dt)#gaint_startのポイント数
        Y=int(ged/Dt)#gaint_endのポイント数
#        X=2000
#        Y=5000


        zero1=np.zeros(X)#時間領域において0~X番目まで0
        zero2=np.zeros(len(f2)-Y)#Y~6402まで0


        #----------------------------------------
        #ゼロと空の配列の準備
        #----------------------------------------
        num=0
        Real2=np.empty(N_Sreal)
        dfSarray1=np.empty(N_Sreal)
        dfSarray2=np.empty(N_Sreal)

#磁場の強度で切る
        for j in range(len(field)):
            Sreal =Real[:,j]
            Simag =Imag[:,j]
            Real2=np.vstack((Real2,Sreal)) #ndarrayの結合
      
            #窓関数を掛ける。
            #このままだと0~fmin GHzとfmin~fmax GHzの境でSパラメータが0から有限の値にジャンプする。
            #そのジャンプがフーリエ変換する上でノイズになるので窓関数を掛ける
            Sreal = Sreal*(0.54-0.46*np.cos(2*np.pi*(freq-freqmin)/(freqmax-freqmin)))
            Simag = Simag*(0.54-0.46*np.cos(2*np.pi*(freq-freqmin)/(freqmax-freqmin)))
            
            #周波数範囲を拡張する
            S0=np.zeros(len(f)-len(Sreal))
            SrealzeroPD=np.append(S0,Sreal)
            SimagzeroPD=np.append(S0,Simag)
            
            Ssum=SrealzeroPD+1j*SimagzeroPD
            
            #-----------------------------------------------------------------------
            #折り返し
            #-------------------------------------------------------------------------
            #詳しくはわからんのだけどfftしてできるこの信号を作ってるっぽい。だから折り返しがf=0GHz対称のものでなくfmax対称としている
            #http://exp1gw.ec.t.kanazawa-u.ac.jp/DSP/Signal-Processing/frequency-domain.html
            Sturn=np.conjugate(Ssum[::-1]) #[::a]で要素をa個毎に切り出す、今回は-1だから折り返しになる
            Ssum=np.append(Ssum,Sturn)
            
            #--------------------------------------------------
            #逆フーリエ変換
            #--------------------------------------------------
            ifS=np.fft.ifft(Ssum)
                
            #--------------------------------------------------
            #gating
            #--------------------------------------------------
            ifSgat=ifS[X:Y]

            ifSgat=np.append(zero1,ifSgat)
            ifSgat=np.append(ifSgat,zero2)#時間「領域においてゲーティングしたSpara
 
            #-------------------------------------------------------
            #フーリエ変換
            #-------------------------------------------------------             
            dfS=np.fft.fft(ifSgat) 
            
            #----------------------------------------------------------
            #前半,ゼロを抽出
            #-----------------------------------------------------------
            dfS=dfS[len(f0):len(f0)+len(freq)]
            
            if True:
                kari1=20*np.log10(abs(ifS))
                kari2=20*np.log10(abs(ifSgat))
#                plt.plot(range(0,len(f2)),kari1)
#                plt.plot(range(0,len(f2)),kari2)
                plt.plot(ifT*1e6,kari1)
                plt.plot(ifT*1e6,kari2)
#                plt.xlim(2300,3600)
                plt.autoscale()
            
                plt.title("time domain")
                plt.xlabel("t / ${\mu}s$ ", fontsize=10)
#                plt.xlabel("points", fontsize=10)
                plt.ylabel("$S_{21}$/dB", fontsize=10)
                plt.xticks(fontsize=10)
                plt.yticks(fontsize=10)
                plt.rcParams["figure.figsize"] = (10, 6)
                plt.show()
                
                dfskari=Pin*abs(dfS)**2
                plt.plot(freq*1e-9,10*np.log10(abs(dfS)))
                plt.rcParams["figure.figsize"] = (10, 6)
                plt.xlim([0.5,10])
                plt.title("frequency domain")
                plt.xticks(fontsize=10)
                plt.yticks(fontsize=10)
                plt.xlabel("$\mathit{f}$ /GHz",fontsize=10)
                plt.ylabel("$S_{21}$/dB", fontsize=10)
                plt.show()
                p_toutatsu=2320 #波が検知用IDTに到達した時間の点数
                velocity = 400e-6/(X*Dt - 61e-6/3940)
                print(velocity)
                stop
      
      
            #---------------------------------------------------------------------
            #窓関数を戻す
            #---------------------------------------------------------------------      
            dfS=dfS/(0.54-0.46*np.cos(2*np.pi*(freq-freqmin)/(freqmax-freqmin)))#gating後の周波数領域Spara
            dfSarray2=np.vstack((dfSarray2,dfS))
            dfSarray1=np.vstack((dfSarray1,dfSarray2[j+1]-dfSarray2[1]))
            
        #--------------------------ゲーティング処理後のSpara!!!-----------------
        dfSarray2= np.delete(dfSarray2, 0, axis=0)#gating処理後のすべてのデータ　0番目のデータはnp.emptyで作った行で0いらないから消す。
        dfSarray1= np.delete(dfSarray1, 0, axis=0)

        #---------------------------------------------------------------
        dfSarray2_real= np.delete(dfSarray2.real, 0, axis=0)#保存用に実部を取り出す 
        dfSarray2_imag= np.delete(dfSarray2.imag, 0, axis=0)#保存用に虚部を取り出す
        dfS_saw=dfSarray2[0,:]#磁場0番目、トレースデータ
        dfS_saw_power=Pin*abs(dfS_saw)**2
        saw_max=np.amax(dfS_saw_power)
        
        
        dfS_trace=dfSarray2-dfS_saw#S-S_saw、行が周波数、列が磁場
        P_SAW=Pin*abs(dfS_saw)**2#SAWエネルギーのS->Pに変換、配列
        P_SAWmax=np.amax(Pin*abs(dfS_saw)**2)#P_sawの最大、値
        DeltaPnorm = Pin*abs(dfSarray1)**2/saw_max
        tracedata = abs(dfSarray2[0,:])**2


#------------------保存
        np.savetxt("DeltaPnorm.txt", DeltaPnorm)
        np.savetxt("tracedata.txt", tracedata)
        
            
            
            
            
        



#plt.plot(freq,Sreduction)
#plt.show()

