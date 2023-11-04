# -----------------------------------------------------------------------
#  Program Algoritma Genetika
#  Refererensi Won Y. Yang dkk
#  Difasilitasi oleh PT. Diamond Startup Indonesia
#  Dibuat Oleh Asrin dan Sidratul Akbar
#  Kendari 2022
# -----------------------------------------------------------------------

import time

import matplotlib.pyplot as plt
import numpy as np
import numpy.matlib
import pandas as pd
from numpy.random import rand

# Menghitung waktu running program
start_time = time.time()
# Input Parameter Model
data = pd.read_excel('Square_Model.xlsx', sheet_name='Sheet1')
grav = np.array(np.genfromtxt('square model.dat'))

# Anomali Gravitasi Data sintetik
f = np.array(grav[:, 1])
f_scpace = np.array(grav[:, 0])

# Parameter Model Blok 2D
nlay = int(data.loc[0, 'nb'])       # Banyak grid
nxg = int(data.loc[0, 'nx'])        # Banyak Blok lateral
nzg = int(data.loc[0, 'nz'])        # Banyak Blok vetikal
dx = int(data.loc[0, 'dx'])         # Dimensi Blok lateral (m)
dh = int(data.loc[0, 'dz'])         # Dimensi Blok vertikal (m)
model = np.array([nxg, nzg, dx, dh])
l = data.loc[1, :]              
l = np.array([l[0:int(nlay)]])
l = l[0]                            # Batas bawah pencarian
u = data.loc[2, :]
u = np.array([u[0:int(nlay)]])
u = u[0]                            # Batas atas pencarian
x0 = l+(u-l)/2                      # Tebakan model awal
Np = 30                             # Jumlah populasi
Nb = data.loc[3, :]                 
Nb = np.array(Nb[0:int(nlay)])      # Jumlah bit
Pc = 1                              # Probabilitas crossover
Pm = 0.01                           # Probabilitas mutasi
eta = 0.1                           # Learning rate and the maximum
kmax = 1000                         # Jumlah iterasi

# ================== Fungsi Algoritma genetika ================== # 
def genetic(f, x0, l, u, Np, Nb, Pc, Pm, eta, kmax, model):
    gmo = forward_grav(x0, model)
    fo = 1/len(gmo)*(np.sum(np.sqrt((gmo-f)**2)))  # Misfit
    X = x0
    for i in range(1, Np):
        x = l+rand(len(x0))*(u-l)
        X = np.append(X, x)
    X = X.reshape(Np, len(x0))
    P = gen_encode(X, Nb, l, u)

    kk = 0
    mfp = []
    while kk < kmax:
        fX = []  # mf
        Xdec = gen_decode(P, Nb, l, u)
        n = 0
        while n < Np:
            gm = forward_grav(Xdec[n, :], model)  # fungsi
            mf = 1/len(gm)*(sum(np.sqrt((gm-f)**2)))  # Misfit
            fX.append(mf)
            n = n+1
        # print('ini ',kk)
        fxb = min(fX)
        mfp.append(fxb)
        nb = fX.index(fxb)
        if fxb < fo:
            fo = fxb
            X = Xdec[nb, :]
        fX1 = max(fX)-fX
        fXm = fX1[nb]
        if fXm < np.finfo(float).eps:
            print('stop')
            return X, mfp  # terminate
        for n in range(Np):
            Xdec[n, :] = Xdec[n, :]+eta * \
                (fXm-fX1[n])/fXm*(Xdec[nb, :]-Xdec[n, :])
        P = gen_encode(Xdec, Nb, l, u)
        iss = shuffle(np.array(range(Np)))
        for n in range(0, Np, 2):
            if np.random.rand(1) < Pc:
                P[iss[n:n+2], :] = crossover(P[iss[n:n+2], :], Nb)

        P = mutation(P, Nb, Pm)
        kk = kk+1
    return X, fo, mfp

def forward_grav(rho, model):
    G = 6.6732e-11   # Kontanta Gravirtasi
    dh = model[2]/2
    VV = np.reshape(rho, (model[1], model[0])).transpose().ravel()

    x = np.arange(0, (model[0]-1)*model[3]+model[3], model[3])+dh
    z = np.array([np.arange(0, ((model[1]+1)*model[2]-10), model[2])])+dh
    nb = model[0]*model[1]

    xxx = np.array(numpy.matlib.repmat(x, model[1], 1)).transpose().ravel()
    zzz = np.array(numpy.matlib.repmat(z, 1, model[0])).ravel()

    AA = []
    for i in range(model[0]):
        for j in range(nb):
            a = 2*G*((x[i]-xxx[j]+dh)*np.log(((((zzz[j]+dh)**2 + (x[i]-xxx[j]+dh)**2)**0.5)*(((zzz[j]-dh)**2 + (x[i]-xxx[j]-dh)**2)**0.5))/((((zzz[j]-dh)**2 + (x[i]-xxx[j]+dh)**2)**0.5)*(((zzz[j]+dh)**2 + (x[i]-xxx[j]-dh)**2)**0.5)))
                     + model[2]*np.log((((zzz[j]+dh)**2 + (x[i]-xxx[j]-dh)**2)**0.5)/(((zzz[j]-dh)**2 + (x[i]-xxx[j]-dh)**2)**0.5)) - (
                         zzz[j]+dh)*((np.arctan((x[i]-xxx[j]-dh)/(zzz[j]+dh)))-(np.arctan((x[i]-xxx[j]+dh)/(zzz[j]+dh))))
                     + (zzz[j]-dh)*((np.arctan((x[i]-xxx[j]-dh)/(zzz[j]-dh)))-(np.arctan((x[i]-xxx[j]+dh)/(zzz[j]-dh)))))
            AA.append(a)
    return np.dot(np.array(AA).reshape(model[0], nb), VV)*1e5

def gen_encode(X, Nb, l, u):
    P = []
    for n in range(np.size(X, 0)):
        for m in range(len(Nb)):
            p = np.binary_repr(
                round(int(((2**Nb[m])-1)*(X[n, m]-l[m])/(u[m]-l[m]))), width=int(Nb[m]))
            pp = list(''.join([(elem)for elem in p]))
            P.append(pp)

    return np.array(P).reshape(Np, len(Nb)*10)

def gen_decode(P, Nb, l, u):
    X = []
    for n in range(np.size(P, 0)):
        b2 = 0
        for m in range(len(Nb)):
            b1 = b2
            b2 = b1+Nb[m]
            X.append(np.array([binaryToDecimal(P[n, int(b1):int(b2)]) * (u[m]-l[m])/(((2**(Nb[m]))-1))+l[m]
                               ]))
    return np.array(X).reshape(Np, len(Nb))

def binaryToDecimal(n):
    n = ''.join([(elem)for elem in n])
    return int(n, 2)

def shuffle(iss):
    N = len(iss)-1
    for n in range(N, 1, -1):
        par = rand(1)*(n-1)
        inn = np.ceil(par)
        tmp = iss[int(inn)]
        iss[int(inn)] = iss[n]
        iss[n] = tmp  # swap the n-th element with the in-th one
    return iss

def crossover(chrms2, Nb):
    b2 = 0
    for m in range(len(Nb)):
        b1 = b2
        bi = b1+np.mod(np.floor(rand(1)*Nb[m]), Nb[m])-1
        b2 = b2+Nb[m]
        tmp = chrms2[0,int(bi):int(b2)]
        chrms2[0,int(bi):int(b2)] = chrms2[1,int(bi):int(b2)]
        chrms2[1,int(bi):int(b2)] = tmp
    return chrms2

def mutation(P, Nb, Pm):
    for n in range(np.size(P, 0)):
        b2 = 0
        for m in range(len(Nb)):
            b1 = b2
            bi = b1+np.mod(np.floor(rand(1)*Nb[m]), Nb[m])
            b2 = b2+Nb[m]
            if rand(1) < Pm:
                P[int(n),int(bi)] = str(1-eval(P[int(n),int(bi)]))
    return P


xo_gen, fo_gen,mfp = genetic(f, x0, l, u, Np, Nb, Pc, Pm, eta, kmax, model) # Algoritma Genetika
gm_cal = forward_grav(xo_gen, model)  
VV = np.array(xo_gen).reshape(nzg, nxg)
zSA1 = np.array(range(0, 40, 10)).transpose()
zSA = zSA1+10/2

# Visualisasi model penampang inversi GA 2D
plt.figure('Inversi GA Gravitasi',figsize=(9,7.5))
plt.subplot(2,1,1).set(xlim=(0,90),ylim=(0,0.4))
plt.plot(f_scpace,f,'-ob',label = 'G-obs')
plt.plot(f_scpace,gm_cal,'-or',label = 'Inversi GA')
plt.legend()
plt.title('Penampang Hasil Inversi GA 2D Gravitasi')
plt.xlabel('Spasi [M]')
plt.ylabel('Anomali Medan Gravitasi [mGal]')
plt.subplot(2,1,2)
a = plt.imshow(VV,  extent = [0, 9, 5, 0], aspect = 'auto')
plt.ylabel('Kedalaman [M]')
plt.colorbar(a, orientation='horizontal')
plt.clim(0,1200)

# Visualisasi misfit rata-rata
plt.figure('Kurva Misfit Rata-rata',figsize=(9,7.5))
plt.subplot(1,1,1).set(ylim=(0,0.4))
plt.plot(mfp,'r-')
plt.title('Kurva Misfit Rata-rata')
plt.xlabel('Iterasi (Generasi)')
plt.ylabel('Misfit Rata-rata [mGal]')

# Output
print('==============================Model===========================\n',VV)
print('=============================Misfit==========================\n',fo_gen)
print("========================== %s seconds =========================" % (time.time() - start_time))
plt.show()
