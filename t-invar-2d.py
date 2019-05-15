
from mpl_toolkits import mplot3d

import numpy as np
import math as m    
import matplotlib.pyplot as plt


## 2D Fourier field for each point in some grid

def fourier_field(kx,ky,box):
    size = int(m.sqrt(np.size(box)))
    size_vec = np.size(kx)
    l = int(m.sqrt(size_vec))
    psi = np.zeros((l,l))
    if size_vec >= 2:
        for i in range(0,l):
            for j in range(0,l):
                if int(kx[i,j]) >= l:
                    psi[i,j] = 0
                elif int(ky[i,j]) >= l:
                    psi[i,j] = 0
                elif int(kx[i,j]) < -l:
                    psi[i,j] = 0
                elif int(ky[i,j]) < -l:
                    psi[i,j] = 0
                else:
                    psi[i,j] = box[int(kx[i,j]),int(ky[i,j])]
        print(psi)
        return psi
    else:
        if kx >= size:
            return 0
        elif ky >= size:
            return 0
        elif kx < -size:
            return 0
        elif ky < -size:
            return 0
        return box[kx,ky]
    
## Transform a non ordered 2D matrix and relabel the points in such a way
## that is like a xOy plot. The point (0,0) will represent the center element
## of the "box" input

def grid_box(box):
    size = int(m.sqrt(np.size(box)))
    hsize = size//2
    grid = np.zeros((size,size))
    tbox = box.transpose()
    box1 = np.zeros((size,size))
    
    for i in range(0,size):
        for j in range(0,size):
            if j <= hsize:
                box1[i,j] = tbox[i,hsize-j]
            else:
                box1[i,j] = tbox[i,j]
                
    
    for k in range(0,size):
        for l in range(0,size):
            if l == size-1:
                grid[l,k] = box1[0,k]
            else:
                grid[l,k] = box1[l+1,k]
                
    return grid

## Makes N odd_number by odd_number random matrices.


def NBox(N,odd_number):
    box = np.random.rand(N,odd_number,odd_number)
    return box

## Computes the A(k,q) estimator for certain vectors k = (x1,y1)
## and q = (x2,y2)

def A_estimator(x1,x2,y1,y2,Nbox):
    oddnumber = int(m.sqrt(np.size(Nbox[0])))
    N = np.size(Nbox)//(oddnumber**2)
    
    grid = np.zeros((N,oddnumber,oddnumber))
    for i in range(0,N):
        grid[i] = grid_box(Nbox[i])
    
    size_vec = np.size(x1)
    l = int(m.sqrt(size_vec))

    if size_vec >= 2:
        psi1 = np.zeros((N,l,l))
        psi2 = np.zeros((N,l,l))
        s = 0
        for k in range(0,N):
            psi1[k] = fourier_field(x1-x2,y1-y2,grid[k])
            psi2[k] = fourier_field(x1+x2,y1+y2,grid[k])
            for i in range(0,l):
                for j in range(0,l):
                    s = s + psi1[k,i,j]*psi2[k,i,j]

        A = s/N
        return A

    else:
                                

        psi1 = np.zeros(N)
        psi2 = np.zeros(N)

        for k in range(0,N):
            psi1[k] = fourier_field(x1-x2,y1-y2,grid[k])
            psi2[k] = fourier_field(x1+x2,y1+y2,grid[k])
            print(k)

        s = 0

        for l in range(0,N):
            s = s + psi1[l]*psi2[l]

        A = s/N

        return A

## Compute how much translation invariant the field is in some vector
## k = (x1,y1) direction


def Sigma_estimator(x1,y1,Nbox):
    oddnumber = int(m.sqrt(np.size(Nbox[0])))
    N = np.size(Nbox)//(oddnumber**2)
    s = 0

    for i in range(-oddnumber+1,oddnumber):
        for j in range(-oddnumber+1,oddnumber):
            s = s + (A_estimator(x1,i,y1,j,Nbox))**2

    return s

##### Plots still in development #####

def TI_measu_plot(Nbox):
    oddnumber = int(m.sqrt(np.size(Nbox[0])))
    length = int(oddnumber/2 - 1/2)
    N = np.size(Nbox)//(oddnumber**2)

    kx = np.linspace(-length,length,oddnumber)
    ky = np.linspace(-length,length,oddnumber)

    X, Y = np.meshgrid(kx, ky)
    Z = Sigma_estimator(X,Y,Nbox)

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.contour3D(X, Y, Z, 50, cmap='binary')
    ax.set_xlabel('kx')
    ax.set_ylabel('ky')
    ax.set_zlabel('Sigma');

def test(x,y):
    return x**2+y**2
    
