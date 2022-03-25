import numpy as np

Nx = input('Nx ? ')
Nx = int(Nx)
Nz = input('Nz ? ')
Nz = int(Nz)
Lx = input('Lx ? ')
Lx = float(Lx)
Lz = input('Lz ? ')
Lz = float(Lz)

Re = input('Re ? ')
Re = float(Re)

kx = 2*np.pi*Nx/Lx
kz = 2*np.pi*Nz/Lz

kcut = max(kx,kz)

freq = kcut**2 / Re

dt  = 1/freq/2 #safety factor

print("dt = %s" % dt)
