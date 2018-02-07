import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.optimize import curve_fit

np.set_printoptions(formatter={'float':'{: .5E}'.format})
Aex = 13e-12
mu0= 4 * np.pi * 1e-7
ms = 1./mu0
alpha = 0.008
V = 2e-9**3
T = 1.
hext = 1.0/mu0
gamma=221276.1488637255
kb=1.38064852e-23
w0 = gamma*hext

def read_m_t_magnumfe(logfile):
  ''' return (t, mx, my, mz) read from magnumfe logfile'''
  x = np.loadtxt(logfile)
  return x[:,(0,1,2,3)]

def fft(t, x, dt):
  # interpolate on equidistant grid
  x_func = interpolate.interp1d(t, x)
  t_interp = np.arange(min(t), max(t), dt)
  x_interp = x_func(t_interp)

  # perform FFT
  XX = np.fft.fftshift(np.fft.fft(x_interp))
  f = np.fft.fftshift(np.fft.fftfreq(len(t_interp), d=dt))
  #X = np.abs(XX)**2 # fangohrs version
  X = (2*np.absolute(XX))**2*dt**2/(max(t)-min(t)) # dieters version
  return f, X

# main routine
x = read_m_t_magnumfe("m.dat")
t, my = x[:,0], x[:,2]
f, MY = fft(t, my, t[1]-t[0])

np.savetxt("psd.dat", np.vstack([f,MY]).T, header="Frequency[Hz], Spectral Density")

plt.subplot(2, 1, 1)
plt.grid()
plt.plot(t,my)
plt.xlabel('Time [s]')
plt.ylabel('Polarization Jx [T]')

plt.subplot(2, 1, 2)
#plt.plot(f,MY)
plt.xlim([0,3*w0/2/np.pi])
plt.ylim([1e-18, 1e-10])
plt.semilogy(f, MY, 'bo', label="simulation")
plt.semilogy(f, MY, 'b-', label=None, alpha=0.2)

func = lambda w,w0,dw,a: a*(1/((w-w0)**2+dw**2) + 1/((w+w0)**2+dw**2))
p0=np.array([w0, alpha*w0, gamma*alpha*kb*T/mu0/ms/V])
popt, pcov = curve_fit(func, 2*np.pi*f, MY, p0=p0)
print "popt:", popt
print "p0:  ", p0
f = np.linspace(0, 3*w0/2/np.pi, num=1000)
plt.semilogy(f, func(2.*np.pi*f, *popt), 'b-', linewidth=2, alpha=0.2)

PSD_analytic = lambda w: gamma*alpha*kb*T/mu0/ms/V/((w-w0)**2+(alpha*w0)**2) + gamma*alpha*kb*T/mu0/ms/V/((w+w0)**2+(alpha*w0)**2)
plt.semilogy(f, [PSD_analytic(2.*np.pi*fi) for fi in f], "g-", label="analytic", linewidth=2, alpha=0.5)
plt.legend()
plt.grid()
plt.xlabel('Frequency [Hz]')
plt.ylabel('FFT of Polarization F{Jy} [T]')
plt.savefig("fft.png")
#plt.show()
