from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np
import pylandau
import math
import csv

# Create fake data with possion error
mpv, eta, sigma, A = 30, 5, 4, 1000
x = np.arange(0, 100, 0.5)
y = pylandau.langau(x, mpv, eta, sigma, A)
yerr = np.random.normal(np.zeros_like(x), np.sqrt(y))
yerr[y < 1] = 1
y += yerr


#read in tab delimited data file
energies=[]
counts=[]
countsErr = np.array([])
with open('decayEnergy_100kevBins_filter2.txt','r') as f:
    next(f) # skip headings
    reader=csv.reader(f,delimiter='\t')
    for energy,count in reader:
        energies.append(energy)
        counts.append(count)
       
       
E = np.array(energies, dtype=float)
C = np.array(counts, dtype=float)

#np.sqrt(E)


        
print(E)
# ('Mark', 'Matt', 'John', 'Jason', 'Matt', 'Frank', 'Frank', 'Frank', 'Frank')
print(C)
# ('32', '29', '67', '45', '12', '11', '34', '65', '78')
#print(countsErr)

# Create fake data with possion error
# parameters below seem to give correct distribution, btw the scaling factor (last in list) does abosolute scaling to value (unlike root langau)
mpv2, eta2, sigma2, A2 = 1240, 64, 45, 5
mpv3, eta3, sigma3, A3 = 1611, 64, 45, 2
x = np.arange(0, 10000, 0.5)
y = pylandau.langau(x, mpv2, eta2, sigma2, A2)
y2 = pylandau.langau(E, mpv2, eta2, sigma2, A2)
y5 = pylandau.langau(x, mpv3, eta3, sigma3, A3)
y3 = pylandau.langau(E, mpv3, eta3, sigma3, A3)
y4 = y2+y3
y6 = y+y5
#yerr = np.random.normal(np.zeros_like(x), np.sqrt(y))
#yerr[y < 1] = 1
#y += yerr


# Fit with constrains
#DH replacing x,y with energies and counts
coeff, pcov = curve_fit(pylandau.langau, E, C,
                        sigma=np.sqrt(C),
                        absolute_sigma=True,
                        p0=(mpv, eta, sigma, A),
                        bounds=(1, 1000))

# Plot
plt.errorbar(E, C, np.sqrt(C), fmt=".")
#plt.plot(E, pylandau.langau(E, *coeff), "-")
plt.plot(x, y6)

chi2 = 0.0
for value1,value2 in zip(C,y4):
    chi2 += math.pow((value1-value2),2)/np.sqrt(value1)

chi2 = chi2/4.
    
print(y2)
print(chi2)

plt.show()
