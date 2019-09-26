import csv
import numpy as np
import matplotlib.pyplot as pl
import matplotlib.animation as ani
import pylandau

outputStem = "outputs/tripletGroundState"

params0 = []
params1 = []
params2 = []
params3 = []
params4 = []
params5 = []

params = []



with open(outputStem+'/position.out','r') as f:
    #next(f) # skip headings
    reader=csv.reader(f,delimiter=' ')
    for parameters in reader:
        #params0.append(param0)
        #params1.append(param1)
        #params2.append(param2)
        #params3.append(param3)
        params.append(parameters)
        
#print(params)
#print(params[0][0])

#remove empty strings from list
params = list(filter(None,params))

#if(params[23]): #this will return false if empty
#print(params[24])

nwalkers = 24
nsteps = 1000

time = []
walker1 = []
walker2 = []
walker3 = []

#param0walk = np.arange(nsteps)
param0walk = []

for j in range(nsteps):
    #print(j)
    time.append(j)
    walker1.append(params[j*nwalkers+2])
    #update_line(hl, j, float(walker1[j][1]))
    #for i in range(nwalkers):
    #pl.plot(float(j),float(walker1[j][1]),'b.')
    param0walk.append(float(walker1[j][0]))

param0walk = [float(i) for i in param0walk]
    
param0walkTest = []
param0walkTest.append(time)
param0walkTest.append(param0walk)
fig1 = pl.figure()

#for creating animation plot


#def update_line(num, new_ydata, hl):
#    hl.set_xdata(np.append(hl.get_xdata(), num))
#    hl.set_ydata(np.append(hl.get_ydata(), new_ydata))
#    return hl
    #pl.draw()

def update_line(num, data, line):
    line.set_data(data[...,:num])
    #line.set_data(data[num])
    return line,

#print(param0walkTest)

hl, = pl.plot([], [], 'r-')
    
line_ani = ani.FuncAnimation(fig1, update_line, nsteps, fargs=(param0walkTest,hl),
                                  interval=50, blit=True)

#def update_line(num, data, line):
#    line.set_data(data[..., :num])
#    return line,

fig1 = pl.figure()

data = np.random.rand(2, 25)
l, = pl.plot([], [], 'r-')
pl.xlim(0, 1)
pl.ylim(0, 1)
pl.xlabel('x')
pl.title('test')

print(data)

#line_ani = ani.FuncAnimation(fig1, update_line, 25, fargs=(data, l),
 #                                  interval=50, blit=True)
                                   

#pl.plot(time,walker1[][0])
#print(walker1)
pl.show()
