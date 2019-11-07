# Comments for Maria are prefaced by "MP:"
# Everything else is Dan and a little Cole

# MP: this imports the .so library that contains my C++ chisq function
# In my case, the library is at: ./exposeAPI/ParameterData.so
import sys
#sys.path.insert(0, 'exposeAPI')
#import ParameterData as pd
import csv
import pylandau
import math
# MP: visualization/numerical libraries, some of which are needed for corner
# plots
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as pl
import numpy as np
import corner
import multiprocessing as mp
from scipy import integrate
from matplotlib import rcParams

rcParams["font.size"] = 16
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
#rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"

# MP: import emcee. The MPIPool routine is used for parallelization within
# emcee; if you don't want to parallelize emcee's walkers, no need to import this
import emcee
#from emcee import PTSampler
#from emcee.utils import MPIPool

# if len(sys.argv) < 3:
#     sys.exit("Error: expected 3 or more arguments, but encountered " +
#             str(sys.argv[0:]))

# nucleus = sys.argv[1]
# outputStem = sys.argv[2]


inputStem = "/home/hoff/markovChainMonteCarlo/markovChainMonteCarlo"
outputStem = "outputs/doubletGroundState/thermo"
#outputStem = "outputs/tripletGroundState/thermo"

#read in tab delimited data file
energies=[]
counts=[]
countsErr = np.array([])
#with open(inputStem+'/decayEnergy_100kevBins_filter2.txt','r') as f:
#with open("decayEnergy_100keVBins_filter2.txt",'r') as f:
#with open("decayEnergy_1keVBins_filter.txt",'r') as f:
with open("decayEnergy_1keVBins_filter2.txt",'r') as f: #removed data points below 1 MeV, which is somewhat justified looking at full dataset. This allowed for no extraneous low energy solutions
    #next(f) # skip headings
    reader=csv.reader(f,delimiter='\t')
    for energy,count in reader:
        energies.append(energy)
        counts.append(count)
       
       
E = np.array(energies, dtype=float)
C = np.array(counts, dtype=float)

data = zip(E,C)

# MP: Specifying where to find the prior distribution limits and my "best guess"
# of the parameter values. In the example I showed at the DOM meeting, no "best
# guess" was needed because I just randomly initialized the parameters within
# the prior.
paramFileName = "parameters/current.inp"
limitFileName = "parameters/prior.txt"


theta = []
thetaLimits = []

with open(paramFileName,'r') as f:
    #next(f) # skip headings
    reader=csv.reader(f,delimiter='\n')
    for param in reader:
        theta.append(param)

with open(limitFileName,'r') as f:
    for line in f:
        #line.replace('\t','\n')
        tempLimits = line.split()
        if len(tempLimits)==0:
            continue
        thetaLimits.append(tempLimits)


#print(theta)
#print(thetaLimits)

#for temperature I need this parameter GLOBAL, yeah I know there is probably a better way to do this...
beta = 0.0


#apparently need function to change global variable... weird
def changeBeta(value):
    global beta
    beta = value

def getBeta():
    global beta
    return beta

def printBeta():
    global beta
    print(beta)

# MP: define the natural-log probability of the prior distribution. Given
# "theta" which holds the prior distribution limits, this function returns 0 if
# the parameters are within the limits, and return -infinity if any of the
# parameters are outside their prior-distribution limits. Practically speaking,
# it means that the walkers are forbidden from exiting the prior-distribution
# boundaries
def lnprior(theta):
    for param, limits in zip(theta, thetaLimits):
        if not (float(limits[0]) < param < float(limits[1])):
            return -np.inf

    return 0.0

# MP: define the natural-log probability of the likelihood distribution (that
# is, the likelihood of a parameter set "theta", given the known experimental data).
# Assuming that the experimental data are all independent of each other, this is
# just -chisq/2, where the chisq is calculated in my C++ method
# "generateChiSquare" after the parameters have been read in. The factor of 1/2
# actually doesn't matter for our purposes because it scales all the chisq the
# same, so it doesn't change the walker trajectory (I think...?).


x1 = np.linspace(0.5,5000.5,5000)

def lnlike(theta):
    #parameters.set_params(pd.pyListToVec(theta))
    #params = pd.vector()
    #params[:] = theta
    #parameters.set_params(params)
    #f = pd.fit(parameters)
    #chisq = f.generateChiSquare()
    #print "chisq = " + str(chisq)

    #below is OK for traditional chisquared calculation
    #should find a way to do this dynamically at some point, based on input parameters.
    #y = pylandau.langau(E, float(theta[0]), 64, 45, float(theta[1])) #units in keV, mpv, width, sigma, scale
    #y2 = pylandau.langau(E,float(theta[2]), 64, 45, float(theta[3]))
    #y3 = pylandau.langau(E,float(theta[4]), 64, 45, float(theta[5]))

    #y4 = y+y2#+y3

    chisq = 0.0
    #for value1,value2 in zip(C,y4):
     #   chisq += math.pow((value1-value2),2)/np.sqrt(value1)

    #I have to figure out ssh for this computer, or put all of this on the server at UML!(DUH OF COURSE)
    #At least for the doublet,I have to reduce scale parameter to relative ratio
    #For triplet, I could proceed with same scale parameter and normalize using 3 parameter
    #could reduce parameters by turning scale between 2 low energy peaks into a ratio
    #ratioLowE = theta[1]/theta[5]

    #changing to correlation between scale factors (ratio dependence)
    #could also introduce scale factor
    #below is for "unbinned" log(likelihood) calculation
    
    #for doublet
    #fractionTot = float(theta[2]+theta[3])
    
    #ratio is better constrained
    y = pylandau.langau(x1, float(theta[0]), 64, 45, float(theta[2])) #units in keV, mpv, width, sigma, scale
    y2 = pylandau.langau(x1,float(theta[1]), 64, 45, 1.-float(theta[2])) #if prior between 0 and 1 this should converge to fractional distribution

    # #remember to redo params and priors dan!
    y4 = y+y2


    #for triplet
    # fractionTot = float(theta[3]+theta[4]+theta[5])
    # #fractionTot = float(theta[5])

    #going to make theta[3] the log of the intensity ratio between low energy peak1 to peak2
    #theta[4] is then fraction that goes to high energy peak (I hope the math below works out for the individual fractions)

    # fraction12 = theta[3]
    # # fraction12 = np.exp(theta[3])

    # y = pylandau.langau(x1, float(theta[0]), 64, 45, float(1.-theta[4])/float(1+fraction12)) #units in keV, mpv, width, sigma, scale
    # y2 = pylandau.langau(x1,float(theta[1]), 64, 45, float(theta[4])) #if prior between 0 and 1 this should converge to fractional distribution
    # y3 = pylandau.langau(x1,float(theta[2]), 64, 45, float(fraction12)*float(1.-theta[4])/float(1+fraction12))

    # #remember to redo params and priors dan!
    # y4 = y+y2+y3




    #original
    #y = pylandau.langau(x1, float(theta[0]), 64, 45, float(theta[1])) #units in keV, mpv, width, sigma, scale
    #y2 = pylandau.langau(x1,float(theta[2]), 64, 45, float(theta[3]))
    #y3 = pylandau.langau(x1,float(theta[4]), 64, 45, float(theta[5]))

    #if I get rid of normalize perhaps I will recover absolute heigh information
    #of course it collapses to largest values since this "maximizes" probability
    normalize = np.trapz(y4)
    #normalize = 1.

    for energy in E: #because of binning of x1, energy is mapped -0.5 from the entry
        chisq += np.log(y4[int(energy-0.5)]/normalize) #for likelihood analysis the the probability that a model is correct for n points is by comparing the model distribution, P, at points, p_n, and multiplying them all P(p_1)*P(p_2)*...*P(p_n). Thus when taking the log these can all be added together. This is the natural way to do these kind of fits, though time consuming.

    #I'm appending the temperature parameter beta to the end of the 

    return chisq
    #return -(chisq/2.0) #good for traditional chisquared calculation where Likelihood = exp(-chisq/2.0)

    
# MP: define the natural-log probability of the posterior distribution.
# According to Bayes formula, posterior = likelihood x prior. So lnprior +
# lnlike = lnposterior. Note that I've added the "np.isfinite" command as a
# catch-all: if my chisq function has numerical problems (for example, I start
# out in a weird point in parameter space and my chisq generates a NaN error), I
# can just return -infinity here, and the walkers will stay away from the weird
# point.
def lnprob(theta):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta)


def lnprobBeta(theta):
    lp = lnprior(theta)
    #global beta
    #print(getBeta())
    #printBeta()
    if not np.isfinite(lp):
        return -np.inf
    return lp + getBeta()*lnlike(theta)



#found code for calculating the Bayes Factor for different dimension
def integrate_posterior_3D(log_posterior, xlim, ylim, zlim, data=data):
    func = lambda theta2, theta1, theta0: np.exp(log_posterior([theta0, theta1, theta2], data))
    return integrate.tplquad(func, xlim[0], xlim[1],
                             lambda x: ylim[0], lambda x: ylim[1],
                             lambda x, y: zlim[0], lambda x, y: zlim[1])

def integrate_posterior_6D(log_posterior, xlim, ylim, zlim, data=data):
    func = lambda theta5, theta4, theta3, theta2, theta1, theta0: np.exp(lnlike([theta0, theta1, theta2, theta3, theta4, theta5], data))
    return integrate.tplquad(func, xlim[0], xlim[1],
                             lambda x: ylim[0], lambda x: ylim[1],
                             lambda x, y: zlim[0], lambda x, y: zlim[1])


# MP: set the total number of dimensions to match the number of parameters I'd
# like to sample, and set 4 walkers for each dimension. emcee requires that
# nwalkers be at least twice the number of dimensions, unless you pass the
# argument "live_dangerously=True" to the EnsembleSampler method below.

#ndim = len(variedParams)

ndim = len(theta)
#nwalkers = ndim*4

nwalkers = ndim*8

#the multiprocessing used below is compatable with Mac OS x
#if you want to burn out your system
pool = mp.Pool(processes=6)
#medium
#pool = mp.Pool(processes=2)
#lite
#pool = mp.Pool(processes=1)

# MP: for parallelization of the walkers. You may not want this
#pool = MPIPool()
# if not pool.is_master():
#    pool.wait()
#    sys.exit(0)
#
# MP: the first command below if for the parallelized version. If not
# using parallelized emcee, use the second (commented) command

#sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, pool=pool) 

#sampler for thermo calculation
#sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprobBeta, pool = pool) 
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprobBeta) 


#sampler=PTSampler(ntemps, nwalkers, ndim, lnprob, lnprior, pool=pool)
#sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob)


# MP: The "normal usage" of this MCMC sampler includes a burn-in period, where
# the walkers are allowed to "loosen up" and walk around parameter space for a
# while, getting "close" to the stationary distribution. Then, after the burn-in
# is over, you can start taking "real" samples that tell you what the
# true posterior distribution looks like. The implementation of this approach
# is on the emcee website and I've largely copied it below.

# Per emcee's website, all the calculations used in the "burn-in" are thrown
# away. In my case, I didn't want this behavior - I wanted to see the walkers' position all
# the way through, from beginning to end. So I've commented out the "burn-in"
# below and I just write out the position of each walker after every step. 
# Then at the end, when I can look at the whole progression of steps, I decide what
# the "burn-in" should be, basically when the walker distribution stops moving around.

# perform burn-in
#burnin_steps = 400
#position_output = open(outputStem+"/position.out", "w")
#position_output.close()
#burnin_pos = []

burnin_steps=0

#
#for i, result in enumerate(sampler.sample(p0, iterations=burnin_steps,
#    storechain=False)):
#    position = result[0]
#    position_output = open(outputStem+"/position.out", "a")
#    for k in range (position.shape[0]):
#        walkerPos = [str(var) for var in position[k]]
#        position_output.write("{}\n".format(" ".join(walkerPos).replace("\n","")))
#    position_output.write("\n")
#    position_output.close()
#    burnin_pos = result 
#    if(i+1) % 1 == 0:
#        sys.stdout.write("\rBurn-in progress: {0:5.1%}".format(float(i+1) / burnin_steps))
#        sys.stdout.flush()
#
##pos, prob, state = sampler.run_mcmc(p0, burnin_steps)
#sampler.reset()

#reset is the function I want for the termperature annealing

# MP: This is where I actually take the samples and write out each walker's
# position to "position.out" file after every step.

#really?
#nsteps = 10

#testing
#nsteps = 50

#testing 1
#nsteps = 52

#more testing
#nsteps = 100

#more more testing
#nsteps = 200

#rare
nsteps = 600

#medium rare
#nsteps = 1200

#medium
#nsteps = 2500

#well done
#nsteps = 5000

#go big or go home baby!
#nsteps = 20000

likelihoodPartialMean = []
likelihoodPrior = []
likelihoodPriorPartialMean = []
likelihoodBeta = []
Betas = []

likelihood = []
elikes = []

ntemps = 8;
maxT = 1.
minT = 0.0001
#minT = 0.5 #for testing
deltaX = np.log(maxT/minT)/float(ntemps);
deltaBeta = minT
n = 0
#following procedure outlined in Goggans & Chi (2004) doi: 10.1063/1.1751356


while beta <= 1:

    sampler.reset() #to make sure that there isn't anything else in there before I start
    Betas.append(beta)

    print("\nbeta: {0:0.6f}".format(getBeta()))

# MP: initialize parameters. The first command uniformly samples each parameter
# within the limits defined in the prior. The second (commented) command samples
# each parameter according to a Gaussian with sigma of 5% of the "best fit"
# parameter value, where I take the "best fit" from the DOM's powell method,
# like Mack and I have been using.

#p0 = [[(float(limits[1])-float(limits[0]))*np.random.random_sample()+float(limits[0]) for limits in thetaLimits] for i in range(nwalkers)]

    #for thermo
    p0 = [[(float(limits[1])-float(limits[0]))*np.random.random_sample()+float(limits[0]) for limits in thetaLimits] for i in range(nwalkers)]


    #need extra set of walkers for each "temperature" in PTSampler
    #p0 = [[[np.random.normal((float(limits[1])+float(limits[0]))/2.,(float(limits[1])-float(limits[0]))/7.) for limits in thetaLimits] for i in range(nwalkers)]for k in range(ntemps)]

    #print(p0)
    #p0 = [[variedParams[j] + np.random.normal(scale=0.05*variedParams[j]) for j in range(ndim)] for i in range(nwalkers)]


    minwj = 1E10
    maxwj = 0.
    maxE = 0.
    minE = 1E10
    position = []
    avgLpos = 0.
    avgLPrior = 0.
    avgElike = 0

    for i, result in enumerate(sampler.sample(p0,iterations=nsteps)):
        #position = result[0]
        position = sampler.get_chain()[-1]
        logprob = sampler.get_log_prob()[-1]
        #print(result)
        #print(i)
        #print(position[i])
        #print(position)
        #print(logprob)
        try:
            position_output = open(outputStem+"/position.out", "a")
        except IOError:
            continue
            
        sum = 0.
        

 #       if(i > 100):#to not get burn in

#            for k in range (position.shape[0]):
                #attempting to track the likelihood calculated for each walker, the running average should converge to the relative Bayes Factor (can directly use ratio if prior probabilities are equal)
                
                #should really generate model comparisons after the fact by looking at the walker positions

                #using "temperature" method, before was using running average of likelihood, I'm not sure what's better or more correct
                #Lpos = np.exp(lnlike(position[k]))
                #avgLpos += Lpos/float(nwalkers)
            
                # Elike = -lnlike(position[k])
                # avgElike +=  Elike/float(nwalkers)

                # wj = np.exp(-deltaBeta*Elike)
                # #print(wj)
                # #walkerPos = [str(var) for var in position[k]]
                # #position_output.write("{}\n".format(" ".join(walkerPos).replace("\n",""))
                # if(Elike > maxE):
                #     maxE = Elike
                # if(Elike < minE):
                #     minE = Elike
                # if(wj > maxwj):
                #     maxwj = wj
                # if(wj < minwj):
                #     minwj = wj

        
                # #avgLpos = avgLpos/nwalkers
                # #likelihood.append(avgLpos)
                # elikes.append(avgElike)

            #likelihoodPartialMean.append(np.mean(likelihood))

            #I'm going to sample the likelihood from the prior distribution to extract true Bayes factor
            #p1 = [(float(limits[1])-float(limits[0]))*np.random.random_sample()+float(limits[0]) for limits in thetaLimits]
            #gaussian priors
            #p1 = [np.random.normal((float(limits[1])+float(limits[0]))/2.,(float(limits[1])-float(limits[0]))/7.) for limits in thetaLimits]
            
            #print(p1)
            #likelihoodPrior.append(np.exp(lnlike(p1)))
            #likelihoodPriorPartialMean.append(np.mean(likelihoodPrior))


        
            #position_output.write("\n")
            #position_output.close()
        if(i+1) % 1 == 0:
            sys.stdout.write("\rSampling progress: {0:5.1%}".format(float(i+1) / nsteps))
            sys.stdout.flush()


    for k in range(nwalkers):
        #Elike = -sampler.get_log_prob(flat=True)[nwalkers*(nsteps-1)+k]
        #Elike = -sampler.get_log_prob()[-1][k] #I can't take log_prob from chain because this has beta multiplied by it, I need actuall prob distribution
        Elike = -lnlike(position[k])
        #print(Elike)
        #ElikeTest = -lnlike(sampler.get_chain(flat=True)[:,k])
        #print(ElikeTest)
        #Elike = -lnlike(position[k])
        avgElike +=  Elike/float(nwalkers)
        
        wj = np.exp(-deltaBeta*Elike)
        #print(wj)
        #walkerPos = [str(var) for var in position[k]]
        #position_output.write("{}\n".format(" ".join(walkerPos).replace("\n",""))
        if(Elike > maxE):
            maxE = Elike
        if(Elike < minE):
            minE = Elike
        if(wj > maxwj):
            maxwj = wj
        if(wj < minwj):
            minwj = wj
                        
        
        #avgLpos = avgLpos/nwalkers
        #likelihood.append(avgLpos)

    elikes.append(avgElike)

    # print("\n")
    # print(maxE)
    # print(minE)
    # print(maxwj)
    # print(minwj)

    print("\n")
    #likelihoodBeta.append(np.mean(likelihood))
    if(beta == 1): #just to make sure that it kicks out
        break

    
    if(beta > 0):
        deltaBeta = np.exp(np.log(minT)+float(n)*deltaX) - np.exp(np.log(minT)+float(n-1)*deltaX);

    #print(deltaBeta)
    #beta = np.exp(np.log(minT)+float(n)*deltaX);

    changeBeta(beta + deltaBeta)
    #print(getBeta())
    #deltaBetaTest = np.log(maxwj/minwj)/(maxE - minE)
    #print(deltaBetaTest)
    #print(np.log(maxwj/minwj)/(maxE-minE))
    n = n + 1
    if(beta > 1): #just to make sure that I get the beta=1
        beta = 1


# MP: This should be between 0.25 and 0.5, according to the emcee website.
print("\n Mean acceptance fraction: {0:0.3f}\n"
        .format(np.mean(sampler.acceptance_fraction)))

# print("\n Likelihood from prior: {0:0.5f}\n"
#       .format(np.mean(likelihoodPrior)))

#print("\n Likelihood of posterior:")# {0:0.5d}\n"
#print(np.mean(likelihood))

#print(likelihoodBeta)
#print(elikes)
#print("\n Bayes Factor from PTSampler:")# {0:.3f}"
#       .format(sampler.thermodynamic_integration_log_evidence()[0]))

#print(sampler.thermodynamic_integration_log_evidence())

#print(np.sum(likelihood)/len(likelihood))
#pool.close()

# MP: below is the code for creating a corner plot. Basically the same as the
# example on the emcee website
samples = sampler.chain[:, burnin_steps:, :].reshape((-1,ndim))
fig = corner.corner(samples)#, labels=variedParamNames, truths=variedParams)
#fig = corner.corner(samples, quantiles=(0.16,0.84) )#, labels=variedParamNames, truths=variedParams)
fig.savefig(outputStem+"/triangle.png",dpi=600)

# MP: this gives 1-D histograms for each dimension
for i in range(ndim):
    #this will draw on bottom corner plot if line below is commented out
    #pl.figure()
    
    #for regular emcee use
    probDensities, bins, bars = pl.hist(sampler.flatchain[:,i], 100, color="k", histtype="step")

    #for PTSampler we need to truncate off the temperatures and ndim
    #probDensities, bins, bars = pl.hist(sampler.flatchain[-1][:,i], 100, color="k", histtype="step")
    #pl.title("Dimension {0:d}".format(i))
    #print(bins)
    #print(probDensities)
    #np.savetxt('test.txt', data, delimiter='\t')
    
    with open(outputStem+"/probDensity"+str(i)+".txt", 'w') as f2:
        writer = csv.writer(f2,delimiter='\t')
        writer.writerows(zip(bins, probDensities))


pl.figure()
pl.plot(Betas, elikes,'bo')
pl.savefig(outputStem+"/annealingLikelihood")

# pl.figure()
# pl.plot(likelihoodPriorPartialMean,'bo')
# pl.savefig(outputStem+"/priorLikelihood")

pl.show()
