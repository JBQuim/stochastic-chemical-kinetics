import matplotlib.pyplot as plt
import math
import numpy as np
import scipy.special
import random
import sys
import configparser

# ------------
# reading from config.txt
# -------------

# runs - scalar - number of runs to be simulated
# segmentCount - scalar - number of segments x axis is divided into
# percent - scalar - percent of the longest run to graph, ends of 100% can be noisy
# as not many runs reach that far

# varNames - 1xs, string - names of each variable
# c - 1xr - stochastic rate constants
# stoichR - rxs - coefficient of reactants involved in elementary step
# stoichP - rxs - change in number of species after reaction
# n0 - 1xs - starting amount of each reactant
# Tf - scalar - final time
# maxEvents - scalar - number of events at which the simulation will be terminated
# should be larger than the time needed to reach Tf

config = configparser.ConfigParser()
config.read('config.txt')

varNames = np.array(config['parameters']['varNames'].split(','))
c = np.array(list(map(float, config['parameters']['c'].split(','))))

s = len(varNames)
r = len(c)

stoichR = np.resize(np.array(list(map(int, config['parameters']['stoichR'].split(',')))), (r, s))
stoichP = np.resize(np.array(list(map(int, config['parameters']['stoichP'].split(',')))), (r, s))
n0 = np.array(list(map(float, config['parameters']['n0'].split(','))))

Tf = float(config['parameters']['Tf'])
maxEvents = int(config['parameters']['maxEvents'])
runs = int(config['parameters']['runs'])
seed = int(config['parameters']['seed'])

# booleans determine whether to:
# draw the average run
# draw the standard deviation shadow
# draw all the points  across all runs
# draw histograms
averages = config['parameters'].getboolean('averages')
deviations = config['parameters'].getboolean('deviations')
scatter = config['parameters'].getboolean('scatter')
histogramAmounts = config['parameters'].getboolean('histogramAmounts')

segmentCount = int(config['parameters']['segmentCount'])
percent = float(config['parameters']['percent'])

binCount = np.array(list(map(int, config['parameters']['bincounts'].split(','))))

# number of individual runs to draw
lines = int(config['parameters']["lines"])

# save image or not?
save = config['parameters'].getboolean('save')


# ----------------
# declaring functions
# -----------------

# calculates propensities from stochastic rate constants and n
def calcPropensities(number):
    # h[j] is the number of combinations of the Rj reactant molecules
    # h[j] is the product of all the indexes of y
    h = np.zeros(r)
    for j in range(0, r):
        y = np.zeros(s)
        for u in range(0, s):
            if stoichR[j][u] > number[u]:
                y[u] = 0
            else:
                y[u] = scipy.special.comb(number[u], stoichR[j][u])
        h[j] = np.prod(y)

    return c * h


# returns the value of j for which a_j-1/a0 < u2 < a_j/a0
def chooseReaction(a, a0):
    u2 = random.random()
    for j in range(0, r):
        u2 -= a[j] / a0
        if u2 <= 0:
            return j

def runSim(finalTime, finalEvents, initialN):
    # initialize simulation variables
    t = 0.0
    n = initialN
    eventCount = 0

    # history of the system
    history = np.full((finalEvents, s + 1), np.nan)
    history[0] = np.insert(n, 0, t)

    # simulate until t reaches final time or events reaches maxEvents
    while not (t > finalTime or eventCount + 1 >= finalEvents):
        a = calcPropensities(n)
        # a0 is the propensity for any reaction to happen
        a0 = np.sum(a)
        # if the propensities sum to 0 then no reaction is possible
        if a0 <= 0:
            break

        # tau is the time to next reaction and is exponentially distributed
        u1 = random.random()
        tau = -math.log(u1) / a0
        t += tau
        if t > finalTime:
            break

        # change system
        reaction = chooseReaction(a, a0)
        n = n + stoichP[reaction]  # warning! n+=stoichP[reaction] changes n0!!!!
        eventCount += 1

        history[eventCount] = np.insert(n, 0, t)

    return history


# takes in the data from various runs, splits into segments and finds the mean and stdev of each segment
def splitData(values, varNumber, segmentNumber, percentile):
    # sort values according to the time they occur
    order = values[0].argsort()
    values = (values.T[order]).T
    notNan = ~np.isnan(values[0])

    # now produce means and stds
    means = np.full((varNumber + 1, segmentNumber), np.nan)
    stdDevs = np.full((varNumber + 1, segmentNumber), np.nan)

    segmentWidth = np.percentile(values[0][notNan], percentile) / segmentNumber
    for h in range(0, segmentNumber):
        # sliceArgs is a boolean array that is true at the elements that are part of the slice
        # values needs to be cleaned up of nans before use, -1 will return false as its not in any slice
        values[0][~notNan] = -1
        sliceArgs = (values[0] < segmentWidth * (h + 1)) & (values[0] >= segmentWidth * h)

        for j in range(0, varNumber + 1):
            means[j][h] = np.mean(values[j][sliceArgs])
            stdDevs[j][h] = np.std(values[j][sliceArgs])

    return means, stdDevs

# generate dataset of the final state of each variable, for histogram generation
def getEndings(values, runCount, varNumber):
    endings = np.full((runCount, varNumber+1), 0)

    for h in range(0,runCount):
        #data without nan
        refinedData = values[h].T[~np.isnan(values[h][0])]
        endings[h] = refinedData[-1]

    return endings.T

# -------------
# initialize system
# ------------

# initial conditions
# the seed is saved and displayed in case the exact run wants to be replicated
if seed == -1:
    seed = random.randrange(sys.maxsize)
random.seed(seed)

# -----------
# actually run simulation
# ----------

# save every run to an array
rawData = np.zeros((runs, s + 1, maxEvents))
for d in range(0, runs):
    rawData[d] = runSim(Tf, maxEvents, n0).T

# generate histogram
if histogramAmounts:
    histogramData = getEndings(rawData, runs, s)

# join all the runs together
data = np.concatenate(rawData, axis=1)

# process data into segments, take averages and standard deviations
processedData = splitData(data, s, segmentCount, percent)

# plot results
if histogramAmounts:
    fig, (ax1,ax2) = plt.subplots(1,2,sharey=True,gridspec_kw={'width_ratios': [3, 1]})
else:
    fig, ax1 = plt.subplots(1,1)

fig.set_size_inches(12, 6)

for i in range(0, s):
    if averages:
        ax1.plot(processedData[0][0], processedData[0][i + 1])
    if deviations:
        ax1.fill_between(processedData[0][0], processedData[0][i + 1] + processedData[1][i + 1],
                         processedData[0][i + 1] - processedData[1][i + 1], alpha=0.5)
    if scatter:
        ax1.scatter(data[0], data[i + 1], s=5, alpha=0.03)

    cmap = plt.get_cmap("tab10")
    for m in range(0, lines):
        for i in range(0, s):
            ax1.plot(rawData[m][0], rawData[m][i + 1], color=cmap(i),alpha=0.5)
    ax1.set(xlabel = "Time",ylabel = 'Number of molecules')
    ax1.legend(varNames)
    ax1.set_title('Runs: ' + str(runs) + '     Seed used: ' + str(seed))


if histogramAmounts:
    for i in range(0,s):
        ax2.hist(histogramData[i+1], binCount[i],orientation="horizontal",alpha=0.8)
    ax2.set_title('Final conditions over all the runs')


if save:
    plt.savefig('Images\\' + str(seed) + '.png')

plt.show()
