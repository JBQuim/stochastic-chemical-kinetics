import matplotlib.pyplot as plt
import math
import numpy as np
from scipy.special import comb
import random
import sys
import configparser
from matplotlib.lines import Line2D
import time

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
if lines > runs: lines = runs
if lines == -1: lines = runs

# save image or not?
save = config['parameters'].getboolean('save')


# ----------------
# declaring functions
# -----------------

def checkError():
    # not an exhaustive list of all possible errors in config file, but a list of the common ones
    error = "Starting..."
    if np.shape(stoichP) != np.shape(stoichR):
        error = "stoichP and stoichR must be the same size"
    elif s != len(binCount):
        error = "bincounts must be same size as varNames"
    elif s != len(n0):
        error = "n0 must be same size as varNames"
    elif not 100 >= percent >= 0:
        error = "percent must be between 0 and 100"
    elif lines == 0 and ~scatter and ~deviations and ~averages:
        error = "Must choose whether to graph lines, deviations, averages or scatter graph"
    print(error)
    if error != "Starting...": quit()


# calculates propensities from stochastic rate constants and n
def calcPropensities(number):
    # h is the number of combinations of the Rj reactant molecules
    # h is the product of all the indexes of y for every reaction
    number = number[:, None]
    y = comb(number, stoichR.T).T
    h = np.prod(y, axis=1)
    return h * c


# chooses an integer between 0 and r, weighted according to propensities
def chooseReaction(a, a0):
    return np.random.choice(range(0, r), p=a / a0)


def runSim(finalTime, finalEvents, initialN, randomNumbers):
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
        u1 = randomNumbers[eventCount]
        tau = -math.log(u1) / a0
        t += tau
        if t > finalTime:
            break

        # change system
        reaction = chooseReaction(a, a0)
        n = n + stoichP[reaction]  # warning! n+=stoichP[reaction] changes n0!!!!
        eventCount += 1

        history[eventCount] = np.insert(n, 0, t)

    # eventCount is used to keep track of how many random numbers were used
    return history, eventCount


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
    endings = np.full((runCount, varNumber + 1), 0)

    for h in range(0, runCount):
        # data without nan
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
checkError()
rawData = np.zeros((runs, s + 1, maxEvents))
# generate random numbers
randomData = np.random.rand(maxEvents * runs)
numbersUsed = 0
for d in range(0, runs):
    data, used = runSim(Tf, maxEvents, n0, randomData[numbersUsed:numbersUsed + maxEvents])
    rawData[d] = data.T
    numbersUsed += used
    print("Run " + str(d + 1) + "/" + str(runs) + " completed")

# generate histogram
if histogramAmounts:
    histogramData = getEndings(rawData, runs, s)

# join all the runs together
data = np.concatenate(rawData, axis=1)

# process data into segments, take averages and standard deviations
processedData = splitData(data, s, segmentCount, percent)

# plot results
if histogramAmounts:
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, gridspec_kw={'width_ratios': [3, 1]})
else:
    fig, ax1 = plt.subplots(1, 1)

fig.set_size_inches(12, 6)

customLines = []
for i in range(0, s):
    cmap = plt.get_cmap("tab10")
    # creates colours for legend
    customLines.append(Line2D([0], [0], color=cmap(i)))

    if averages:
        ax1.plot(processedData[0][0], processedData[0][i + 1], color=cmap(i))
    if deviations:
        ax1.fill_between(processedData[0][0], processedData[0][i + 1] + processedData[1][i + 1],
                         processedData[0][i + 1] - processedData[1][i + 1], alpha=0.5, color=cmap(i))
    if scatter:
        ax1.scatter(data[0], data[i + 1], s=5, alpha=0.03, color=cmap(i))

    for m in range(0, lines):
        ax1.plot(rawData[m][0], rawData[m][i + 1], color=cmap(i), alpha=0.7)

ax1.legend(customLines, varNames)
ax1.set(xlabel="Time", ylabel='Number of molecules')
ax1.set_title('Runs: ' + str(runs) + '     Seed used: ' + str(seed))

if histogramAmounts:
    for i in range(0, s):
        ax2.hist(histogramData[i + 1], binCount[i], orientation="horizontal", alpha=0.8)
    ax2.set_title('Final conditions over all the runs')

if save:
    plt.savefig('Images\\' + str(seed) + '.png')

plt.show()
