#!/usr/bin/env python2

#
# Fitting experimental SWAXS curves to calcualted curves
#
# Written by Jochen Hub in 2013 and 2014, all rights reserved.
#
# If you use this script for your research, please cite:
#
# Christopher J Knight, and Jochen S Hub
# WAXSiS - a web server for the calcualtion of SAXS/WAXS curves based on explicit-solvent molecular dynamics
# Nucleic Acid Research (2015)
#
# This software comes with ABSOLUTELY NO WARRANTY, to the extent permitted by applicable law.
#

import sys, math, numpy, optparse, random
from scipy.optimize import fmin_powell

#
# Some simple IO to xmgrace format
#
def read_xvg(fn, qmin = -1, qmax = -1):
    fp   = open(fn, 'r')
    x    = []
    y    = []
    yerr = []
    for line in fp.readlines():
        if (not (line[0] == '#' or line[0] == '@' or line[0] == '&' or line == '\n')):
            l = line.split()
            xp = float(l[0])
            bInRange = ((qmin < 0 or xp >= qmin) and (qmax < 0 or xp <= qmax))
            if (bInRange):
                x.append(xp)
                y.append(float(l[1]))
                if (len(l) > 2):
                    yerr.append(float(l[2]))
    fp.close()
    return [x, y, yerr]
    
def write_to_xvg(fp, x, y, yerr):
    if (len(yerr) > 0):
        print >> fp, '@type xydy'
        for i in range(len(x)):
            print >> fp, '%g %g %g' % (x[i], y[i], yerr[i])
    else:
        print >> fp, '@type xy'
        for i in range(len(x)):
            print >> fp, '%g %g' % (x[i], y[i])


############################################################################################
# List of functions used to compute Chi2
############################################################################################

def intensityDiff_fac(pos, *args):
    # Workaround since I don't know how to pass the list of lists with the star at *args
    if (len(args) == 1):
        args = args[0]
    y        = args[0]
    ysig     = args[1]   # is ignored in this function
    ysam     = args[2]
    ysamsig  = args[3]
    ybuf     = args[4]
    ybufsig  = args[5]
    bUseWeights = args[6]
    alpha0   = args[7]  # 1 - volume fraction of the solute
    sigAlpha = args[8]  # Uncertainty of the parameter a, see eq. 4 Grishaev and Bax, JACS 2010, 132:15484-15486

    if (len(y) != len(ysam) or len(y) != len(ybuf) or len(y) != len(ysamsig) or len(y) != len(ybufsig) ):
        printf >> sys.stderr, 'Different lengths of y, ysam, or ybuf'
        sys.exit(1)

                 # I(exp) = f * ( ysam - a*ybuf + c)
    f = pos[0]   # The overall scaling of the experimental intensity
    a = pos[1]   # Correction factor for buffer intensity
    c = pos[2]   # constant added to experimental intensity
                 # print 'Inside logIntensityDiff: f / c = %g / %g' % (f,c)
    #print 'Now f / a / c = %g %g %g' % (f,a,c)
    #print 'Now alpha0 sig bUseWeights= %g %g %d' % (alpha0, sigAlpha, bUseWeights)
    chi2 = 0
    for i in range(len(y)):
        ynet = f*(ysam[i] - a*ybuf[i] + c)
        if (ynet < -1e5):
            # don't allow large negative I(q)
            return 1e20
        tmp = (y[i] - ynet)**2
        if (bUseWeights):
            ynetSigma2 = f**2 * (ysamsig[i]**2 + (a*ybufsig[i])**2)
        else:
            ynetSigma2 = 1.
        chi2 += tmp/ynetSigma2

    chi2 /= len(y)
    # Regularization factor (a - a0)^2/sigAlpha^2
    if (sigAlpha > 0):
        chi2 += ((a - alpha0)/sigAlpha)**2
        #print ((a - alpha0)/sigAlpha)**2
        
    return chi2
 
############################################################################################

def IntensityDiff_fc(pos, *args):
    y1    = args[0]
    y1sig = args[1]   # is ignored in this function
    y2    = args[2]
    y2sig = args[3]
    bUseWeights = args[4]
    if (len(y1) != len(y2)):
        printf >> sys.stderr, 'Different lengths of y1 and y2'
        sys.exit(1)

    if (bUseWeights and (len(y2sig) != len(y2))):
        printf >> sys.stderr, 'Error in IntensityDiff_fc(): requested weighted sum, but sigma arrays are empty'
        sys.exit(1)


                 # I(exp) = f * y2 + c
    f = pos[0]   # The overall scaling of the experimental intensity
    c = pos[1]   # constant added to experimental intensity
                 # print 'Inside logIntensityDiff: f / c = %g / %g' % (f,c)
    sumChi2 = 0.
    for i in range(len(y1)):
        thisy2 = f*y2[i] + c
        tmp    = (y1[i]-thisy2)**2
        if (bUseWeights):
            sig2 = (f*y2sig[i])**2
        else:
            sig2 = 1.
        sumChi2 += tmp/sig2
        # print '%d  / %12g %12g - chi2 %12g' % (i, tmp, sig2, sumChi2)
    # print 'In function, returning %g' % (sumChi2/len(y1))
    return sumChi2/len(y1)


def logIntensityDiff_fc(pos, *args):
    y1    = args[0]
    y1sig = args[1]   # is ignored in this function
    y2    = args[2]
    y2sig = args[3]
    bUseWeights = args[4]
    if (len(y1) != len(y2)):
        printf >> sys.stderr, 'Different lengths of y1 and y2'
        sys.exit(1)

                 # I(exp) = f * y2 + c
    f = pos[0]   # The overall scaling of the experimental intensity
    c = pos[1]   # constant added to experimental intensity
                 # print 'Inside logIntensityDiff: f / c = %g / %g' % (f,c)
    sumChi2 = 0.
    sig2    = 1.
    for i in range(len(y1)):
        thisy2    = f*y2[i] + c
        # doing the difference between the logarithms of the intensities
        if (thisy2 < 0):
            # don't allow negative I(q)
            return 1e20
        tmp  = (math.log(y1[i])-math.log(thisy2))**2
        if (len(y2sig) > 0 and bUseWeights):
            # Note sigma[log(I)] = sig(I)/I
            sig2 = (f*y2sig[i]/thisy2)**2

        sumChi2 += tmp/sig2
        # print '%d  / %12g %12g - chi2 %12g' % (i, tmp, sig2, sumChi2)
    # print 'In function, returning %g' % (sumChi2/len(y1))
    return sumChi2/len(y1)

def intensityDiff_SAXS(pos, *args):
    y1    = args[0]
    y2    = args[2]
    y2sig = args[3]
    bUseWeights = args[4]
    if (len(y1) != len(y2)):
        printf >> sys.stderr, 'Different lengths of y1 and y2'
        sys.exit(1)
    f = pos[0]   # The overall scaling of the experimental intensity
    print 'Inside intensityDiff_SAXS: f = %g' % (f)
    sumChi2 = 0.
    sig2    = 1.
    for i in range(len(y1)):
        thisy2   = f*y2[i]
        tmp      = (y1[i]-thisy2)**2
        if (len(y2sig) > 0 and bUseWeights):
            sig2 = (f*y2sig[i])**2

        sumChi2 += tmp/sig2
        # print '%d  / %12g %12g - chi2 %12g' % (i, tmp, sig2, sumChi2)
    print 'In function, returning %g' % (sumChi2/len(y1))
    return sumChi2/len(y1)


def logIntensityDiff_SAXS(pos, *args):
    y1    = args[0]
    y2    = args[2]
    y2sig = args[3]
    bUseWeights = args[4]
    if (len(y1) != len(y2)):
        printf >> sys.stderr, 'Different lengths of y1 and y2'
        sys.exit(1)
    f = pos[0]   # The overall scaling of the experimental intensity
    print 'Inside logIntensityDiff: f = %g' % (f)
    sumChi2 = 0.
    for i in range(len(y1)):
        thisy2    = f*y2[i]
        # doing the difference between the logarithms of the intensities
        if (thisy2 < 0):
            # don't allow negative I(q)
            return 1e20
        tmp  = (math.log(y1[i])-math.log(thisy2))**2
        sumChi2 += tmp
        # print '%d  / %12g %12g - chi2 %12g' % (i, tmp, sig2, sumChi2)
    print 'In function, returning %g' % (sumChi2/len(y1))
    return sumChi2/len(y1)


# Difference between sample and buffer
def sampleBufferDiff(xsam, ysam, ysamerr, xbuf, ybuf, ybufSig, alpha=1):
    # First make sure that the sample x-range is within the buffer x-range
    while (xsam[0] < (xbuf[0]-0.0001)):
        xsam = xsam[1:]
        ysam = ysam[1:]
        if (len(ysamerr) > 0):
            ysamerr = ysamerr[1:]
    while (xsam[-1] > (xbuf[-1]+0.0001)):
        xsam = xsam[:-1]
        ysam = ysam[:-1]
        if (len(ysamerr) > 0):
            ysamerr = ysamerr[:-1]
    
    x0 = xsam[0]
    x1 = xsam[-1]
    nx = len(xsam)

    ybufInterp    = numpy.interp(xsam, xbuf, ybuf)
    if (len(ybufSig) > 0 and len(ysamerr) > 0):
        ybufErrInterp = numpy.interp(xsam, xbuf, ybufSig)
    xNet    = []
    yNet    = []
    yNetErr = []
    for i in range(nx):
        xNet.append(xsam[i])
        yNet.append(ysam[i] - alpha*ybufInterp[i])
        if (len(ybufSig) > 0 and len(ysamerr) > 0):
            sig = math.sqrt(ysamerr[i]**2 + alpha**2 * ybufErrInterp[i]**2)
            yNetErr.append(sig)
    return [xNet, yNet, yNetErr]


def chi_free(x, ycal, ycalErr, yexp, yexpErr, D, nfree, bLogScale, bUseWeights):
    ns = int( round(1.0*(x[-1]-x[0])*options.D/math.pi) )

    nX = 1.0 * len(x) / ns
    print 'qmax = %8.3g 1/nm  -- qmin = %8.3g 1/nm -- D = %8.3g nm -- ns = %d' % (x[-1], x[0], D, ns)
    print 'Having %8.3g data points per bin, %d data points total' % (nX, len(x))

    if (D < 0):
        print >> sys.stderr, 'Error, cannot get chi(free) without knowing the diameter'
        sys.exit(1)

    f    = []
    c    = []
    chi2 = []
    
    for k in range(nfree):
        xT       = []
        yT       = []
        yTerr    = []
        yTexp    = []
        yTexpErr = []
        for i in range(ns):
            j = int( nX * (i + random.random()) )
            if (j >= len(x)):
                print >> sys.stderr, 'WARNING, got j = %d of %d, reducing' % (j, len(x))
                j = len(x)-1
            xT.      append(x      [j])
            yT.      append(ycal   [j])
            yTerr.   append(ycalErr[j])
            yTexp.   append(yexp   [j])
            yTexpErr.append(yexpErr[j])
            
            pos  = [y[0]/ynet[0], 1.]
        if (bLogScale == False):
            fminOut = fmin_powell(IntensityDiff_fc,    pos, args=(yT, yTerr, yTexp, yTexpErr, bUseWeights), full_output=True, disp=0)
        else:
            fminOut = fmin_powell(logIntensityDiff_fc, pos, args=(yT, yTerr, yTexp, yTexpErr, bUseWeights), full_output=True, disp=0)
        f.append(fminOut[0][0])
        c.append(fminOut[0][1])
        chi2.append(fminOut[1])

        # print '%12.4g  %12.4g  %12.4g' % (fminOut[0][0], fminOut[0][1], chi2[-1])

    # print chi2
    chi2Median = numpy.median(chi2)
    chiMedian  = math.sqrt(chi2Median)
    chi2Av     = sum(chi2)/len(chi2)
    tmp = 0
    for i in range(len(chi2)):
        tmp += (chi2[i]-chi2Av)**2
    chi2Stddev = math.sqrt(tmp/len(chi2))

    print 'Average of chi2      = %8.4g' % (chi2Av)
    print 'Stddev  of chi2      = %8.4g' % chi2Stddev
    print 'Median of chi2(free) = %8.4g' % (chi2Median)
    print 'Median of chi (free) = %8.4g' % (chiMedian)



#####################################
# MAIN PROGRAM ######################

FIT_F   = 0
FIT_FC  = 1
FIT_FAC = 2

parser = optparse.OptionParser(conflict_handler="error")

parser.add_option("--nofit",
                  action="store_const", const=FIT_F , dest="fitType",
                  help="Only fit the overall scaling parameter f for the net intensity: I[exp,fit] = f * [I[sam] - I[buf]]", default=1)

parser.add_option("--fc",
                  action="store_const", const=FIT_FC, dest="fitType",
                  help="Fit I[exp,fit] = f * [I[sam] - I[buf] + c] to experient, using NO experimental errors.",
                  default=1)

parser.add_option("--fac",
                  action="store_const", const=FIT_FAC, dest="fitType",
                  help="Fit I[exp,fit] = f *[I[sam] - a*I[buf] + c] to experient, using experimental errors.",
                  default=1)

parser.add_option("-c", "--calc",
                  action="store", type="string", dest="fnSim", help="Calculated (simulation) intensity")

parser.add_option("-n", "--net",
                  action="store", type="string", dest="fnExpNetIntensity", help="Experiental net intensity, I[sam]-I[buf]")

parser.add_option("-s", "--sample",
                  action="store", type="string", dest="fnSample", help="Experiental sample intensity, only for --fac")

parser.add_option("-b", "--buffer",
                  action="store", type="string", dest="fnBuffer", help="Experiental buffer intensity, only for --fac")

parser.add_option("--qmin", action="store", type="float", dest="qmin", help="Minimum q to fit", default='0')

parser.add_option("-q", "--qmax", action="store", type="float", dest="qmax", help="Maximum q to fit", default='-1')

parser.add_option("-w", "--use-weights", action="store_true", dest="bUseWeights",
                  help="Use experimental errors as inverse fitting weights", default=False)

parser.add_option("--log", action="store_true", dest="bLogScale",
                  help="Fit on a logarithmic scale", default=False)

parser.add_option("--volfrac", action="store", type="float", dest="volSolute",
                      help="Fractional volume of solute, only relevant for --fac: Buffer subtracted as I[exp] = I[sam]-(1-v)*I[buf]", default='0')

parser.add_option("--volfracerr", action="store", type="float", dest="volSoluteSig",
                  help="Uncertainty in fractional volume of solute", default='0.01')

parser.add_option("--exp-ang", action="store_true", dest="bExpAngstroem",
                  help="Inverse Angstrom units in experimental files", default=False)

parser.add_option("-o", "--fit-overlay",
                  action="store", type="string", dest="fnOverlay", help="xvg file overlaying calculated and fitted",
                  default='fittedOverlay.xvg')

parser.add_option("--fit-overlay-expQ",
                  action="store", type="string", dest="fnOverlayExpQ", help="xvg file overlaying calculated and fitted with experimental q-points")

parser.add_option("--Icalc-at-expQ",
                  action="store", type="string", dest="fnIcalcExpQ", help="calculated I(q) at exprimental q-points")

parser.add_option("-f", "--fitted",
                  action="store", type="string", dest="fnFitted", help="xvg file of fitted",
                  default='fittedExp.xvg')

parser.add_option("-D",
                  action="store", type="float", dest="D", help="Diameter of solute, used only for -chi-free analysis", default=-1)

parser.add_option("--nfree",
                  action="store", type="int", dest="nfree", help="Number of rounds for chi(free) calculation (Note: chi-free analysis is still controversial.)", default=0)


parser.add_option("-v", "--verbose", action="store_true", dest="bVerbose",
                  help="Be verbose", default=False)

(options, args) = parser.parse_args()


if (options.fnSim == None):
    print >> sys.stderr, 'You must provide the calculated intensity curve'
    sys.exit(1)
if ((options.fnExpNetIntensity) == None and (options.fnSample == None or options.fnBuffer == None)):
    print >> sys.stderr, 'You must provide either sample AND buffer, or privide the experimental net intensity'
    sys.exit(1)


startFit    = options.qmin
endFit      = options.qmax
bVerbose    = options.bVerbose
bUseWeights = options.bUseWeights
bLogScale   = options.bLogScale
xsam1 = []
xbuf1 = []

# Read calculated intensity
[xcal1, ycal1, ycalerr1] = read_xvg(options.fnSim)

# Read experimental net intensity, or sample+buffer
if (options.fnExpNetIntensity != None):
    [xnet1, ynet1, yneterr1] = read_xvg(options.fnExpNetIntensity)
else:
    [xsam1, ysam1, ysamerr1] = read_xvg(options.fnSample)
    [xbuf1, ybuf1, ybuferr1] = read_xvg(options.fnBuffer)
    if (xsam1[0] != xbuf1[0] or xsam1[-1] != xbuf1[-1]):
        print >> sys.stderr, 'q range mismatch between sample and buffer intensity files'

if (options.fitType == FIT_F or options.fitType == FIT_FC):
    bNeedNetIntensity = True
else:
    # We need sample and buffer separately
    bNeedNetIntensity = False
    if (options.fnSample == None or options.fnBuffer == None):
         'With --Aac, you must provide sample AND buffer intensities'

if (options.fnExpNetIntensity == None):
    [xnet1, ynet1, yneterr1] = sampleBufferDiff(xsam1, ysam1, ysamerr1, xbuf1, ybuf1, ybuferr1, alpha=1-options.volSolute)


# Scale x by 10 if experimental files were given in inverse Angstroem
if (options.bExpAngstroem):
    for i in range(len(xnet1)):
        xnet1[i] *= 10
    if (len(xsam1) > 0):
        for i in range(len(xsam1)):
            xsam1[i] *= 10         
    if (len(xbuf1) > 0):
        for i in range(len(xbuf1)):
            xbuf1[i] *= 10


# Get a list x[] where we will fit (within net intensity x-range and within [startfit, endfit]
if (endFit == -1):
    endFit = xnet1[-1]
x = xcal1
while (x[0] < xnet1[0] or x[0] < startFit):
    x = x[1:]
while (x[-1] > xnet1[-1] or x[-1] > endFit):
    x = x[:-1]


# interpolate at N points between Min and Max
N = len(x)
print 'Will fit at %d points between q = %g and %g ' % (N, x[0], x[-1])
y    = numpy.interp(x, xcal1, ycal1)
# yerr = numpy.interp(x, xcal1, ycalerr1)

if (len(ycalerr1) > 0):
    yerr = numpy.interp(x, xcal1, ycalerr1)
else:
    yerr = []

if (options.fitType == FIT_F or options.fitType == FIT_FC):
    # Using net intensity in fit
    ynet = numpy.interp(x, xnet1, ynet1)
    if (len(yneterr1) > 0):
        yneterr = numpy.interp(x, xnet1, yneterr1)
    else:
        yneterr = []
elif (options.fitType == FIT_FAC):
    # Using sample and buffer in fit
    ysam = numpy.interp(x, xsam1, ysam1)
    if (len(ysamerr1) > 0):
        ysamerr = numpy.interp(x, xsam1, ysamerr1)
    ybuf = numpy.interp(x, xbuf1, ybuf1)
    if (len(ybuferr1) > 0):
        ybuferr = numpy.interp(x, xbuf1, ybuferr1)


###################################        
# DO THE FIT ######################
###################################
        
if (options.fitType == FIT_FC):
    # Fitting only the absolute scale and offset : I(fit) = f I(net) + c
    pos  = [y[0]/ynet[0], 1.]
    # Initial values
    # print 'Start with %g / %g' % (pos[0], pos[1])
    # args = [y1int, y2int]
    if (bLogScale == False):
        fminOut = fmin_powell(IntensityDiff_fc, pos, args=(y, yerr, ynet, yneterr, bUseWeights), full_output=True, disp=0)
    else:
        fminOut = fmin_powell(logIntensityDiff_fc, pos, args=(y, yerr, ynet, yneterr, bUseWeights), full_output=True, disp=0)
        print 'After: func = %g' % (fminOut[1])
    xopt    = fminOut[0]
    funcopt = fminOut[1]
    f = fminOut[0][0]
    c = fminOut[0][1]        
    
elif (options.fitType == FIT_F):
    # Fitting only the absolute scale : I(fit) = f I(net) 
    lenSAXS = 0
    pos = [y[0]/ynet[0]]
    if (bLogScale == False):
        fminOut = fmin_powell(intensityDiff_SAXS, pos, args=(y, yerr, ynet, yneterr, bUseWeights), full_output=True, disp=0)
    else:
        fminOut = fmin_powell(logIntensityDiff_SAXS, pos, args=(y, yerr, ynet, yneterr, bUseWeights), full_output=True, disp=0)
    # print 'After: func = %g' % (fminOut[1])
    print fminOut
    xopt    = fminOut[0]
    funcopt = fminOut[1]
    f = fminOut[0]
    c = 0
    
elif (options.fitType == FIT_FAC):
    # Fitting the AXES method: I(fit) = f [ I(sam) - a I(buf) + c]
    alpha0   = 1 - options.volSolute
    alphasig = options.volSoluteSig
    pos      = [y[0]/(ysam[0]-ybuf[0]), 1., 1.]
    fminOut = fmin_powell(intensityDiff_fac, pos,
                          args=(y, yerr, ysam, ysamerr, ybuf, ybuferr, bUseWeights, alpha0, alphasig),
                          full_output=True, disp=0)
    print fminOut
    xopt    = fminOut[0]
    funcopt = fminOut[1]
    f = fminOut[0][0]
    a = fminOut[0][1]
    c = fminOut[0][2]


#
# Compute the fitted intensity
#

if (options.fitType == FIT_FC or options.fitType == FIT_F):
    # Calculate the fitted net intensity from ynet
    for i in range(len(ynet1)):
        ynet1[i] = f * ynet1[i] + c
    for i in range(len(yneterr1)):
        yneterr1[i] = f * yneterr1[i]
    chi2 = funcopt
    print 'f    = %12g' % f
    print 'c    = %12g' % c
    print 'chi2 = %12g' % funcopt
    print 'chi  = %12g' % (math.sqrt(funcopt))

elif (options.fitType == FIT_FAC):
    # Calculate the fitted intensity from sample and buffer
    for i in range(len(ynet1)):
        ynet1[i] = f * (ysam1[i] - a*ybuf1[i] + c)
    for i in range(len(yneterr1)):
        yneterr1[i] = math.sqrt( f**2 * (ysamerr1[i]**2 + (a*ybuferr1[i])**2))

    # To get chi2, we need to recalculate chi2 withtin the regularization on alpha
    pos  = [f, a, c]
    args = [y, yerr, ysam, ysamerr, ybuf, ybuferr, bUseWeights, alpha0, alphasig]
    chi2 = intensityDiff_fac(pos, args)
    print 'f    = %12g' % f
    print 'a    = %12g' % a
    print 'c    = %12g' % c
    print 'chi2 = %12g' % (chi2)
    print 'chi  = %12g' % (math.sqrt(chi2))


# Write overlay between calcualated and fitted experimental curve
fp = open(options.fnOverlay, 'w')
write_to_xvg(fp, xcal1, ycal1, ycalerr1)
write_to_xvg(fp, xnet1, ynet1, yneterr1)
fp.close()

# Write only fitted experimental curve
fp = open(options.fnFitted, 'w')
write_to_xvg(fp, xnet1, ynet1, yneterr1)
fp.close()

#
# Compute Chi(free) following Rambo & Tainer
#
if (options.nfree):
    chi_free(x, y, yerr, ynet, yneterr, options.D, options.nfree, bLogScale, bUseWeights)
    

#
# Write the calculated I(q) at the experimental q-values via interpolation
#
# Get q points of experiment, but only within the calculated q-range
#
xExpForInterp = xnet1
while (xExpForInterp[0] < xcal1[0]):
    xExpForInterp = xExpForInterp[1:]
while (xExpForInterp[-1] > xcal1[-1]):
    xExpForInterp = xExpForInterp[:-1]

if (options.fnOverlayExpQ != None or options.fnIcalcExpQ != None):
    yCalcAtExpQ    = numpy.interp(xExpForInterp, xcal1, ycal1)
    yCalcErrAtExpQ = numpy.interp(xExpForInterp, xcal1, ycalerr1)
    
    yExpAtExpQ     = numpy.interp(xExpForInterp, xnet1, ynet1)
    yExpErrAtExpQ  = numpy.interp(xExpForInterp, xnet1, yneterr1)

if (options.fnOverlayExpQ != None):
    fp = open(options.fnOverlayExpQ, 'w')
    print >> fp, '#\n# Fitted experimental I(q) and calculated I(q) in multi-column format\n#'
    print >> fp, '# Columns 1 though 5 are:'
    print >> fp, '#   1) Q (from experimental data)'
    print >> fp, '#   2) Fitted experimental I(q)'
    print >> fp, '#   3) Experimental error after fit'
    print >> fp, '#   4) Calculated/model I(q) (no fit)'
    print >> fp, '#   5) Error of calcucated I(q)\n#'
    for i in range(len(xExpForInterp)):
        print >> fp, '%12g    %12g %12g    %12g %12g' % (xExpForInterp[i], yExpAtExpQ[i], yExpErrAtExpQ[i], yCalcAtExpQ[i], yCalcErrAtExpQ[i])
    fp.close()

if (options.fnIcalcExpQ != None):
    fp = open(options.fnIcalcExpQ, 'w')
    print >> fp, '# Calculated I(q) in this file is the WAXSiS result linearly interpolated to the experimental q values.'
    print >> fp, '# Columns: q, I(q), sigma(q)'
    write_to_xvg(fp, xExpForInterp, yCalcAtExpQ, yCalcErrAtExpQ)
    fp.close()


    
