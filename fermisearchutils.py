"""
fermisearchutils.py contains the functions to search for periodic sources (BW/RBs) within 3-sigma of Fermi position.
Input: 
PTF_*.ctlg
IPAC SExtractor catalogs with zero point correction
Output:
A list of periodic candidates and their light curve plots.
    
Sumin Tang, May 2014
    
"""


import sys
import os
import pylab
import numpy as np
import pyfits
from glob import glob
import matplotlib.pyplot as plt
import matplotlib
from scipy.stats import percentileofscore
import time

# Import Lomb-Scargle w/ polynomial detrending, treating meas. error
from lomb_scargle_cc_order import lomb as lombo
# Import multi-harmonic version
from lomb_scargle_refine import lomb as lombr

# plotting parameters
matplotlib.rcParams.update({'font.size': 14})
matplotlib.rc('axes', linewidth=1.0)
matplotlib.rcParams['figure.figsize'] = 8, 10

# === ctlgcombine: build a big table including all sources within 3*sma
def ctlgcombine(key = 'f01', ferminame = 'J1544.5-1126', ra0 = 236.139, dec0 = -11.446, sma = 5.4):
    BASE_DIR = '/scr4/stang/fermi/data/' + ferminame + '/'
    fnames = glob(BASE_DIR + 'PTF*_' + key + '_*.ctlg')
    outfile = key + '_combined.tbl'
    g = open(BASE_DIR + outfile, 'w')
    g.write('# ra, dec, mag, magerr, MHJD\n')
    for fname in fnames:
        hdulist = pyfits.open(fname)
        HJD = hdulist[2].header['HJD']
        MHJD = HJD - 2400000.5

        phtcalex = hdulist[2].header['PHTCALEX']
        #phtcalfl = hdulist[2].header['PHTCALFL']
        #print phtcalex, phtcalfl

        if (phtcalex == 1):
            data = hdulist[1].data
            ra = data['X_WORLD']
            dec = data['Y_WORLD']
            magzero = data['ZEROPOINT']
            mag = data['MAG_AUTO']#[:,2] 
            magerr = data['MAGERR_AUTO']

            dsep = ((ra - ra0)**2*np.cos(dec0*np.pi/180.)**2 + (dec - dec0)**2)**0.5*60. # in arcmin
            ix = np.where(dsep <= 3*sma)[0]
            for i in ix:
                g.write('%f  %f  %.3f  %.3f  %f\n' % (ra[i], dec[i], mag[i] + magzero[i], magerr[i], MHJD))

    g.close()
    print 'Catalogs combined for %s; output: %s' % (ferminame, outfile)

# === match sources and split light curves
def lcsplit(infile = 'f01_combined.tbl', ferminame = 'J1544.5-1126'):
    BASE_DIR = '/scr4/stang/fermi/data/' + ferminame + '/'
    outdir = BASE_DIR + 'lc/'
    os.system("mkdir -p %s" % outdir) # mkdir if not exist

    key = infile[0:3]
    listfile = 'lc_' + key + '.list'
    g = open(outdir + listfile, 'w')
    g.write('# starid, ra, dec, ndet, mag\n')

    # load the data
    ra, dec, mag, magerr, MHJD = np.loadtxt(BASE_DIR + infile, unpack = True, skiprows = 1)
    ix = np.argsort(ra)
    ra2 = ra[ix]
    dec2 = dec[ix]
    mag2 = mag[ix]
    magerr2 = magerr[ix]
    mjd2 = MHJD[ix]

    ra3 = ra2[:]
    dec3 = dec2[:]   
    nvar = 0
    for i in range(len(ra3)):
        ra1 = ra3[i]
        dec1 = dec3[i]    
        dist = (((ra2-ra1)*np.cos(dec1*np.pi/180.))**2 + (dec2-dec1)**2 )**0.5*3600.0 # angular separation in arcsec
        ix = np.where(dist<2)[0]
        nmatch = len(ix)

        if nmatch > 5: # more than 5 detections
            nvar = nvar + 1
            myra = ra2[ix]
            mydec = dec2[ix]        
            mymag = mag2[ix]
            mymagerr = magerr2[ix]
            mymjd = mjd2[ix]

            g.write('%i  %f  %f  %i  %.3f\n' % (nvar, np.median(myra), np.median(mydec), len(myra), np.median(mymag)))

            h = open(outdir  +  'lc_' + key + '_star_' + np.str(nvar) + '.txt', 'w')
            h.write("# ra, dec, mag, magerr, mhjd\n")

            for j in range(nmatch):
                h.write('%f  %f  %.3f  %.3f  %f \n' % (myra[j], mydec[j], mymag[j], mymagerr[j], mymjd[j])
)
            h.close()       

        # drop the matched ones from ra2
        ix2 = np.where(dist>=3)[0]
        ra2 = ra2[ix2]
        dec2 = dec2[ix2]
        mag2 = mag2[ix2]
        magerr2 = magerr2[ix2]
        mjd2 = mjd2[ix2]

    g.close()
    print "Total number of objects: %i" % nvar
    return nvar


# === combine two bands
def combgnr(file1 = 'lc_f01.list', file2 = 'lc_f02.list', ferminame = 'J1544.5-1126'):
    BASE_DIR = '/scr4/stang/fermi/data/' + ferminame + '/'
    indir = BASE_DIR + 'lc/'
    outfile = 'lc_f01-f02.list'
    g = open(indir + outfile, 'w')
    g.write('# gid, gra, gdec, gndet, gmag, rid, rra, rdec, rndet, rmag\n')
    
    num_lines1 = sum(1 for line in open(indir + file1))
    num_lines2 = sum(1 for line in open(indir + file2))
    if num_lines1>10:
        id1, ra1, dec1, ndet1, mag1 = np.loadtxt(indir + file1, unpack = True, skiprows = 1)
    if num_lines2>10:
        id2, ra2, dec2, ndet2, mag2 = np.loadtxt(indir + file2, unpack = True, skiprows = 1)

    # match the two lists if f01/g data is not empty
    if num_lines1>10:
        # print all f01 list
        for i in range(len(id1)):
            ra3 = ra1[i]
            dec3 = dec1[i]    
            dist = (((ra2-ra3)*np.cos(dec3*np.pi/180.))**2 + (dec2-dec3)**2 )**0.5*3600.0 # angular separation in arcsec
            ix = np.where(dist<1)[0]
            nmatch = len(ix)

            if nmatch == 1:
                g.write('%i  %f  %f  %i  %.3f  %i  %f  %f  %i  %.3f\n' % (id1[i], ra1[i], dec1[i], ndet1[i], mag1[i], id2[ix], ra2[ix], dec2[ix], ndet2[ix], mag2[ix]))
            else:
                g.write('%i  %f  %f  %i  %.3f  %i  %f  %f  %i  %.3f\n' % (id1[i], ra1[i], dec1[i], ndet1[i], mag1[i], 0, 0, 0, 0, 0))

        # print all f02 list which are not matched with f01
        for i in range(len(id2)):
            ra3 = ra2[i]
            dec3 = dec2[i]    
            dist = (((ra1-ra3)*np.cos(dec3*np.pi/180.))**2 + (dec1-dec3)**2 )**0.5*3600.0 # angular separation in arcsec
            ix = np.where(dist<1)[0]
            nmatch = len(ix)

            if nmatch == 0:
                g.write('%i  %f  %f  %i  %.3f  %i  %f  %f  %i  %.3f\n' % (0, 0, 0, 0, 0, id2[i], ra2[i], dec2[i], ndet2[i], mag2[i]))

    # if f01/g is empty
    else:
        for i in range(len(id2)):
            g.write('%i  %f  %f  %i  %.3f  %i  %f  %f  %i  %.3f\n' % (0, 0, 0, 0, 0, id2[i], ra2[i], dec2[i], ndet2[i], mag2[i]))

    g.close()
    print "f01 and f02 lists combined: %s%s" % (indir, outfile)



# === relative photometry
def relphot(key = 'f01', ferminame = 'J1544.5-1126'):
    BASE_DIR = '/scr4/stang/fermi/data/' + ferminame + '/'
    inlist = BASE_DIR + 'lc/' + 'lc_' + key + '.list'
    intbl = BASE_DIR + key + '_combined.tbl'
    outfile = BASE_DIR + 'lc/' + key + 'relphot.list'
    g = open(outfile, 'w')
    g.write('# MHJD, dmag, dmagerr\n')

    # load the data
    ra, dec, mag, magerr, MHJD = np.loadtxt(intbl, unpack = True, skiprows = 1)
    id1, ra1, dec1, ndet1, mag1 = np.loadtxt(inlist, unpack = True, skiprows = 1)
    
    # select objects in 16-18th mag
    ix = np.where((mag1>16) & (mag1<18))[0]
    ra2 = ra1[ix]
    dec2 = dec1[ix]
    mag2 = mag1[ix]
    nstar = len(ra2)
    print 'Number of stars between 16-18 mag: %i' % nstar
    
    umjd = np.unique(MHJD)
    print "Number of images: %i" % len(umjd)
    for mjd in umjd:
        iy = np.where(abs(MHJD-mjd)<1./24./60.)[0]
        ra3 = ra[iy]
        dec3 = dec[iy]
        mag3 = mag[iy]
        
        dmag = []
        for i in range(nstar):
            myra = ra2[i]
            mydec = dec2[i]
            mymag = mag2[i]
            
            dist = (((myra-ra3)*np.cos(mydec*np.pi/180.))**2 + (mydec-dec3)**2 )**0.5*3600.0 # angular separation in arcsec
            iz = np.where(dist<1)[0]
            nmatch = len(iz)
            if nmatch == 1:
                mydmag = mymag - mag3[iz]
                dmag.append(mydmag)

        dmag2 = np.median(np.asarray(dmag).flatten())
        dmag2err = np.std(np.asarray(dmag).flatten())
        g.write('%f  %.3f  %.3f\n' % (mjd, dmag2, dmag2err))
    g.close()
            
        
    # apply to all the lc files in the given filter
    relmjd, reldmag, reldmagerr = np.loadtxt(outfile, unpack = True, skiprows = 1)
            
    fnames = glob(BASE_DIR + 'lc/lc_' + key + '_star*.txt')
    for fname in fnames:
        ra, dec, mag, magerr, mhjd = np.loadtxt(fname, unpack = True, skiprows = 1)
        
        fname2 = fname[:-4] + '_relphot.dat'
        h = open(fname2, 'w')
        h.write('# ra, dec, mag, magerr, mhjd, magrel, magrelerr\n')
        
        for j in range(len(ra)):
            mymjd = mhjd[j]
            ix = np.where(abs(relmjd-mymjd)<1./24./60.)[0]
            if len(ix) == 1:
                mymagrel = mag[j] + reldmag[ix]
                mymagrelerr = (magerr[j]**2 + reldmagerr[ix]**2)**0.5
                h.write('%f  %f  %.3f  %.3f  %f  %.3f  %.3f\n' % (ra[j], dec[j], mag[j], magerr[j], mymjd, mymagrel, mymagrelerr))
        
        h.close()






# === calculate lc statistics, including periodic info
def lcstat(infile = 'lc_f01-f02.list', ferminame = 'J1544.5-1126'):
    BASE_DIR = '/scr4/stang/fermi/data/' + ferminame + '/'
    indir = BASE_DIR + 'lc/'
    id1, ra1, dec1, ndet1, gmag, id2, ra2, dec2, ndet2, rmag = np.loadtxt(indir + infile, unpack = True, skiprows = 1)

    outfile = indir + 'lc_f01-f02_stats.txt'
    g = open(outfile, 'w')
    g.write('# gid, gra, gdec, gndet, glcmedian, grawlcrms, glcrms, glcerr, grchi2, glcamp0, glcamp2, glcamp5, gS, gK, gJB, gBC, gbestpsd, gbestp, gbestt0, rid, rra, rdec, rndet, rlcmedian, rrawlcrms, rlcrms, rlcerr, rrchi2, rlcamp0, rlcamp2, rlcamp5, rS, rK, rJB, rBC, rbestpsd, rbestp, rbestt0\n')

    for i in range(len(id1)):
        # === f01
        myid = int(id1[i])
        key = 'f01'
        nall = ndet1[i]
        if (myid>0) & (nall>1):
            fname = 'lc_' + key + '_star_' + str(myid) + '_relphot.dat' 
            ra, dec, magraw, magrawerr, mhjd, mag, magerr = np.loadtxt(indir + fname, unpack = True, skiprows = 1)
            nall = len(mag)
 
            lcmedian = np.median(mag)
            rawlcrms = np.std(magraw)
            lcrms = np.std(mag)
            lcerr = np.median(magerr)
            dmag = mag - lcmedian
            rchi2 = np.sum(dmag**2/magerr**2)/(nall-1)
            lcamp0 = mag.max() - mag.min()
            lcamp2 = np.percentile(mag, 98) - np.percentile(mag, 2)
            lcamp5 = np.percentile(mag, 95) - np.percentile(mag, 5)
            # skewness, kurtosis and JB test
            S = sum(dmag**3.) / nall / lcrms**3. 
            K = sum(dmag**4.) / nall / lcrms**4. - 3.
            JB = nall/6 *(S**2. + ((K-3.)**2)/4)
            BC = (S**2 + 1)/(K + 3*(nall-1)**2/(nall-2)/(nall-3))
            # http://www.i3.psychologie.uni-wuerzburg.de/fileadmin/06020300/user_upload/Janczyk/pfister_et_al_2013_bimodality.pdf
            # BC>0.555 => bimodal distribution
            # didn't work

            # period search
            x = mhjd
            y = mag
            dy = magerr
            Xmax = x.max()-x.min()
            f0 = 1./3.; df = 0.1/Xmax; fe = 50. # from 0.02 to 3 days
            numf = int((fe-f0)/df)
            freqin = f0 + df*np.arange(numf,dtype='float64')
            # Fit a constant plus sinusoid at every frequency
            detrend_order=0
            psd = lombo(x,y,dy,f0,df,numf,detrend_order=detrend_order)
            bestpsd = psd.max()
            bestp = 1./freqin[psd.argmax()]
        else:
            lcmedian, rawlcrms, lcrms, lcerr, rchi2, lcamp0, lcamp2, lcamp5, S, K, JB, BC, bestpsd, bestp = (0,)*14

        # === f02
        myid = int(id2[i])
        key = 'f02'
        nall = ndet2[i]
        if (myid>0) & (nall>1):
            fname = 'lc_' + key + '_star_' + str(myid) + '_relphot.dat' 
            ra, dec, magraw, magrawerr, mhjd, mag, magerr = np.loadtxt(indir + fname, unpack = True, skiprows = 1)
            nall = len(mag)
 
            rlcmedian = np.median(mag)
            rrawlcrms = np.std(magraw)
            rlcrms = np.std(mag)
            rlcerr = np.median(magerr)
            dmag = mag - lcmedian
            rrchi2 = np.sum(dmag**2/magerr**2)/(nall-1)
            rlcamp0 = mag.max() - mag.min()
            rlcamp2 = np.percentile(mag, 98) - np.percentile(mag, 2)
            rlcamp5 = np.percentile(mag, 95) - np.percentile(mag, 5)
            # skewness, kurtosis and JB test
            rS = sum(dmag**3.) / nall / rlcrms**3. 
            rK = sum(dmag**4.) / nall / rlcrms**4. 
            rJB = nall/6 *(rS**2. + ((rK-3.)**2)/4)
            rBC = (rS**2 + 1)/(rK + 3*(nall-1)**2/(nall-2)/(nall-3))


            # period search
            x = mhjd
            y = mag
            dy = magerr
            Xmax = x.max()-x.min()
            f0 = 1./3.; df = 0.1/Xmax; fe = 50. # from 0.02 to 3 days
            numf = int((fe-f0)/df)
            freqin = f0 + df*np.arange(numf,dtype='float64')
            # Fit a constant plus sinusoid at every frequency
            detrend_order=0
            psd = lombo(x,y,dy,f0,df,numf,detrend_order=detrend_order)
            rbestpsd = psd.max()
            rbestp = 1./freqin[psd.argmax()]
            # constant plus sinusoid and 7 harmonics at every frequency; takes 10x time of lombo
            """
            tic = time.time()
            psdr, res2 = lombr(x,y,dy,f0,df,numf,detrend_order=detrend_order)
            print 'Time for lombr: %f' % (time.time() - tic)
            # p1 = 1./res2['freq']
            rbestt0 = res2['time0']
            """
        else:
            rlcmedian, rrawlcrms, rlcrms, rlcerr, rrchi2, rlcamp0, rlcamp2, rlcamp5, rS, rK, rJB, rBC, rbestpsd, rbestp = (0,)*14

        g.write('%i  %f  %f  %i  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %f  %i  %f  %f  %i  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %f\n' % (id1[i], ra1[i], dec1[i], ndet1[i], lcmedian, rawlcrms, lcrms, lcerr, rchi2, lcamp0, lcamp2, lcamp5, S, K, JB, BC, bestpsd, bestp, id2[i], ra2[i], dec2[i], ndet2[i], rlcmedian, rrawlcrms, rlcrms, rlcerr, rrchi2, rlcamp0, rlcamp2, rlcamp5, rS, rK, rJB, rBC, rbestpsd, rbestp))

    g.close()
    print 'light curve statistics calculated: %s' % outfile


# === Plot statistics 
def statplot(infile = 'lc_f01-f02_stats.txt', ferminame = 'J1544.5-1126', ra0 = 236.139, dec0 = -11.446, sma = 5.4):
    BASE_DIR = '/scr4/stang/fermi/data/' + ferminame + '/'
    outdir = BASE_DIR + 'plots/'
    os.system("mkdir -p %s" % outdir) # mkdir if not exist
    indir = BASE_DIR + 'lc/'
    id1, ra1, dec1, ndet1, lcmedian, rawlcrms, lcrms, lcerr, rchi2, lcamp0, lcamp2, lcamp5, S, K, JB, BC, bestpsd, bestp, id2, ra2, dec2, ndet2, rlcmedian, rrawlcrms, rlcrms, rlcerr, rrchi2, rlcamp0, rlcamp2, rlcamp5, rS, rK, rJB, rBC, rbestpsd, rbestp = np.loadtxt(indir + infile, unpack = True, skiprows = 1)
    
    os.system("rm -f %s*.eps" % outdir)
    
    figname = ferminame + '_lcmedian-rms.eps'
    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax.plot(lcmedian, rawlcrms, 'g.')
    ax.plot(rlcmedian, rrawlcrms, 'r.')
    ax.set_xlabel('median mag', fontsize=16)
    ax.set_ylabel('rms (SEcat)', fontsize=16)
    ax.set_xlim([12, 22])
    ax.set_ylim([0, 0.8])
    ax.set_title(ferminame + ', ra=%.4f, dec=%.4f' % (ra0, dec0))
    ax = fig.add_subplot(212)
    ax.plot(lcmedian, lcrms, 'g.')
    ax.plot(rlcmedian, rlcrms, 'r.')
    ax.set_xlabel('median mag', fontsize=16)
    ax.set_ylabel('rms (relphot)', fontsize=16)
    ax.set_xlim([12, 22])
    ax.set_ylim([0, 0.4])
    # ax.set_ylim(ax.get_ylim()[::-1])
    pylab.savefig(outdir + figname)
    plt.close()

    figname = ferminame + '_psd-period.eps'
    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax.plot(bestp, bestpsd, 'g.')
    ax.plot(rbestp, rbestpsd, 'r.')
    pct = 95.
    prange = [0.01, 3.]
    ax.plot(prange, [np.percentile(bestpsd, pct), np.percentile(bestpsd, pct)], 'g--')
    ax.plot(prange, [np.percentile(rbestpsd, pct), np.percentile(rbestpsd, pct)], 'r--')
    ax.set_xlabel('Period (days)', fontsize=18)
    ax.set_ylabel('Power', fontsize=18)
    ax.set_xlim(prange)
    ax.set_xscale('log')
    #ax.set_ylim([0, 1])
    # ax.set_ylim(ax.get_ylim()[::-1])
    ax.set_title(ferminame + ', ra=%.4f, dec=%.4f' % (ra0, dec0))
    ax = fig.add_subplot(223)
    ax.plot(bestpsd, rbestpsd, 'k.')
    ax.set_xlabel('Power (g)', fontsize=18)
    ax.set_ylabel('Power (r)', fontsize=18)
    ax = fig.add_subplot(224)
    ax.plot(bestp, rbestp, 'k.')
    ax.set_xlabel('Period (g)', fontsize=18)
    ax.set_ylabel('Period (r)', fontsize=18)
    ax.yaxis.set_label_position("right")
    ax.plot(prange, prange, 'b--')
    ax.set_xlim(prange)
    ax.set_ylim(prange)
    ax.set_xscale('log')
    ax.set_yscale('log')
    pylab.savefig(outdir + figname)
    plt.close()





# === search for periodic objects
def varselect(infile = 'lc_f01-f02_stats.txt', ferminame = 'J1544.5-1126', ra0 = 236.139, dec0 = -11.446, sma = 5.4):
    BASE_DIR = '/scr4/stang/fermi/data/' + ferminame + '/'
    outdir = BASE_DIR + 'plots/'
    os.system("mkdir -p %s" % outdir) # mkdir if not exist
    
    # output list
    outlist = BASE_DIR + ferminame +'_periodic_candidates.list'
    h = open(outlist, 'w')
    h.write('# ra, dec, gid, rid, period(days)\n')
    
    
    indir = BASE_DIR + 'lc/'
    id1, ra1, dec1, ndet1, lcmedian, rawlcrms, lcrms, lcerr, rchi2, lcamp0, lcamp2, lcamp5, S, K, JB, BC, bestpsd, bestp, id2, ra2, dec2, ndet2, rlcmedian, rrawlcrms, rlcrms, rlcerr, rrchi2, rlcamp0, rlcamp2, rlcamp5, rS, rK, rJB, rBC, rbestpsd, rbestp = np.loadtxt(indir + infile, unpack = True, skiprows = 1)
    
    pct1 = 98.
    pct2 = 90.
    psdth1 = np.percentile(bestpsd[np.nonzero(bestpsd)], pct1)
    rpsdth1 = np.percentile(rbestpsd[np.nonzero(rbestpsd)], pct1)
    psdth2 = np.percentile(bestpsd[np.nonzero(bestpsd)], pct2)
    rpsdth2 = np.percentile(rbestpsd[np.nonzero(rbestpsd)], pct2)
    
    # Select candidates:
    # 1. only if ndet1>10 & ndet2>20 & rrchi2>2
    # 2. R-band psdscore>0.9, rrchi2>2
    # 3. at least one psdscore>0.98 
    # 4. P not near 1.0 or 0.5 days
    ixcand = np.where( (ndet1>10) & (ndet2>20) & (rrchi2>2.) & (rbestpsd>rpsdth2) & ( (bestpsd>psdth1) | (rbestpsd>rpsdth1)) & (abs(rbestp - 1.)>0.02) & (abs(rbestp - 0.5)>0.01))[0]
   
    print 'Number of candidates: %i (out of %i)' % (len(ixcand), len(ra1))     
    
    for i in ixcand:        
        # pick period
        psdscore1 = percentileofscore(bestpsd[np.nonzero(bestpsd)], bestpsd[i])
        psdscore2 = percentileofscore(rbestpsd[np.nonzero(rbestpsd)], rbestpsd[i])
        if psdscore1>psdscore2:
            myp = bestp[i]
            mypsdscore = psdscore1
        else:
            myp = rbestp[i]
            mypsdscore = psdscore2        
        
        # plot the light curves
        figname = ferminame + '_g-%i_r-%i.eps' % (id1[i], id2[i])
        fig = plt.figure()
        ax = fig.add_subplot(211)
        
        # g-band
        gid = int(id1[i])
        gkey = 'f01'
        mjdall = []
        if gid>0:
            fname = 'lc_' + gkey + '_star_' + str(gid) + '_relphot.dat' 
            gra, gdec, gmagraw, gmagrawerr, gmhjd, gmag, gmagerr = np.loadtxt(indir + fname, unpack = True, skiprows = 1)
            plt.errorbar(gmhjd, gmag, yerr = gmagerr, fmt = 'g.', markersize = 3)  
            myra = np.median(gra)
            mydec = np.median(gdec)
            mjdall.append(gmhjd)
        
        # R-band
        rid = int(id2[i])
        rkey = 'f02'
        if rid>0:
            fname = 'lc_' + rkey + '_star_' + str(rid) + '_relphot.dat' 
            rra, rdec, rmagraw, rmagrawerr, rmhjd, rmag, rmagerr = np.loadtxt(indir + fname, unpack = True, skiprows = 1)
            plt.errorbar(rmhjd, rmag, yerr = rmagerr, fmt = 'r.', markersize = 3)
            myra = np.median(rra)
            mydec = np.median(rdec)
            mjdall.append(rmhjd)
                
        dsep = (((myra - ra0)*np.cos(dec0*np.pi/180.))**2 + (mydec-dec0)**2 )**0.5*60.0 # angular separation in arcmin
        
        ax.set_xlabel('MJD')
        mjdall2 = np.asarray(np.hstack(mjdall))
        mjdmin = np.percentile(mjdall2, 10) - 10.
        mjdmax = np.percentile(mjdall2, 90) + 10.                
        ax.set_xlim([mjdmin, mjdmax])
        ax.set_ylabel('Mag')
        ax.set_ylim(ax.get_ylim()[::-1])
        ax.set_title(figname[:-4] + "; ra=%.5f, dec=%.5f, dsep=%.2f'" % (myra, mydec, dsep), fontsize=12)
    
        ax = fig.add_subplot(212)
        msize = 10
        if gid>0:
            gphase = np.mod(gmhjd, myp)/myp
            ax.plot(gphase, gmag-np.mean(gmag), 'g.', markersize=msize)
            ax.plot(gphase+1., gmag-np.mean(gmag), 'g.', markersize=msize)
        if rid>0:
            rphase = np.mod(rmhjd, myp)/myp
            ax.plot(rphase, rmag-np.mean(rmag), 'r.', markersize=msize)
            ax.plot(rphase+1., rmag-np.mean(rmag), 'r.', markersize=msize)

        ax.set_xlabel('Phase; P = %f hr, PSD score = %.2f (g), %.2f (R)' % (myp*24., psdscore1, psdscore2))
        ax.set_ylabel('Mag')
        ax.set_ylim(ax.get_ylim()[::-1])
        ax.set_xlim([0, 1.5])

        pylab.savefig(outdir + figname)
        plt.close()

        h.write('%f  %f  %i  %i  %f\n' % (myra, mydec, gid, rid, myp))
    h.close()
    # combine the plots to a single pdf
    if len(ixcand)>1:
        plotdir = '/scr4/stang/fermi/data/plots/'
        os.system('/usr/bin/gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=%s%s_allplots.pdf %s%s*.eps' %(plotdir, ferminame, outdir, ferminame))


if __name__ == "__main__":
    tic = time.time()

    # ferminame = 'J0212.1+5318'
    # ra0 = 33.040
    # dec0 = 53.305
    # sma = 3.0 # sma in arcmin

    ferminame = 'J1544.5-1126'  
    ra0 = 236.139  
    dec0 = -11.446    
    sma = 5.4 

    BASE_DIR = '/scr4/stang/fermi/data/' + ferminame + '/'
    
    """
    ctlgcombine(key = 'f01')
    ctlgcombine(key = 'f02')
    lcsplit(infile = 'f01_combined.tbl')
    lcsplit(infile = 'f02_combined.tbl')
    combgnr(file1 = 'lc_f01.list', file2 = 'lc_f02.list')

    relphot(key = 'f01')
    relphot(key = 'f02')

    lcstat(infile = 'lc_f01-f02.list')
    """
    statplot(infile = 'lc_f01-f02_stats.txt')
    varselect(infile = 'lc_f01-f02_stats.txt')

    print 'Time elapsed: %.3f' % (time.time()-tic)
