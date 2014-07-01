"""
fermisearch.py is a pipeline to automatically calibrate the brightness of stars in PTF images, and search for periodic variable stars (BW/RBs candidates) within 3-sigma of Fermi LAT position.
Input: 
PTF_*.ctlg
IPAC SExtractor catalogs with zero point correction
Output:
A list of periodic candidates and their light curve plots
    
Sumin Tang, May 2014
    
"""


from fermisearchutils import *

tic0 = time.time()


infile = '/scr4/stang/fermi/ptfcoverage/fermigoodcoverage.txt'
dt = [('name', '|S20'), ('ra', 'f'), ('dec', 'f'), ('l', 'f'), ('b', 'f'), ('sma', 'f'), ('fieldid', 'i'), ('ccdid', 'i'), ('nexpr', 'i'), ('nexpg', 'i')]
data0 = np.loadtxt(infile, dtype = dt, unpack = False, skiprows = 2)

# for sma<=10 ones
ix = np.where(data0['sma']<=10.)[0]
data = data0[ix]

fldnames = np.unique(data['name'])
ufid, indices = np.unique(data['fieldid'], return_index=True)
print "Number of Fermi fields: %i" % len(fldnames)
print "Number of PTF/iPTF images in the R-band: %i" % sum(data['nexpr'][indices])
print "Number of PTF/iPTF images in the g-band: %i" % sum(data['nexpg'][indices])
print "Number of PTF/iPTF CCD-images in the R-band: %i" % sum(data['nexpr'])
print "Number of PTF/iPTF CCD-images in the g-band: %i" % sum(data['nexpg'])

for i in range(len(fldnames)):
    ferminame = fldnames[i]
    #for ferminame in fldnames:
    ix = np.where(data['name'] == ferminame)[0]
    print len(ix)
    ra0 = data['ra'][ix[0]]
    dec0 = data['dec'][ix[0]]
    sma = data['sma'][ix[0]]
    
    print "=== Start runing Fermi periodic search for %s; sma=%.1f ===" % (ferminame, sma)
    tic = time.time()
    
    ctlgcombine(key = 'f01', ferminame = ferminame, ra0 = ra0, dec0 = dec0, sma = sma)    
    ctlgcombine(key = 'f02', ferminame = ferminame, ra0 = ra0, dec0 = dec0, sma = sma)
    
    nvar1 = lcsplit(infile = 'f01_combined.tbl', ferminame = ferminame)
    nvar2 = lcsplit(infile = 'f02_combined.tbl', ferminame = ferminame)
    
    if (nvar1>50) & (nvar2>50):
        combgnr(file1 = 'lc_f01.list', file2 = 'lc_f02.list', ferminame = ferminame)
        relphot(key = 'f01', ferminame = ferminame)
        relphot(key = 'f02', ferminame = ferminame)
        lcstat(infile = 'lc_f01-f02.list', ferminame = ferminame)
        statplot(infile = 'lc_f01-f02_stats.txt', ferminame = ferminame, ra0 = ra0, dec0 = dec0, sma = sma)
        varselect(infile = 'lc_f01-f02_stats.txt', ferminame = ferminame, ra0 = ra0, dec0 = dec0, sma = sma)
    else:
        print "Missing data in at least 1 band! Terminating the process, and move to the next field."

    print 'Time elapsed for field %s: %.2f' % (ferminame, time.time()-tic)


print 'Total time elapsed: %.2f' % (time.time()-tic0)
