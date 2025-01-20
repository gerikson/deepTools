import pytest
from deeptools.hp import r_bamcoverage
import os.path
import filecmp
from os import unlink

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
BAMFILE_A = ROOT + "testA.bam"
BAMFILE_B = ROOT + "testB.bam"


def test_r_bamcoverage():
    bamifile = BAMFILE_A
    ofile = ROOT + "output.bedgraph"
    ofiletype = "bedgraph"
    norm = "None"
    effectivegenomesize = 0
    scalefactor = 1.0
    mnase = False
    offset = [1, -1]  # Adjusted to a regular Python list
    extendreads = 0
    centerreads = False
    filterrnastrand = "none"
    blacklist = ""
    ignorechr = []
    skipnoncovregions = False
    smoothlength = 0
    binsize = 50
    ignoreduplicates = False
    minmappingquality = 0
    samflaginclude = 0
    samflagexclude = 0
    minfraglen = 0
    maxfraglen = 0
    nproc = 1
    regions = []
    verbose = False


    # Call the Rust function
    r_bamcoverage(
        bamifile, ofile, ofiletype, norm, effectivegenomesize, scalefactor, mnase, offset, extendreads, 
        centerreads, filterrnastrand, blacklist, ignorechr, skipnoncovregions, smoothlength, 
        binsize, ignoreduplicates, minmappingquality, samflaginclude, samflagexclude, minfraglen, maxfraglen, 
        nproc, regions, verbose
    )
    # Add assertions to verify the expected behavior
    expected = ['3R\t0\t100\t0\n', '3R\t100\t200\t1\n', 'chr_cigar\t0\t50\t1\n', 'chr_cigar\t50\t200\t0\n']
    _foo = open(ofile, 'r')
    resp = _foo.readlines()
    _foo.close()
    assert f"{resp}" == f"{expected}", f"{resp} != {expected}"
    unlink(ofile)

def test_r_bamcoverage_RPKM():
    bamifile = BAMFILE_A
    ofile = ROOT + "output.bedgraph"
    ofiletype = "bedgraph"
    norm = "RPKM"
    effectivegenomesize = 0
    scalefactor = 1.0
    mnase = False
    offset = [1, -1]  # Adjusted to a regular Python list
    extendreads = 0
    centerreads = False
    filterrnastrand = "none"
    blacklist = ""
    ignorechr = []
    skipnoncovregions = False
    smoothlength = 0
    binsize = 50
    ignoreduplicates = False
    minmappingquality = 0
    samflaginclude = 0
    samflagexclude = 0
    minfraglen = 0
    maxfraglen = 0
    nproc = 1
    regions = []
    verbose = False


    # Call the Rust function
    r_bamcoverage(
        bamifile, ofile, ofiletype, norm, effectivegenomesize, scalefactor, mnase, offset, extendreads, 
        centerreads, filterrnastrand, blacklist, ignorechr, skipnoncovregions, smoothlength, 
        binsize, ignoreduplicates, minmappingquality, samflaginclude, samflagexclude, minfraglen, maxfraglen, 
        nproc, regions, verbose
    )
    # Add assertions to verify the expected behavior
    expected = ['3R\t0\t100\t0\n', '3R\t100\t200\t6666666.5\n', 'chr_cigar\t0\t50\t6666666.5\n', 'chr_cigar\t50\t200\t0\n']
    _foo = open(ofile, 'r')
    resp = _foo.readlines()
    _foo.close()
    assert f"{resp}" == f"{expected}", f"{resp} != {expected}"
    unlink(ofile)

def test_r_bamcoverage_RPKM_with_effectivegenomesize():
    bamifile = BAMFILE_B
    ofile = ROOT + "testB.bedgraph"
    ofiletype = "bedgraph"
    norm = "RPKM"
    effectivegenomesize = 142573017
    scalefactor = 1.0
    mnase = False
    offset = [1, -1]  # Adjusted to a regular Python list
    extendreads = 0
    centerreads = False
    filterrnastrand = "none"
    blacklist = ""
    ignorechr = []
    skipnoncovregions = False
    smoothlength = 0
    binsize = 50
    ignoreduplicates = False
    minmappingquality = 0
    samflaginclude = 0
    samflagexclude = 0
    minfraglen = 0
    maxfraglen = 0
    nproc = 1
    regions = []
    verbose = False


    # Call the Rust function
    r_bamcoverage(
        bamifile, ofile, ofiletype, norm, effectivegenomesize, scalefactor, mnase, offset, extendreads, 
        centerreads, filterrnastrand, blacklist, ignorechr, skipnoncovregions, smoothlength, 
        binsize, ignoreduplicates, minmappingquality, samflaginclude, samflagexclude, minfraglen, maxfraglen, 
        nproc, regions, verbose
    )
    # Add assertions to verify the expected behavior
    expected = ['3R\t0\t50\t0\n', '3R\t50\t150\t5000000\n', '3R\t150\t200\t10000000\n']
    _foo = open(ofile, 'r')
    resp = _foo.readlines()
    _foo.close()
    assert f"{resp}" == f"{expected}", f"{resp} != {expected}"
    unlink(ofile)