import pytest
from deeptools.bamCompare2 import r_bamcompare
import os.path
import filecmp
from os import unlink

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
BAMFILE_A = ROOT + "test1.bam"
BAMFILE_B = ROOT + "test2.bam"


def test_r_bamcompare():
    bamifile1 = BAMFILE_A
    bamifile2 = BAMFILE_B
    ofile = ROOT + "r_bamcompare_output.bedgraph"
    ofiletype  = "bedgraph"
    norm  = "None"
    effective_genome_size = 0
    scalefactorsmethod = "None"
    operation = "ratio"
    pseudocount = 0
    # filtering options
    ignoreduplicates = False
    minmappingquality = 0
    samflaginclude = 0
    samflagexclude = 0
    minfraglen = 0
    maxfraglen = 0
    nproc = 1
    _ignorechr = []
    binsize = 500
    regions = []
    verbose = False

    # Call the Rust function
    r_bamcompare(
        bamifile1, bamifile2, ofile, ofiletype, norm, effective_genome_size, scalefactorsmethod, operation, 
        pseudocount, ignoreduplicates, minmappingquality, samflaginclude, samflagexclude, minfraglen, maxfraglen, nproc, 
        _ignorechr, binsize, regions, verbose
    )
    # Add assertions to verify the expected behavior
    expected = ['3R\t0\t500\t0.6213592\n',
            '3R\t500\t1000\t0.8252427\n',
            '3R\t1000\t1500\t0.2\n',]
    _foo = open(ofile, 'r')
    resp = _foo.readlines()
    _foo.close()
    assert f"{resp}" == f"{expected}", f"{resp} != {expected}"
    unlink(ofile)

def test_r_bamcompare_RPKM():
    bamifile1 = BAMFILE_A
    bamifile2 = BAMFILE_B
    ofile = ROOT + "r_bamcompare_output.bedgraph"
    ofiletype  = "bedgraph"
    norm  = "RPKM"
    effective_genome_size = 0
    scalefactorsmethod = "None"
    operation = "ratio"
    pseudocount = 0
    # filtering options
    ignoreduplicates = False
    minmappingquality = 0
    samflaginclude = 0
    samflagexclude = 0
    minfraglen = 0
    maxfraglen = 0
    nproc = 1
    _ignorechr = []
    binsize = 500
    regions = []
    verbose = False

    # Call the Rust function
    r_bamcompare(
        bamifile1, bamifile2, ofile, ofiletype, norm, effective_genome_size, scalefactorsmethod, operation, 
        pseudocount, ignoreduplicates, minmappingquality, samflaginclude, samflagexclude, minfraglen, maxfraglen, nproc, 
        _ignorechr, binsize, regions, verbose
    )
    # Add assertions to verify the expected behavior
    expected = ['3R\t0\t500\t0.6213592\n',
            '3R\t500\t1000\t0.8252427\n',
            '3R\t1000\t1500\t0.2\n',]
    _foo = open(ofile, 'r')
    resp = _foo.readlines()
    _foo.close()
    assert f"{resp}" == f"{expected}", f"{resp} != {expected}"
    unlink(ofile)