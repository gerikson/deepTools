import deeptools.alignmentSieve as sieve
import os.path
from os import unlink
import hashlib
import pysam



ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
BAMFILE_FILTER = ROOT + "test_filtering.bam"
BEDFILE_FILTER = ROOT + "test_filtering.blacklist.bed"
PAIREDBAMFILE_FILTER = ROOT + "test_paired.bam"
BAMFILE_FILTER25 = ROOT + "test_sieve25.bam"

def test_sieve_ignoreDups():
    """
    Test how the ignore duplication works
    """    
    outfile = '/tmp/test_sieve25.bam'
    outfiltered = '/tmp/test_sieveFiltered.bam'
    outlog = '/tmp/test_sieve.log'
    args = '-b {} --smartLabels --ignoreDuplicates  -o {} --filterMetrics {} --filteredOutReads {}'.format(BAMFILE_FILTER25, outfile, outlog, outfiltered).split()
    sieve.main(args)

    _foo = open(outlog, 'r')
    resp = _foo.readlines()
    _foo.close()

    expected = ['#bamFilterReads --filterMetrics\n',
                '#File\tReads Remaining\tTotal Initial Reads\n',
                'test_sieve25\t238\t268\n']
    assert f"{resp}" == f"{expected}", f"{resp} != {expected}"
    unlink(outlog)
    h = hashlib.md5(pysam.view(outfile).encode('utf-8')).hexdigest()
    expectedh = 'ef17f211ed2a72007f7d2665852c3937'
    assert f'{h}' == f'{expectedh}'
    unlink(outfile)

    h1 = hashlib.md5(pysam.view(outfiltered).encode('utf-8')).hexdigest()
    expectedh = 'd41d8cd98f00b204e9800998ecf8427e'
    assert f"{h1}" == f"{expectedh}"
    unlink(outfiltered)



def test_sieve():
    """
    Test filtering a BAM file by MAPQ, flag, and blacklist
    """
    outfile = '/tmp/test_sieve.bam'
    outfiltered = '/tmp/test_sieveFiltered.bam'
    outlog = '/tmp/test_sieve.log'
    args = '-b {} --smartLabels --minMappingQuality 10 --samFlagExclude 512 -bl {} -o {} --filterMetrics {} --filteredOutReads {}'.format(BAMFILE_FILTER, BEDFILE_FILTER, outfile, outlog, outfiltered).split()
    sieve.main(args)

    _foo = open(outlog, 'r')
    resp = _foo.readlines()
    _foo.close()

    expected = ['#bamFilterReads --filterMetrics\n',
                '#File\tReads Remaining\tTotal Initial Reads\n',
                'test_filtering\t5\t193\n']
    assert f"{resp}" == f"{expected}", f"{resp} != {expected}"
    unlink(outlog)
    h = hashlib.md5(pysam.view(outfile).encode('utf-8')).hexdigest()
    expectedh = 'acbc4443fb0387bfd6c412af9d4fc414'
    assert f'{h}' == f'{expectedh}'
    unlink(outfile)

    h1 = hashlib.md5(pysam.view(outfiltered).encode('utf-8')).hexdigest()
    expectedh = 'b90befdd5f073f14acb9a38661f301ad'
    assert f"{h1}" == f"{expectedh}"
    unlink(outfiltered)


def test_sieve_BED():
    """
    Test alignmentSieve with the --BED option
    """
    outfile = '/tmp/test_sieve.bed'
    args = '-b {} --minMappingQuality 10 --BED -o {}'.format(PAIREDBAMFILE_FILTER, outfile).split()
    sieve.main(args)

    _foo = open(outfile, 'r')
    resp = _foo.readlines()
    _foo.close()

    expected = ['chr2\t5000026\t5000390\n',
                'chr2\t5000303\t5000711\n',
                'chr2\t5000384\t5000531\n',
                'chr2\t5000384\t5000531\n',
                'chr2\t5000559\t5000941\n',
                'chr2\t5000736\t5001171\n',
                'chr2\t5000819\t5001228\n',
                'chr2\t5000821\t5001158\n',
                'chr2\t5000821\t5001158\n',
                'chr2\t5000821\t5001158\n',
                'chr2\t5000834\t5001249\n',
                'chr2\t5000855\t5001277\n',
                'chr2\t5000867\t5001218\n',
                'chr2\t5000925\t5001023\n',
                'chr2\t5000925\t5001023\n',
                'chr2\t5000937\t5001338\n',
                'chr2\t5001010\t5001176\n',
                'chr2\t5001025\t5001431\n',
                'chr2\t5001050\t5001436\n',
                'chr2\t5001114\t5001413\n',
                'chr2\t5001115\t5001269\n',
                'chr2\t5001115\t5001269\n',
                'chr2\t5001226\t5001603\n',
                'chr2\t5001491\t5001527\n',
                'chr2\t5001700\t5001736\n']

    assert f"{resp}" == f"{expected}", f"{resp} != {expected}"
    unlink(outfile)


def test_sieve_BED_shift():
    """
    Test alignmentSieve --BED --shift
    """
    outfile = '/tmp/test_sieve_shift.bed'
    args = '-b {} --minMappingQuality 10 --BED -o {} --shift 1 -2 3 -4'.format(PAIREDBAMFILE_FILTER, outfile).split()
    sieve.main(args)

    _foo = open(outfile, 'r')
    resp = _foo.readlines()
    _foo.close()

    expected = ['chr2\t5000027\t5000388\n',
                'chr2\t5000307\t5000708\n',
                'chr2\t5000388\t5000528\n',
                'chr2\t5000385\t5000529\n',
                'chr2\t5000560\t5000939\n',
                'chr2\t5000737\t5001169\n',
                'chr2\t5000823\t5001225\n',
                'chr2\t5000825\t5001155\n',
                'chr2\t5000825\t5001155\n',
                'chr2\t5000825\t5001155\n',
                'chr2\t5000835\t5001247\n',
                'chr2\t5000859\t5001274\n',
                'chr2\t5000868\t5001216\n',
                'chr2\t5000929\t5001020\n',
                'chr2\t5000929\t5001020\n',
                'chr2\t5000941\t5001335\n',
                'chr2\t5001011\t5001174\n',
                'chr2\t5001026\t5001429\n',
                'chr2\t5001054\t5001433\n',
                'chr2\t5001118\t5001410\n',
                'chr2\t5001119\t5001266\n',
                'chr2\t5001119\t5001266\n',
                'chr2\t5001230\t5001600\n']

    assert f"{resp}" == f"{expected}", f"{resp} != {expected}"
    unlink(outfile)

