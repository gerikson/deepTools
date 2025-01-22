import pytest
from deeptools.alignmentSieve2 import r_alignmentsieve
import os.path


ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
BAMFILE_A = ROOT + "test_paired.bam"
BAMFILE_B = ROOT + "test_paired2.bam"

def test_r_alignmentsieve_basic():
    bamifile = BAMFILE_A
    ofile = BAMFILE_B
    nproc = 4
    filter_metrics = ROOT + "filter_metrics_basic.txt"
    filtered_out_readsfile = ROOT + "filtered_out_reads.bam"
    verbose = False
    shift = []
    _bed = False
    filter_rna_strand = "None"
    min_mapping_quality = 30
    sam_flag_incl = 0
    sam_flag_excl = 0
    _blacklist = ""
    min_fragment_length = 0
    max_fragment_length = 1000
    extend_reads = 0
    center_reads = False

    result = r_alignmentsieve(
        bamifile, ofile, nproc, filter_metrics, filtered_out_readsfile, verbose, shift, _bed,
        filter_rna_strand, min_mapping_quality, sam_flag_incl, sam_flag_excl, _blacklist,
        min_fragment_length, max_fragment_length, extend_reads, center_reads
    )

    expected_content = [
        "#bamFilterReads --filterMetrics",
        "#File\tReads\tRemaining Total\tInitial Reads",
        "test_paired.bam\t49\t49"
    ]

    with open(filter_metrics, 'r') as f:
        lines = f.readlines()
    
    # Remove the file path prefix from the third line
    lines[2] = os.path.basename(lines[2])

    # Strip newline characters from the lines
    lines = [line.strip() for line in lines]

    assert lines == expected_content


def test_r_alignmentsieve_with_shift():
    bamifile = BAMFILE_A
    ofile = BAMFILE_B
    nproc = 4
    filter_metrics = ROOT + "filter_metrics_shift.txt"
    filtered_out_readsfile = ROOT + "filtered_out_reads.bam"
    verbose = False
    shift = []
    _bed = False
    filter_rna_strand = "None"
    min_mapping_quality = 30
    sam_flag_incl = 0
    sam_flag_excl = 0
    _blacklist = ""
    min_fragment_length = 0
    max_fragment_length = 1000
    extend_reads = 0
    center_reads = False

    result = r_alignmentsieve(
        bamifile, ofile, nproc, filter_metrics, filtered_out_readsfile, verbose, shift, _bed,
        filter_rna_strand, min_mapping_quality, sam_flag_incl, sam_flag_excl, _blacklist,
        min_fragment_length, max_fragment_length, extend_reads, center_reads
    )
    
    expected_content = [
        "#bamFilterReads --filterMetrics",
        "#File\tReads\tRemaining Total\tInitial Reads",
        "test_paired.bam\t49\t49"
    ]

    with open(filter_metrics, 'r') as f:
        lines = f.readlines()
    
    # Remove the file path prefix from the third line
    lines[2] = os.path.basename(lines[2])

    # Strip newline characters from the lines
    lines = [line.strip() for line in lines]

    assert lines == expected_content


def test_r_alignmentsieve_with_filtering():
    bamifile = BAMFILE_A
    ofile = BAMFILE_B
    nproc = 4
    filter_metrics = ROOT + "filter_metrics_withFiltering.txt"
    filtered_out_readsfile = ROOT +"filtered_out_reads.bam"
    verbose = False
    shift = []
    _bed = False
    filter_rna_strand = "forward"
    min_mapping_quality = 30
    sam_flag_incl = 0
    sam_flag_excl = 0
    _blacklist = ""
    min_fragment_length = 100
    max_fragment_length = 500
    extend_reads = 0
    center_reads = False

    result = r_alignmentsieve(
        bamifile, ofile, nproc, filter_metrics, filtered_out_readsfile, verbose, shift, _bed,
        filter_rna_strand, min_mapping_quality, sam_flag_incl, sam_flag_excl, _blacklist,
        min_fragment_length, max_fragment_length, extend_reads, center_reads
    )
    
    expected_content = [
        "#bamFilterReads --filterMetrics",
        "#File\tReads\tRemaining Total\tInitial Reads",
        "test_paired.bam\t27\t49"
    ]

    with open(filter_metrics, 'r') as f:
        lines = f.readlines()
    
    # Remove the file path prefix from the third line
    lines[2] = os.path.basename(lines[2])

    # Strip newline characters from the lines
    lines = [line.strip() for line in lines]

    assert lines == expected_content
