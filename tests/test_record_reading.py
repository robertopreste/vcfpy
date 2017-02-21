# -*- coding: utf-8 -*-
"""Tests for reading records

The tests in this suite read records from a temporary VCF file.  The main
point under test here is that the Records read from a VCF file with the
Cython-htslib binding have a compatible API to the pure Python types.
"""

import os
import textwrap

from vcfpy import Reader, CYHTSLIB_ENABLED, Substitution

import pytest


@pytest.fixture
def vcf_no_id():
    """Return VCF that has one record with empty ID value"""
    return textwrap.dedent("""
        ##fileformat=VCFv4.3
        ##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003
        20\t14370\t.\tG\tA\t29\t.\t.\tGT\t0/0\t1/0\t1/1
        """).lstrip()


@pytest.fixture
def vcf_one_id():
    """Return VCF that has one record with one ID value"""
    return textwrap.dedent("""
        ##fileformat=VCFv4.3
        ##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003
        20\t14370\trs1\tG\tA\t29\t.\t.\tGT\t0/0\t1/0\t1/1
        """).lstrip()


@pytest.fixture
def vcf_two_ids():
    """Return VCF that has one record with two ID values"""
    return textwrap.dedent("""
        ##fileformat=VCFv4.3
        ##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003
        20\t14370\trs1;rs2\tG\tA\t29\t.\t.\tGT\t0/0\t1/0\t1/1
        """).lstrip()


@pytest.fixture
def vcf_no_filter():
    """Return VCF that has one record with empty FILTER value"""
    return textwrap.dedent("""
        ##fileformat=VCFv4.3
        ##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FILTER=<ID=PASS,Description="All filter passed">
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003
        20\t14370\t.\tG\tA\t29\t.\t.\tGT\t0/0\t1/0\t1/1
        """).lstrip()


@pytest.fixture
def vcf_pass_filter():
    """Return VCF that has one record with PASS FILTER value"""
    return textwrap.dedent("""
        ##fileformat=VCFv4.3
        ##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FILTER=<ID=PASS,Description="All filter passed">
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003
        20\t14370\t.\tG\tA\t29\tPASS\t.\tGT\t0/0\t1/0\t1/1
        """).lstrip()


@pytest.fixture
def vcf_two_filters():
    """Return VCF that has one record with two FILTER values"""
    return textwrap.dedent("""
        ##fileformat=VCFv4.3
        ##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FILTER=<ID=PASS,Description="All filter passed">
        ##FILTER=<ID=q10,Description="Quality below 10">
        ##FILTER=<ID=q20,Description="Quality below 20">
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003
        20\t14370\trs1;rs2\tG\tA\t29\tq10;q20\t.\tGT\t0/0\t1/0\t1/1
        """).lstrip()


@pytest.fixture
def vcf_all_types_info():
    """Return VCF that shows the possible types of values for the INFO field
    """
    return textwrap.dedent("""
        ##fileformat=VCFv4.3
        ##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
        ##FILTER=<ID=PASS,Description="All filter passed">
        ##FILTER=<ID=q10,Description="Quality below 10">
        ##FILTER=<ID=q20,Description="Quality below 20">
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003
        20\t14370\t.\tG\tA\t29\tPASS\t.\tGT\t0/0\t1/0\t1/1
        """).lstrip()


def create_file(vcf_txt, tmpdir):
    """Helper for creating VCF with text content in temporary directory
    """
    p = tmpdir.join('test.vcf')
    p.write(vcf_txt)
    return str(p)


def test_read_vcf_base_fields(vcf_no_id, tmpdir):
    """Test reading of basic fields from VCF"""
    path = create_file(vcf_no_id, tmpdir)
    reader = Reader.from_path(path)
    it = iter(reader)
    record = next(it)
    assert record.CHROM == '20'
    assert record.POS == 14370
    assert record.REF == 'G'
    assert len(record.ALT) == 1
    assert record.ALT[0] == Substitution('SNV', 'A')
    assert record.QUAL == 29
    with pytest.raises(StopIteration):  # at end
        next(it)


def test_read_vcf_no_id(vcf_no_id, tmpdir):
    """Test reading VCF that has a record with one ID"""
    path = create_file(vcf_no_id, tmpdir)
    reader = Reader.from_path(path)
    it = iter(reader)
    record = next(it)
    assert record.ID == []
    with pytest.raises(StopIteration):  # at end
        next(it)


def test_read_vcf_one_id(vcf_one_id, tmpdir):
    """Test reading VCF that has a record with one ID"""
    path = create_file(vcf_one_id, tmpdir)
    reader = Reader.from_path(path)
    it = iter(reader)
    record = next(it)
    assert record.ID == ['rs1']
    with pytest.raises(StopIteration):  # at end
        next(it)


def test_read_vcf_two_ids(vcf_two_ids, tmpdir):
    """Test reading VCF that has a record with two IDs"""
    path = create_file(vcf_two_ids, tmpdir)
    reader = Reader.from_path(path)
    it = iter(reader)
    record = next(it)
    assert record.ID == ['rs1', 'rs2']
    with pytest.raises(StopIteration):
        next(it)


def test_read_vcf_no_filter(vcf_no_filter, tmpdir):
    """Test reading VCF that has a record with no filter"""
    path = create_file(vcf_no_filter, tmpdir)
    reader = Reader.from_path(path)
    it = iter(reader)
    record = next(it)
    assert record.FILTER == []
    with pytest.raises(StopIteration):
        next(it)


def test_read_vcf_pass_filter(vcf_pass_filter, tmpdir):
    """Test reading VCF that has a record with PASS filter"""
    path = create_file(vcf_pass_filter, tmpdir)
    reader = Reader.from_path(path)
    it = iter(reader)
    record = next(it)
    assert record.FILTER == ['PASS']
    with pytest.raises(StopIteration):
        next(it)


def test_read_vcf_two_filters(vcf_two_filters, tmpdir):
    """Test reading VCF that has a record two filter entries"""
    path = create_file(vcf_two_filters, tmpdir)
    reader = Reader.from_path(path)
    it = iter(reader)
    record = next(it)
    assert record.FILTER == ['q10', 'q20']
    with pytest.raises(StopIteration):
        next(it)


def test_read_vcf_all_types_info(vcf_all_types_info, tmpdir):
    """Test reading VCF that shows all possible types of values for the INFO
    field
    """
    path = create_file(vcf_all_types_info, tmpdir)
    reader = Reader.from_path(path)
    assert reader is not None
