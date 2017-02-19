# -*- coding: utf-8 -*-
"""Tests for cyhtslib"""

import os

import pytest

from vcfpy import cyhtslib


def test_read_text():
    path = os.path.join(os.path.dirname(__file__), '../vcfs/full_vcf43.vcf')
    r = cyhtslib.VCFFile(path)
    # assert len(r.header.lines) == 18
    assert r.samples
    assert r.samples.names == ['NA00001', 'NA00002', 'NA00003']
    records = []
    for record in r:
        records.append(record)
    assert len(records) == 5


def test_read_bgzip():
    path = os.path.join(os.path.dirname(__file__), '../vcfs/full_vcf43.vcf.gz')
    r = cyhtslib.VCFFile(path)
    # assert len(r.header.lines) == 18
    assert r.samples
    assert r.samples.names == ['NA00001', 'NA00002', 'NA00003']
    records = []
    for record in r:
        records.append(record)
    assert len(records) == 5
