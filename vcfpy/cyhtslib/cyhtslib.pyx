"""Cython-based wrapper around htslib, enabling vcfpy to use faster I/O

Based on code from cyvcf2.
"""

from __future__ import print_function

# Python imports

import locale
import os
import sys

from vcfpy.header import SamplesInfos

# Cython imports

from libc cimport stdlib
from cython cimport view

#: Encoding to use for byte/string conversion
ENC = locale.getpreferredencoding()

# Overcome lack of __file__ in cython
import inspect
if not hasattr(sys.modules[__name__], '__file__'):
    __file__ = inspect.getfile(inspect.currentframe())


cdef to_bytes(s, enc=ENC):
    """Helper function for ensuring ``bytes`` type"""
    if not isinstance(s, bytes):
        return s.encode(enc)
    return s


cdef from_bytes(s):
    """Helper function for ensuring ``string`` type"""
    if isinstance(s, bytes):
        try:
            return s.decode(ENC)
        except UnicodeDecodeError:
            return s.decode('utf8')
    return s


cdef class VCFFile(object):
    """Representation of BCF/VCF file"""
    cdef htsFile *hts
    cdef bcf_hdr_t *hdr
    cdef tbx_t *idx
    cdef hts_idx_t *hidx
    cdef int n_samples
    cdef int PASS
    cdef bytes fname
    cdef bint lazy
    #: vcfpy.headre.SamplesInfo object
    cdef public object samples
    cdef list _seqnames
    # holds a lookup of format field -> type.
    cdef dict format_types

    def __init__(self, fname, mode='r', lazy=False, samples=None, threads=None):
        # Prepare file name, open file, balk out in the case of errors
        if fname == b'-' or fname == '-':
            fname = b'/dev/stdin'
        if not os.path.exists(fname):
            raise Exception("bad path: %s" % fname)
        fname, mode = to_bytes(fname), to_bytes(mode)
        self.hts = hts_open(fname, mode)
        if self.hts == NULL:
            raise IOError("Error opening %s" % fname)
        if self.hts.format.format != vcf and self.hts.format.format != bcf:
            raise IOError("%s if not valid bcf or vcf" % fname)

        # Read BCF header
        cdef bcf_hdr_t *hdr
        hdr = self.hdr = bcf_hdr_read(self.hts)
        # Set samples to be pulled out, if only limited to a sub set
        # TODO: pure python implementation is lacking this
        if samples is not None:
            self.set_samples(samples)
        self.n_samples = bcf_hdr_nsamples(self.hdr)
        # Extract SamplesInfos from header, possibly after limiting to subset
        if samples is None:
            samples = self.get_samples()
        self.samples = SamplesInfos(samples)
        # Initialize members
        self.PASS = -1
        self.fname = to_bytes(fname)
        self.lazy = lazy
        self._seqnames = []
        self.format_types = {}
        if threads is not None:
            self.set_threads(threads)

    def set_threads(self, int n):
        """Sets number of reader/writer threads in this object's htsfile"""
        v = hts_set_threads(self.hts, n)
        if v < 0:
            raise Exception("error setting number of threads: %d" % v)

    def set_samples(self, samples):
        """Set the samples to be pulled from the VCF; this must be called before any iteration.
        Parameters
        ----------
        samples: list
            list of samples to extract.
        """
        if samples is None:
            samples = "-".encode()
        if isinstance(samples, list):
            samples = to_bytes(",".join(samples))
        else:
            samples = to_bytes(samples)

        ret = bcf_hdr_set_samples(self.hdr, <const char *>samples, 0)
        assert ret >= 0, ("error setting samples", ret)
        if ret != 0 and samples != "-":
            s = samples.split(",")
            if ret < len(s):
                # TODO: update name, makes no sense here
                sys.stderr.write("warning: not all samples in PED found in VCF\n")

    def get_samples(self):
        """Return list of samples from file
        """
        return [from_bytes(self.hdr.samples[i])
                for i in range(0, bcf_hdr_nsamples(self.hdr))]

    def __iter__(self):
        return self

    def __next__(self):
        cdef bcf1_t *b = bcf_init()
        cdef int ret
        with nogil:
            ret = bcf_read(self.hts, self.hdr, b)
        if ret >= 0:
            return newRecord(b, self)
        else:
            bcf_destroy(b)
        raise StopIteration

    def __dealloc__(self):
        """Deallocation for VCFFile

        - deallocate header struct
        - close HTS file
        - free any index-related memory
        """
        if self.hdr != NULL:
            bcf_hdr_destroy(self.hdr)
            self.hdr = NULL
        if self.hts != NULL:
            hts_close(self.hts)
            self.hts = NULL
        if self.idx != NULL:
            tbx_destroy(self.idx)
        if self.hidx != NULL:
            hts_idx_destroy(self.hidx)


cdef class Record(object):
    #: Pointer to the C struct with the BCF record
    cdef bcf1_t *b
    #: Reference to the owning VCFFile
    cdef VCFFile vcf

    def __init__(self, *args, **kwargs):
        raise TypeError("Variant object cannot be instantiated directly.")

    def __cinit__(self):
        self.b = NULL

    def __repr__(self):
        return "Record(%s:%d %s/%s)" % (self.CHROM, self.POS, self.REF, ",".join(self.ALT))

    def __str__(self):
        cdef kstring_t s
        s.s, s.l, s.m = NULL, 0, 0
        vcf_format(self.vcf.hdr, self.b, &s)
        try:
            return s.s[:s.l].decode()
        finally:
            stdlib.free(s.s)

    def __dealloc__(self):
        """Perform deallocation

        - free associated BCF record
        - free all allocated buffers
        """
        if self.b is not NULL:
            bcf_destroy(self.b)
            self.b = NULL

    property CHROM:
        """Return string with the chromosome name"""
        def __get__(self):
            return bcf_hdr_id2name(self.vcf.hdr, self.b.rid).decode()

    property POS:
        """``int`` with 1-based start position of variant"""
        def __get__(self):
            return self.b.pos + 1

    property ID:
        """Return value of ID from the VCF field"""
        def __get__(self):
            cdef char *id = self.b.d.id
            if id == b'.':
                return []
            else:
                return id.decode().split(';')

    property REF:
        """Return ``str`` with reference allele"""
        def __get__(self):
            return self.b.d.allele[0].decode()

    property ALT:
        """Alternative alleles, list of ``str`` for now"""
        def __get__(self):
            cdef int i
            return [self.b.d.allele[i].decode() for i in range(1, self.b.n_allele)]

    property QUAL:
        """The quality value, can be ``None``"""
        def __get__(self):
            cdef float q = self.b.qual
            if bcf_float_is_missing(q):
                return None
            else:
                return q

    property FILTER:
        """Value of the FILTER field from VCF, as list of strings"""
        def __get__(self):
            cdef int i
            cdef int n = self.b.d.n_flt
            if n == 1:
                if self.vcf.PASS != -1:
                    if self.b.d.flt[0] == self.vcf.PASS:
                        return []
                else:
                    v = bcf_hdr_int2id(self.vcf.hdr, BCF_DT_ID, self.b.d.flt[0])
                    if v == b'PASS':
                        self.vcf.PASS = self.b.d.flt[0]
                        return ['PASS']
                    return v
            if n == 0:
                return []
            return b';'.join(bcf_hdr_int2id(self.vcf.hdr, BCF_DT_ID, self.b.d.flt[i]) for i in range(n))


cdef inline Record newRecord(bcf1_t *b, VCFFile vcf):
    """Construct new Record object in ``bcf1_t`` ``*b`` and with ``VCFFile`` ``vcf``
    """
    cdef Record rec = Record.__new__(Record)
    rec.b = b
    if not vcf.lazy:
        with nogil:
            bcf_unpack(rec.b, 15)
    else:
        with nogil:
            bcf_unpack(rec.b, 1|2|4)
    rec.vcf = vcf
    return rec

