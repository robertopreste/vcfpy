# -*- coding: utf-8 -*-

# TODO: use Python dict in case of dict, move this to "compat" module
try:
    from cyordereddict import OrderedDict
except ImportError:
    from collections import OrderedDict

#: Boolean that gives whether cytslib is working (available and enabled)
CYHTSLIB_ENABLED = True
try:
    from .cyhtslib import ENABLED as _CYHTSLIB_ENABLED
    CYHTSLIB_ENABLED = _CYHTSLIB_ENABLED
except ImportError as e:
    CYHTSLIB_ENABLED = False


from .exceptions import VCFPyException, InvalidHeaderException, \
    InvalidRecordException, IncorrectVCFFormat, HeaderNotFound

from .header import Header, HeaderLine, SimpleHeaderLine, \
    AltAlleleHeaderLine, ContigHeaderLine, FilterHeaderLine, MetaHeaderLine, \
    PedigreeHeaderLine, SampleHeaderLine, CompoundHeaderLine, InfoHeaderLine, \
    FormatHeaderLine, SamplesInfos, FieldInfo, header_without_lines

from .record import Record, Call, AltRecord, Substitution, BreakEnd, \
    SingleBreakEnd, SymbolicAllele

from .record import SNV, MNV, DEL, INS, INDEL, SV, BND, SYMBOLIC, MIXED
from .record import HOM_REF, HET, HOM_ALT
from .record import FIVE_PRIME, THREE_PRIME, FORWARD, REVERSE

from .reader import Reader

from .writer import Writer

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
