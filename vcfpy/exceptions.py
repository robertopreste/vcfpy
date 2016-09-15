# -*- coding: utf-8 -*-
"""Exceptions for the vcfpy module
"""

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'


class VCFPyException(RuntimeError):
    """Base class for module's exception"""


class InvalidHeaderException(VCFPyException):
    """Raised in the case of invalid header formatting"""


class InvalidRecordException(VCFPyException):
    """Raised in the case of invalid record formatting"""


class IncorrectVCFFormat(VCFPyException):
    """Raised on problems parsing VCF"""