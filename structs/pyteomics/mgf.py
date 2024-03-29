"""
mgf - read and write MS/MS data in Mascot Generic Format
========================================================

Summary
-------

`MGF <http://www.matrixscience.com/help/data_file_help.html>`_ is a simple
human-readable format for MS/MS data. It allows storing MS/MS peak lists and
exprimental parameters.

This module provides minimalistic infrastructure for access to data stored in
MGF files. The most important function is :py:func:`read`, which
reads spectra and related information as saves them into human-readable
:py:class:`dicts`.
Also, common parameters can be read from MGF file header with
:py:func:`read_header` function. :py:func:`write` allows creation of MGF
files.

Functions
---------

  :py:func:`read` - iterate through spectra in MGF file. Data from a
  single spectrum are converted to a human-readable dict.

  :py:func:`chain` - read multiple files at once.

  :py:func:`chain.from_iterable` - read multiple files at once, using an
  iterable of files.

  :py:func:`read_header` - get a dict with common parameters for all spectra
  from the beginning of MGF file.

  :py:func:`write` - write an MGF file.


Dependencies
------------

This module requires :py:mod:`numpy`.

-------------------------------------------------------------------------------
"""

#   Copyright 2012 Anton Goloborodko, Lev Levitsky
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

from . import auxiliary as aux
import numpy as np
import itertools as it

_comments = '#;!/'

@aux._file_reader()
def read(source=None, use_header=True):
    """Read an MGF file and return entries iteratively.

    Read the specified MGF file, **yield** spectra one by one.
    Each 'spectrum' is a :py:class:`dict` with four keys: 'm/z array',
    'intensity array', 'charge array' and 'params'. 'm/z array' and
    'intensity array' store :py:class:`numpy.ndarray`'s of floats,
    'charge array' is a masked array (:py:class:`numpy.ma.MaskedArray`) of ints,
    and 'params' stores a :py:class:`dict` of parameters (keys and values are
    :py:class:`str`, keys corresponding to MGF, lowercased).

    Parameters
    ----------

    source : str or file or None, optional
        A file object (or file name) with data in MGF format. Default is
        :py:const:`None`, which means read standard input.

    use_header : bool, optional
        Add the info from file header to each dict. Spectrum-specific parameters
        override those from the header in case of conflict.
        Default is :py:const:`True`.

    Returns
    -------

    out : FileReader
    """
    header = read_header(source)
    reading_spectrum = False
    params = {}
    masses = []
    intensities = []
    charges = []
    if use_header: params.update(header)
    for line in source:
        if not reading_spectrum:
            if line.strip() == 'BEGIN IONS':
                reading_spectrum = True
            # otherwise we are not interested; do nothing, just move along
        else:
            if not line.strip() or any(
                line.startswith(c) for c in _comments):
                    pass
            elif line.strip() == 'END IONS':
                reading_spectrum = False
                if 'pepmass' in params:
                    try:
                        pepmass = tuple(map(float, params['pepmass'].split()))
                    except ValueError:
                        raise aux.PyteomicsError('MGF format error: cannot parse '
                                'PEPMASS = {}'.format(params['pepmass']))
                    else:
                        params['pepmass'] = pepmass + (None,)*(2-len(pepmass))
                if isinstance(params.get('charge'), str):
                    params['charge'] = aux._parse_charge(params['charge'], True)
                out = {'params': params,
                       'm/z array': np.array(masses),
                       'intensity array': np.array(intensities),
                       'charge array': np.ma.masked_equal(charges, 0)}
                yield out
                del out
                params = dict(header) if use_header else {}
                masses = []
                intensities = []
                charges = []
            else:
                l = line.split('=', 1)
                if len(l) > 1: # spectrum-specific parameters!
                    params[l[0].lower()] = l[1].strip()
                elif len(l) == 1: # this must be a peak list
                    l = line.split()
                    if len(l) >= 2:
                        try:
                            masses.append(float(l[0]))            # this may cause
                            intensities.append(float(l[1]))       # exceptions...
                            charges.append(aux._parse_charge(l[2]) if len(l) > 2 else 0)
                        except ValueError:
                            raise aux.PyteomicsError(
                                 'Error when parsing %s. Line:\n%s' %
                                 (source, line))

@aux._keepstate
def read_header(source):
    """
    Read the specified MGF file, get search parameters specified in the header
    as a :py:class:`dict`, the keys corresponding to MGF format (lowercased).

    Parameters
    ----------

    source : str or file
        File name or file object representing an file in MGF format.

    Returns
    -------

    header : dict
    """
    with aux._file_obj(source, 'r') as source:
        header = {}
        for line in source:
            if line.strip() == 'BEGIN IONS':
                break
            l = line.split('=')
            if len(l) == 2:
                key = l[0].lower()
                val = l[1].strip()
                header[key] = val
        if 'charge' in header:
            header['charge'] = aux._parse_charge(header['charge'], True)
        return header

_key_order = ['title', 'pepmass', 'rtinseconds', 'charge']

def _pepmass_repr(k, pepmass):
    outstr = k.upper() + '='
    if not isinstance(pepmass, (str, int, float)): # assume iterable
        try:
            outstr += ' '.join(str(x) for x in pepmass if x is not None)
        except TypeError:
            raise aux.PyteomicsError(
                    'Cannot handle parameter: PEPMASS = {}'.format(pepmass))
    else:
        outstr += str(pepmass)
    return outstr

def _charge_repr(k, charge):
    return '{}={}'.format(k.upper(), aux._parse_charge(str(charge)))

def _default_repr(key, val):
    return '{}={}'.format(key.upper(), val)

_value_converters = {'pepmass': _pepmass_repr, 'charge': _charge_repr}

def _key_value_line(key, val):
    return _value_converters.get(key, _default_repr)(key, val) + '\n'

def write(spectra, output=None, header=''):
    """
    Create a file in MGF format.

    Parameters
    ----------

    spectra : iterable
        A sequence of dictionaries with keys 'm/z array', 'intensity array',
        and 'params'. 'm/z array' and 'intensity array' should be sequences of
        :py:class:`int`, :py:class:`float`, or :py:class:`str`. Strings will
        be written 'as is'. The sequences should be of equal length, otherwise
        excessive values will be ignored.

        'params' should be a :py:class:`dict` with keys corresponding to MGF
        format. Keys must be strings, they will be uppercased and used as is,
        without any format consistency tests. Values can be of any type allowing
        string representation.

        .. note:: Key order is defined by :py:data:`_key_order` variable. Value
            representation is configured via :py:data:`_value_converters`.

        'charge array' can also be specified.

    output : str or file or None, optional
        Path or a file-like object open for writing. If an existing file is
        specified by file name, it will be opened for appending. In this case
        writing with a header can result in violation of format conventions.
        Default value is :py:const:`None`, which means using standard output.

    header : dict or (multiline) str or list of str, optional
        In case of a single string or a list of strings, the header will be
        written 'as is'. In case of dict, the keys (must be strings) will be
        uppercased.

    Returns
    -------

    output : file
    """
    with aux._file_obj(output, 'a') as output:

        if isinstance(header, dict):
            head_dict = header.copy()
            head_lines = [_key_value_line(k, v) for k, v in header.items()]
            head_str = '\n'.join(head_lines)
        else:
            if isinstance(header, str):
                head_str = header
                head_lines = header.split('\n')
            else:
                head_lines = list(header)
                head_str = '\n'.join(header)
            head_dict = {}
            for line in head_lines:
                if not line.strip() or any(
                    line.startswith(c) for c in _comments):
                   continue
                l = line.split('=')
                if len(l) == 2:
                    head_dict[l[0].lower()] = l[1].strip()
        if head_str:
            output.write(head_str + '\n\n')

        for spectrum in spectra:
            output.write('BEGIN IONS\n')
            found = set()
            for key in it.chain(_key_order, spectrum['params']):
                if key not in found and key in spectrum['params']:
                    found.add(key)
                    val = spectrum['params'][key]
                    if val != head_dict.get(key):
                        output.write(_key_value_line(key, val))

            try:
                for m, i, c in zip(spectrum['m/z array'],
                        spectrum['intensity array'],
                        spectrum.get('charge array', it.cycle((None,)))):
                    output.write('{} {} {}\n'.format(
                        m, i,
                        (c if c not in (None, np.nan, np.ma.masked) else '')))
            except KeyError:
                raise aux.PyteomicsError("'m/z array' and 'intensity array'"
                        " must be present in all spectra.")
            output.write('END IONS\n\n')
        return output

chain = aux._make_chain(read, 'read')
