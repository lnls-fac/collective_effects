"""."""

import gzip as _gzip
import os as _os
import pickle as _pickle
from collections import namedtuple as _namedtuple


def get_namedtuple(name, field_names, values=None):
    """Return an instance of a namedtuple Class.

    Inputs:
        - name:  Defines the name of the Class (str).
        - field_names:  Defines the field names of the Class (iterable).
        - values (optional): Defines field values . If not given, the value of
            each field will be its index in 'field_names' (iterable).

    Raises ValueError if at least one of the field names are invalid.
    Raises TypeError when len(values) != len(field_names)
    """
    if values is None:
        values = range(len(field_names))
    field_names = [f.replace(' ', '_') for f in field_names]
    return _namedtuple(name, field_names)(*values)


def save_pickle(data, fname, overwrite=False, makedirs=False, compress=False):
    """Save data to file in pickle format.

    Args:
        data (any builtin type): python object to be saved
        fname (str): name of the file to be saved. With or without ".pickle"."
        overwrite (bool, optional): whether to overwrite existing file.
            Defaults to False.
        makedirs (bool, optional): create dir, if it does not exist.
            Defaults to False.
        compress (bool, optional): If True, the file will be saved in
            compressed format, using gzip library. Defaults to False.

    Raises:
        FileExistsError: in case `overwrite` is `False` and file exists.

    """
    if not fname.endswith(('.pickle', '.pkl')):
        fname += '.pickle'

    if not overwrite and _os.path.isfile(fname):
        raise FileExistsError(f'file {fname} already exists.')

    if makedirs:
        dirname = _os.path.dirname(fname)
        if not _os.path.exists(dirname):
            _os.makedirs(dirname)

    func = _gzip.open if compress else open
    with func(fname, 'wb') as fil:
        _pickle.dump(data, fil)


def load_pickle(fname):
    """Load ".pickle" file.

    Args:
        fname (str): Name of the file to load. May or may not contain the
            ".pickle" extension.

    Returns:
        data (any builtin type): content of file as a python object.

    """
    if not fname.endswith(('.pickle', '.pkl')):
        fname += '.pickle'

    func = _gzip.open if is_gzip_file(fname) else open

    with func(fname, 'rb') as fil:
        data = _pickle.load(fil)
    return data


def is_gzip_file(fname):
    """Check if file is compressed with gzip.

    Args:
        fname (str): filename.

    Returns:
        bool: whether file is compressed with gzip.
    """
    # thanks to https://stackoverflow.com/questions/3703276/how-to-tell-if-a-file-is-gzip-compressed
    with open(fname, 'rb') as fil:
        return fil.read(2) == b'\x1f\x8b'
