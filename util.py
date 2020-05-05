
from collections import OrderedDict
import numpy as np
import re
import string




class ParamReaderMixin(object):
    """Mixin which adds parse(), read(), write() and asstring() to a dictionary"""

    def __copy__(self):
        return self.copy()

    def parse(self, s):
        # why could some key_values and key_quoted_values end in -e???????
        # key_quoted_value = re.compile(r'([A-Za-z_]+[A-Za-z0-9_]*)\s*=\s*(?:\'|["\{\}])([^"\{\}]+)(?:\'|["\{\}e+-])\s*')
        key_quoted_value = re.compile(r'([A-Za-z_]+[A-Za-z0-9_]*)\s*=\s*(?:\'|["\{\}])([^"\{\}]+)(?:\'|["\{\}+-])\s*')
        # key_value = re.compile(r'([A-Za-z_]+[A-Za-z0-9_]*)\s*=\s*([-0-9A-Za-z_.:\[\]()e+-]+)\s*')
        key_value = re.compile(r'([A-Za-z_]+[A-Za-z0-9_]*)\s*=\s*([-0-9A-Za-z_.:\[\]()+-]+)\s*')
        key_re = re.compile(r'([A-Za-z_]+[A-Za-z0-9_]*)\s*')

        s = s.strip()

        while 1:
            # Match quoted string first, then fall through to plain key=value
            m = key_quoted_value.match(s)
            if m is None:
                m = key_value.match(s) # removes matched key=value pair
                if m is not None:
                    s = key_value.sub('', s, 1)
                else:
                    # Just a key with no value
                    m = key_re.match(s)
                    if m is not None:
                        s = key_re.sub('', s, 1)
            else:
                s = key_quoted_value.sub('', s, 1)

            if m is None: break # No more matches

            key = m.group(1)
            try:
                value = m.group(2)
            except IndexError:
                # default value is None
                value = None

            # Try to convert to (list of) floats, ints
            try:
                numvalue = []
                for x in value.split():
                    if x.find('.') == -1:
                        numvalue.append(int(float(x)))
                    else:
                        numvalue.append(float(x))
                if len(numvalue) == 1:
                    numvalue = numvalue[0] # Only one number
                elif len(numvalue) == 3:
                    numvalue = np.array(numvalue) # 3-vector
                elif len(numvalue) == 9:
                    # 3x3 matrix, fortran ordering
                    numvalue = np.array(numvalue).reshape((3,3), order='F')
                else:
                    raise ValueError
                value = numvalue
            except (AttributeError, ValueError):
                pass

            # Parse boolean values, e.g 'T' -> True, 'F' -> False, 'T T F' -> [True, True, False]
            if isinstance(value, str):
                str_to_bool  = {'T':True, 'F':False}

                if len(value.split()) > 1:
                    if all([x in list(str_to_bool.keys()) for x in value.split() ]):
                        value = [str_to_bool[x] for x in value.split()]
                elif value in str_to_bool:
                    value = str_to_bool[value]

            self[key] = value

    def read(self, f):
        if isinstance(f, str):
            self.parse(f)
        elif hasattr(f, 'keys') and hasattr(f, '__getitem__'):
            self.update(f)
        else:
            try:
                for line in f:
                    self.parse(line)
            except TypeError:
                raise TypeError("Don't know how to read from object - "+str(f))


    def __repr__(self):
        return "%s('%s')" % (self.__class__.__name__, str(self))

    def __str__(self):
        return self.asstring()

    def asstring(self, sep=' '):
        if len(self) == 0: return ''

        type_val_map = {(bool,True): None,
                        (bool,False): 'F',
                        (np.bool_,True): None,
                        (np.bool_,False): 'F'}

        s = ''
        for key in list(self.keys()):
            val = self[key]

            if hasattr(val, '__iter__'):
                val = ' '.join( str(type_val_map.get((type(x),x), x))
                        for x in np.array(val).reshape(np.array(val).size, order='F'))
                val.replace('[', '')
                val.replace(']', '')
            else:
                val = type_val_map.get((type(val),val), val)

            if val is None:
                s = s + '%s%s' % (key, sep)
            elif type(val) == type('') and ' ' in val:
                s = s + '%s="%s"%s' % (key, val, sep)
            else:
                s = s + '%s=%s%s' % (key, str(val), sep)

        return s.strip()

    def write(self, f):
        f.write(self.asstring(sep='\n'))


class PuPyDictionary(OrderedDict, ParamReaderMixin):
    """Subclass of OrderedDict for reading key/value pairs from strings or files.
       The original order of items is maintained. Values that looks like floats or ints
       or lists of floats or ints are automatically converted on reading."""

    def __init__(self, source=None):
        OrderedDict.__init__(self)
        if source is not None:
            self.read(source)

    def __repr__(self):
        return ParamReaderMixin.__repr__(self)

    def __str__(self):
        return ParamReaderMixin.__str__(self)

    def copy(self):
        return PuPyDictionary(OrderedDict.copy(self))




def parse_params(s):
   """Read key=value pairs from a string or list of string and return a standard Python dictionary"""
   p = PuPyDictionary(s)
   return dict(p)
