import time
import os

from parutils import pprint
from globals import alpha, comm

def sz(base64_input):
    """ Return size in bytes of base64 encoded input. 

    Ref: https://blog.aaronlenoir.com/2017/11/10/get-original-length-from-base-64-string/
    """

    b_padded = base64_input.split(str.encode(":"))[1]
    pad_size = b_padded.count(str.encode("="))
    b_len_without_pad = len(b_padded)-4
    byte_len = (b_len_without_pad *3)/4 +(3-pad_size)-1
    return byte_len

def statusstr(status):
    return "(ok)" if status else "(not ok)"

def expand_scientific_notation(numstr):
    was_neg = False
    if not ("e" in numstr):
        return numstr
    if numstr.startswith('-'):
        numstr = numstr[1:]
        was_neg = True 
    str_vals = str(numstr).split('e')
    coef = float(str_vals[0])
    exp = int(str_vals[1])
    return_val = ''
    if int(exp) > 0:
        return_val += str(coef).replace('.', '')
        return_val += ''.join(['0' for _ in range(0, abs(exp - len(str(coef).split('.')[1])))])
    elif int(exp) < 0:
        return_val += '0.'
        return_val += ''.join(['0' for _ in range(0, abs(exp) - 1)])
        return_val += str(coef).replace('.', '')
    if was_neg:
        return_val='-'+return_val
    return return_val

def fmt(num):
    num_truncated_precision = f'{num:.2}'
    return expand_scientific_notation(num_truncated_precision)

class timer:
    def __init__(self, name, verbosity=0):
        self.name = name
        self.verbosity = verbosity
        self.fill = "  " * self.verbosity 
    def __enter__(self):
        self.start = time.perf_counter()
        return self
    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.end = time.perf_counter()
        if int(os.environ['verbosity']) >= self.verbosity:
            pprint("%stime (%s): %s s" % (self.fill, self.name, fmt((self.end-self.start)) ))

class timerperEA:
    """ Timer class that times each authority's time. Simulated by dividing 
    the time taken by the wrapped function by alpha, the number of authorities. """
    def __init__(self, name, verbosity=0):
        self.name = name
        self.verbosity = verbosity
        self.fill = "  " * self.verbosity 
    def __enter__(self):
        self.start = time.perf_counter()
        return self
    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.end = time.perf_counter()
        if int(os.environ['verbosity']) >= self.verbosity:
            pprint("%stime (%s): %s s (per EA)" % (self.fill, self.name, fmt((self.end-self.start)/ alpha) ))


def timed(f):
    def timed_f(*args, **kwargs):
        with timer(f.__name__) as t:
            ret = f(*args, **kwargs)
            comm.barrier()
        return ret

    timed_f.__name__ = f.__name__
    return timed_f

def timedperEA(f):
    def timedperEA_f(*args, **kwargs):
        with timerperEA(f.__name__) as t:
            ret = f(*args, **kwargs)
            comm.barrier()
        return ret
    timedperEA_f.__name__ = f.__name__
    return timedperEA_f

def retval(f):
    def retval_f(*args, **kwargs):
        ret = f(*args, **kwargs)
        pprint("retval (%s): %s" % (f.__name__, ret))
        return ret

    retval_f.__name__ = f.__name__
    return retval_f
