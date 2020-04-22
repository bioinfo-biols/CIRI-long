import os
import sys
import re
import itertools
import threading
import _thread as thread


def quit_function(fn_name):
    sys.stderr.write('{0} took too long\n'.format(fn_name))
    sys.stderr.flush()
    thread.interrupt_main()


def exit_after(s):
    """
    use as decorator to exit process if
    function takes longer than s seconds
    """
    def outer(fn):
        def inner(*args, **kwargs):
            timer = threading.Timer(s, quit_function, args=[fn.__name__])
            timer.start()
            try:
                result = fn(*args, **kwargs)
            finally:
                timer.cancel()
            return result
        return inner
    return outer


def check_file(file_name):
    if os.path.exists(file_name) and os.path.isfile(file_name):
        return os.path.abspath(file_name)
    else:
        sys.exit('File: {}, not found'.format(file_name))


def check_dir(dir_name):
    if os.path.exists(dir_name):
        if os.path.isdir(dir_name):
            # Output directory already exists
            pass
        else:
            sys.exit('Directory: {}, clashed with existed files'.format(dir_name))
    else:
        os.mkdir(dir_name)
    return os.path.abspath(dir_name)


def to_str(bytes_or_str):
    # Return Instance of str
    if isinstance(bytes_or_str, bytes):
        value = bytes_or_str.decode('utf-8')
    else:
        value = bytes_or_str
    return value


def to_bytes(bytes_or_str):
    # Return Instance of bytes
    if isinstance(bytes_or_str, str):
        value = bytes_or_str.encode('utf-8')
    else:
        value = bytes_or_str
    return value


def empty_iter(iterable):
    try:
        first = next(iterable)
    except StopIteration:
        return None
    return itertools.chain([first], iterable)


def grouper(iterable, n, fillvalue=None):
    from itertools import zip_longest
    """
    Collect data info fixed-length chunks or blocks
    grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    """

    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=None)


def pairwise(iterable):
    from itertools import tee
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def tree():
    from collections import defaultdict
    return defaultdict(tree)


def flatten(x):
    """
    Flatten list of lists
    """
    import itertools

    flatted_list = list(itertools.chain(*x))
    return flatted_list
