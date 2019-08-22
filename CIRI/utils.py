import itertools
from datetime import datetime


def log(text):
    time_stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print('{} - {}'.format(time_stamp, text))


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
