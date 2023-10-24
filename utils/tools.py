from collections import defaultdict
import re

def defaultdict_factory():
    return defaultdict(defaultdict_factory)

def dictify(dd):
    dd = dict(dd)
    for k in dd.keys():
        if isinstance(dd[k], defaultdict):
            dd[k] = dictify(dd[k])
    return dd

def rgx_match(rgx, mystr):
    match = re.match(rgx, mystr)
    return match is not None and match.group(0) == mystr