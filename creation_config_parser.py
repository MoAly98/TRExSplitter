from schema import Schema, And, Use, Optional
import json


def keys_to_lower(mydict):
    newdict = {}
    for k in mydict.keys():
        if isinstance(mydict[k], dict):
            newdict[k.lower()] = keys_to_lower(mydict[k])
        else:
            newdict[k.lower()] = mydict[k]
    return newdict

def validate(cfg):
    mydict = json.load(cfg)
    mydict = keys_to_lower(mydict)
    the_schema = Schema(
                        # use str since json returns unicode
                        { 'systematics':  [str],
                            'samples':  [str],
                            'regions': [str],
                            Optional('exclude', default = ['']):  [str],
                            'general':
                                {
                                    'trexconfig':   str,
                                    'outdir':       str,
                                    Optional('othersysts', default = True):   bool,
                                    Optional('othersamples', default = True): bool,
                                    Optional('splitbysys', default = None):   [str],
                                    Optional('notsplitbysys', default = None):  [str],
                                    Optional('maxjobsperchunk', default = None):    int,
                                    Optional('splitby', default = None):    str,
                                    Optional('removefromsuff', default = [None]): [str],
                                    Optional('onlykeep', default = None): [str],
                                    Optional('bootstrap', default = 0): int,
                                    Optional('bootstrapnominalconfig', default = None): str,
                                    Optional('bootstrapparam', default = None): str,
                                    Optional('injectfromcfg', default = None): str,
                                    Optional('likelihoodx', default = [0,0]): [int,int],
                                    Optional('likelihoody', default = [0,0]): [int,int],
                                    Optional('usecache', default = False): bool,
                                }
                            })

    validated =  the_schema.validate(mydict)
    return validated

def process(cfgp):
    with open(cfgp, "r") as f:  validated = validate(f)

    if validated['general']['bootstrap'] != 0 and (validated['general']['bootstrapnominalconfig'] is None or validated['general']['bootstrapparam'] is None):
        raise ValueError("You must specify a bootstrapnominalconfig and bootstrapped param if you want to bootstrap")

    return validated
