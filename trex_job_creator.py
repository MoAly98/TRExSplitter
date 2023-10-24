'''
This script is responsible to take in a configuration file from
trex-fitter and a job configuration and creates trex-fitter n-step
arguments that are split in a way specified in job config
'''
# System Imports
import sys, os
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
# Python Imports
from argparse import ArgumentParser
import re
from typing import Dict, List, Union
import numpy as np
# thqtools imports
from utils.logger import ColoredLogger
from creation_config_parser import process as parse_jobs_config
from utils.config_reader import process_config as parse_trex_config
from utils.tools import rgx_match

logger = ColoredLogger()

def get_args():
    '''
     Method to retreive arguments from command line
     Returns:
        An `ArgumentParser` object with the parsed values
    '''
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('jobs_cfg',     )
    parser.add_argument('--step',       default = 'n')
    parser.add_argument('-s', '--suff', default = '')

    return parser.parse_args()

class Bootstrapper(object):
    def __init__(self, n, nom_hfolder, nom_jobname, param, step):
        self.n = n
        self.nom_hfolder = nom_hfolder
        self.nom_jobname = nom_jobname
        self.param = param
        self.step = step

def main():
    '''
    This is the main function in the script, steering the workflow
    '''
    # ====== Retreive various arguments
    args = get_args()
    step = args.step

    # ====== Parse the jobs config
    options = parse_jobs_config(args.jobs_cfg)
    cfg = options["general"]["trexconfig"]
    outdir = options["general"]["outdir"]
    os.makedirs(outdir, exist_ok = True)



    # ====== Parse the trex-fitter config
    retrieve = {'Sample':  {'attrs': ['UseSystematics','Type']},
                'Region':   {},
                'Systematic': {'attrs': ['Exclude','Samples'], 'split_by_sc': True},
                'OutputDir':    {},
                'InputFolder': {},
                'NormFactor': {},
            }
    ( jobname, parsed_names, parsed_attrs) = parse_trex_config(cfg, **retrieve)

    cfg_samples, cfg_samp_attrs      = parsed_names['Sample'], parsed_attrs['Sample']
    cfg_regions, cfg_regions_attrs   = parsed_names['Region'], parsed_attrs['Region']
    cfg_systematics, cfg_systs_attrs = parsed_names['Systematic'], parsed_attrs['Systematic']
    cfg_normfactors = None
    if step == 'r':
        cfg_normfactors                  = parsed_names['NormFactor']

    histfolder =  next(iter(parsed_names['InputFolder']), '').strip('"')
    if histfolder == '':    histfolder = next(iter(parsed_names['OutputDir']), '').strip('"')
    if histfolder == '':    logger.warning("No histogram folder found in trex-fitter config")

    # ====== Remove Ghost samples from list of samples
    cfg_samples = [s for i, s in enumerate(cfg_samples) if cfg_samp_attrs[i]['Type']!='GHOST']
    cfg_samp_attrs =  [s for i, s in enumerate(cfg_samp_attrs) if cfg_samp_attrs[i]['Type']!='GHOST']
    # ====== Identify the data sample from the list of samples
    data =  [s for i, s in enumerate(cfg_samples) if cfg_samp_attrs[i]['Type']=='DATA']

    # ====== Prepare for bootstrapping
    bootstrap = Bootstrapper(0, None, None, None, None)
    n_bootstraps = options["general"]["bootstrap"]
    if n_bootstraps > 0:
        retrieve = {'InputFolder':    {},}
        ( nom_jobname, parsed_names, parsed_attrs) = parse_trex_config(options["general"]["bootstrapnominalconfig"], **retrieve)
        nom_histfolder =  next(iter(parsed_names['InputFolder']), '').strip('"')
        if nom_histfolder == '':    logger.error("No histogram folder found in BootstrapNominalConfig")
        if not os.path.exists(nom_histfolder): logger.error(f"BootstrapNominalConfig histogram folder {nom_histfolder} does not exist")

        bootstrap = Bootstrapper(n_bootstraps, nom_histfolder, nom_jobname, options["general"]["bootstrapparam"], step)

    # ====== Prepare for injection of fit results
    injectedfitcfg = options["general"]["injectfromcfg"]
    if injectedfitcfg is not None:
        retrieve = {'OutputDir':    {},}
        ( injectcfg_jobname, parsed_names, parsed_attrs) = parse_trex_config(injectedfitcfg, **retrieve)
        injectcfg_outdir =  next(iter(parsed_names['OutputDir']), '').strip('"')
        if injectcfg_outdir == '':    logger.error("No output directory found in InjectFromCfg")
        if not os.path.exists(injectcfg_outdir): logger.error(f"InjectFromCfg output directory {injectcfg_outdir} does not exist")
        injectedfitdir = f'{injectcfg_outdir}/{injectcfg_jobname}/Fits/'
        injectedfitfile = f'{injectedfitdir}/{injectcfg_jobname}.txt'
    else:
        injectedfitfile = None

    # ====== Make the list of jobs to run
    jobs, exp_files, regions = make_jobs(cfg_samples,cfg_regions, cfg_systematics, data, jobname, options, injectedfitfile, bootstrap, cfg_normfactors, run_ranking = step == 'r')

    numjobs = len(jobs) #n_bootstraps*len(jobs) if n_bootstraps > 0 else len(jobs)
    logger.info(f"You have generated {numjobs} jobs")

    if len(jobs) == 0:  logger.error("No jobs to run. Exiting...")

    # Prepare an hupdate file for merging by region
    if args.step in ['n','b']:
        if regions != '**': hupdate_file_maker(regions, jobname, options['general']['outdir'],histfolder, bootstrap)
        if regions == '**': hupdate_file_maker(cfg_regions, jobname, options['general']['outdir'],histfolder, bootstrap)

    # Chunk up the jobs and save them into txt files
    stepsize = len(jobs) if options["general"]["maxjobsperchunk"] is None else options["general"]["maxjobsperchunk"]

    expected_histogram_files = open(f"{outdir}histogram_files.txt", "w")
    job_ctr = 0
    for i, chunk in enumerate(split(jobs, stepsize)):
            jobs_file = open(f"{outdir}arguments_{i}{args.suff}.txt", "w")
            logger.info(f"Jobs will now be saved in file: {outdir}arguments_{i}{args.suff}.txt")
            for job in chunk:
                jobs_file.write(f'trex-fitter {step} {cfg} {job}')
                if step == 'n':
                    for f in exp_files[job_ctr]:    expected_histogram_files.write(f+'\n')
                    job_ctr += 1
            jobs_file.close()

    expected_histogram_files.close()


def split(jobs: List[str], stepsize: int):
    '''
    Method to generate a range of jobs based on a split size
    Args:
        jobs (list):   list of jobs that need to be split
        stepsize (int): Maximum number of jobs in a chunk
    Returns:
        Generator with the chunks of jobs according to stepsize
    '''
    for i in range(0, len(jobs), stepsize):
        yield jobs[i:i + stepsize]

def make_jobs(
    cfg_samples:  List[str],
    cfg_regions:  List[str],
    cfg_systs:  List[str],
    data: List[str],
    jobname: str,
    options:  Dict[str, Union[str, Dict[str, str]]],
    injectedfitfile: str,
    bootstrap: Bootstrapper,
    cfg_normfactors: List[str],
    run_ranking: bool = False,
):
    '''
    Method responsible for constructing the correct trex-fitter command
    that will be saved to be ran by user.
    Args:
        cfg_samples (List[str]):    List of non-ghost samples identified from trex-fitter config
        cfg_regions (List[str]):    List of regions identified from trex-fitter config
        cfg_systs (List[str]):      List of systematics identified from trex-fitter config
        data (List[str]):           The names of the data samples
        options (Dict[str, Union[str, Dict[str, str]]]):  The dictionary of settings parsed from the jobs config
        injectedfitfile (str):       The path to the injected fit results
        bootstrap (Bootstrapper):   The object containing the bootstrap settings

    Returns:
        List of trex-fitter job settings (Regions=x:Samples=y:Suffix=z)
    '''
    # ===============================================================
    # =================== Initialization ============================
    # ===============================================================
    # Declare a jobs list
    jobs = []
    all_exp_files = []

    # =================== Sugar ======================================
    # Should a particular part of the job-combo be removed from Suffix
    removefromsuff = options["general"]["removefromsuff"]
    assert all(suff_to_remove in [None, "SAMPLE", "REGION","SYSTEMATIC"] for suff_to_remove in removefromsuff), logger.error('removefromsuff option can be any/all of ["SAMPLE", "REGION","SYSTEMATIC"]')

    # =================== Object Retrieval ============================
    # Grab the regions, samples, systematics and excluded blocks
    # from the jobs config provided by user
    regions = options["regions"]
    samples = options["samples"]
    systematics = options["systematics"]
    exclude = options["exclude"]

    # Make sure the user can only use * in their config,
    # not .* (this is reserved for the code to use regex)
    if any(".*" in sample for sample in samples):
        logger.error("Cannot use .* in your sample names")
    if any(".*" in region for region in regions):
        logger.error("Cannot use .* in your region names")
    if any(".*" in excl for excl in exclude):
        logger.error("Cannot use .* in your excludes names")

    # Special cases if user wants to split jobs by
    # all samples/region/systematic:
    # i.e. +1 jobs per config sample/region/systematic
    if samples ==     ["*"]:    samples     = cfg_samples
    if regions ==     ["*"]:    regions     = cfg_regions
    if systematics == ["*"]:    systematics = cfg_systs

    # =================== Exclusions ============================
    # make regex-ready list of excludes
    exclude = [opt.replace("*",".*") for opt in options["exclude"]]
    # make commaa separted list of excluded things, to pass to job line
    exclude_str = ",".join(exclude).replace(".*","*")

    # Remove any regions/samples/systematics that need to be excluded
    x_excluded = lambda x: any(re.search('^'+excluded+'$', x) is not None for excluded in exclude)
    regions =     [r for r in regions     if not x_excluded(r)]
    samples =     [s for s in samples     if not x_excluded(s)]
    systematics = [s for s in systematics if not x_excluded(s)]


    # Make a list that excludes all samples/systematics
    # to be used for othersamples/othersysts
    excl_samp = ','.join(samples).replace(".*","*")
    not_none_systs = [sys for sys in systematics if sys != 'NONE']
    excl_syst = ','.join(not_none_systs).replace(".*","*")

    # =========================================================================
    # =================== Collect remains, or not? ============================
    # =========================================================================
    # Should a job be added that combines all systematics not specified in config
    othersysts  = options["general"]["othersysts"]
    # Should a job be added that combines all samples not specified in config
    othersamples = options["general"]["othersamples"]

    # If user wants to split by all samples/not split by samples
    # then cannot use othersamples as well
    if othersamples and samples == ["**"]:         logger.error("Cannot have othersamples when using ** in Samples list")
    if othersamples and samples == cfg_samples:    logger.error("Cannot have othersamples when using * in Samples list")

    # If user wants to split by all systematics/not split by systematics
    # then cannot use othersysts as well
    if othersysts and systematics == ["**"]:         logger.error("Cannot have othersysts when using ** in Systematics list")
    if othersysts and systematics == cfg_systs:      logger.error("Cannot have othersysts when using * in Systematics list")

    # =========================================================================
    # =================== To Split or not to Split ============================
    # =========================================================================
    # A list of samples to split or not-to-split by systematics
    samples_to_split_by_sys = options["general"]["splitbysys"]
    sample_not_split_by_sys = options["general"]["notsplitbysys"]

    if sample_not_split_by_sys is None:
        sample_not_split_by_sys = ['-99']

    # if not specified, all samples are split by systematics provided in config
    if samples_to_split_by_sys is None:
        samples_to_split_by_sys = [s.replace("*",".*") for s in samples if s not in sample_not_split_by_sys] + ['othersamples'] if othersamples else [s.replace("*",".*") for s in samples if s not in sample_not_split_by_sys]


    # =========================================================================
    # ================== Start processing the jobs!  ===========================
    # =========================================================================

    def get_combos():
        for region in regions:
            # Make sure the requested region exists in the trex-fitter config
            # ** will be valid as it wildcards everything
            reg_valid = x_in_trex_xs(region, cfg_regions)
            if not reg_valid:   logger.error(f"Region {region} does not match any regions in trex-fitter config")

            for sample in (samples+['othersamples'] if othersamples else samples):
                # Make sure the requested samples exists in the trex-fitter config
                # ** will be valid as it wildcards everything
                samp_valid = x_in_trex_xs(sample, cfg_samples)
                if not samp_valid and sample != 'othersamples': logger.error(f"Sample {sample} does not match any samples in trex-fitter config")

                # Workout if we need to split this sample by systematics or not
                sample_rgx = sample.replace("*",".*")
                if sample != "**":
                    syst_split_this_sample      =  (any(rgx_match(sample_rgx, split_samp) for split_samp in samples_to_split_by_sys))
                    dont_syst_split_this_sample =  (any(rgx_match(sample_rgx, not_split_samp) for not_split_samp in sample_not_split_by_sys))
                    sample_is_data              =  (any(rgx_match(sample_rgx, datum) for datum in data))

                    if syst_split_this_sample and dont_syst_split_this_sample and sample :
                        logger.error(f"Sample {sample} cannot be both split and not split by systematics")

                    # This will split othersamples by default as it is added in samples_to_split_by_sys if user did not specify samples_to_split_by_sys
                    syst_split_this_sample = syst_split_this_sample and (not dont_syst_split_this_sample) and (not sample_is_data)
                else:
                    syst_split_this_sample = True

                # If we will split by systematics, then we need to loop over all requested systematics
                if syst_split_this_sample:  syst_loop = systematics + ["othersysts"] if othersysts else systematics
                # If we will not split by systematics, then we need to loop over only one systematic named ** and get_jobline will take care of the rest
                else:   syst_loop = ["**"]

                if cfg_normfactors is not None: syst_loop += cfg_normfactors

                for systematic in syst_loop:
                    # Make sure the requested systematic exists in the trex-fitter config
                    syst_valid = x_in_trex_xs(systematic, cfg_systs) or x_in_trex_xs(systematic, cfg_normfactors)
                    if not syst_valid and systematic == 'NONE':  pass
                    elif not syst_valid and systematic != 'othersysts' : logger.error(f"Systematic {systematic} does not match any systematics in trex-fitter config")

                    yield (region, sample, systematic)

    for region, sample, systematic in get_combos():

        suf = get_sufffix(region, sample, systematic ,removefromsuff)
        if options['general']['onlykeep'] is not None:
            if not any(keepsuff in suf for keepsuff in options['general']['onlykeep']): continue

        if bootstrap.n > 0:
            for bs_idx in range(bootstrap.n):
                if injectedfitfile is not None:
                    injectedfitfile = injectedfitfile.replace('.txt', f'_Bootstrap__{bs_idx}.txt')

                if bootstrap.step not in ['n','b']:
                    # Check if we need a botostrap idx in suffix based on which step is being ran (need it in not n)
                    suf = get_sufffix(region, sample, systematic ,removefromsuff, bs_idx=bs_idx)

                job = get_jobline(sample, region, systematic, exclude_str, excl_samp, excl_syst, suf, bs_idx=bs_idx, injectedfitfile=injectedfitfile, runranking=run_ranking)
                jobs.append(job)
        else:
            job = get_jobline(sample, region, systematic, exclude_str, excl_samp, excl_syst, suf, bs_idx=None, injectedfitfile=injectedfitfile,runranking=run_ranking)
            jobs.append(job)

        expected_files = get_histo_name(job, options["general"]["trexconfig"], cfg_regions, cfg_samples, jobname, excl_samp, exclude)
        all_exp_files.append(expected_files)


    jobs = sorted(list(set(jobs)))
    return jobs, all_exp_files, regions


def get_sufffix(
    region:     str,
    sample:     str,
    systematic: str,
    removefromsuff: [str, None],
    bs_idx: int = None,
) -> str:
    '''
    Method that constructs the trex-fitter command suffix
    Args:
        suf_sample (str):       The name of sample from job config
        region (str):       The name of region from job config
        systematic (str):   The name of systematic from job config
        removefromsuff (str):   Part of the suffix to be removed
        bootstrap (Bootstrapper):   The object containing the bootstrap settings

    Returns:
        The suffix part of the trex-fitter job
    '''
    # The object name to be used in suffix (no *)
    suf_sample = sample.replace("*","")
    suf_syst = systematic.replace("*", "").rstrip('_')
    suf_reg = region.replace("*", "")
    suf_bootstrap = f"Bootstrap__{bs_idx}"

    # The template part
    savesuffix = 'SaveSuffix=_'
    suf = f'{savesuffix}!REGION!_!SAMPLE!_!SYSTEMATIC!_!BOOTSTRAP!'
    # If either region/sample/systematic are specified to be removed from suffix, remove them
    if  'SAMPLE'     in removefromsuff or  sample    == '**':  suf = suf.replace('!SAMPLE!','')
    if  'REGION'     in removefromsuff or region     == '**':  suf = suf.replace('!REGION!','')
    if  'SYSTEMATIC' in removefromsuff or systematic == '**':  suf = suf.replace('!SYSTEMATIC!','')
    if  bs_idx is None:                                          suf = suf.replace('!BOOTSTRAP!','')

    # Replace template parts with values
    suf = suf.replace('!REGION!',   suf_reg)
    suf = suf.replace('!SAMPLE!',   suf_sample)
    suf = suf.replace('!SYSTEMATIC!',suf_syst)
    suf = suf.replace('!BOOTSTRAP!',suf_bootstrap)

    # Cleanup the "_"
    compiled = re.compile(r'_{2,}')
    suf = compiled.sub('_', suf)
    suf = suf.rstrip(r'_')
    # If all parts are removed, no suffix should be used
    if suf+'_' == savesuffix:  suf = ''
    return suf

def get_jobline(
    sample: str,  region: str, systematic: str,
    excluded: str, excl_samp: str, excl_syst: str,
    suf: str,
    bs_idx: int = None, injectedfitfile: str = None,
    runranking: bool = False,
) -> str:
    '''
    Method that constructs the trex-fitter command object specification
    i.e. Regions=X:Samples=Y:....
    Args:
        sample (str):       The name of sample from job config
        region (str):       The name of region from job config
        systematic (str):   The name of systematic from job config
        excluded (str):     The comma-separated string of excluded things from job config
        excl_samp (str):    The comma-separated string of all samples not in excluded, used for othersamples
        excl_syst (str):    The comma-separated string of all systs not in excluded, used for othersamples
        suf (str):          The ready-made suffix part of the job
        bs_idx (int):       The bootstrap index, if any
        injectedfitfile (str):  The path to the injected fit results, if any
        runranking (bool):  Whether to run ranking or not

    Returns:
        The object specification part of the trex-fitter job
    '''

    job = f'!REGION!_!SAMPLE!_!SYSTEMATIC!_!EXCLUDE!_!BOOTSTRAP!_!INJECTION!_!SUFFIX!'

    # Remove parts of the templates based on settings
    if sample == '**' or sample == '':  job = job.replace('!SAMPLE!','')
    if region == '**' or sample == '':  job = job.replace('!REGION!','')
    if systematic == '**' or systematic == '':                job = job.replace('!SYSTEMATIC!','')
    if excluded == '':                                        job = job.replace('!EXCLUDE!','')
    if bs_idx is None:                                        job = job.replace('!BOOTSTRAP!','')
    if injectedfitfile is None:                               job = job.replace('!INJECTION!','')
    if not runranking:                                        job = job.replace('!RANKING!','')
    if sample == 'othersamples':        job = job.replace('!SAMPLE!','!EXCLUDESAMPLES!')
    if systematic == 'othersysts':      job = job.replace('!SYSTEMATIC!','!EXCLUDESYSTEMATIC!')


    # Order job string such that the excludes come last
    job_pieces, suf_piece = job.split('_')[:-1],  job.split('_')[-1]
    exclude_pieces, non_exclude_pieces = [], []
    for piece in job_pieces:
        if "EXCLUDE" in piece:  exclude_pieces.append(piece)
        else:   non_exclude_pieces.append(piece)
    job = ':'.join(non_exclude_pieces+exclude_pieces+[suf_piece])

    # The order here matters (thanks to the ordering, i.e. Excludes come last in job string)
    # If we have 3 excludes, replace with 1 Exclude= statement
    job = job.replace('!EXCLUDESAMPLES!:!EXCLUDESYSTEMATIC!:!EXCLUDE!', f'Exclude={excl_samp},{excl_syst},{excluded}')
    # If the above template doesn't exist, we then can have any combo of 2 excludes
    # We replace 2 excludes with 1 Exclude= statement
    job = job.replace('!EXCLUDESAMPLES!:!EXCLUDE!', f'Exclude={excl_samp},{excluded}')
    job = job.replace('!EXCLUDESYSTEMATIC!:!EXCLUDE!', f'Exclude={excl_syst},{excluded}')
    job = job.replace('!EXCLUDESAMPLES!:!EXCLUDESYSTEMATIC!', f'Exclude={excl_samp},{excl_syst}')
    # If the above fail, then we have only one exclude
    # Order doesn't matter, just replace the exclude appropriately
    job = job.replace('!EXCLUDE!', f'Exclude={excluded}')
    job = job.replace('!EXCLUDESAMPLES!', f'Exclude={excl_samp}')
    job = job.replace('!EXCLUDESYSTEMATIC!', f'Exclude={excl_syst}')

    # Prepare the different parts that will enter the job
    job_reg =  f"Regions={region}"
    job_samp = f"Samples={sample}" if sample != 'othersamples' else ''
    job_syst = f"Systematics={systematic}" if systematic != 'othersamples' else ''
    job_inj =  f"NPValuesFromFitResults={injectedfitfile}"
    job_boot = f"BootstrapIdx={bs_idx}"

    # If any part of the job template is not there
    # replace will not do anything
    job = job.replace('!REGION!', job_reg)
    job = job.replace('!SAMPLE!', job_samp)
    job = job.replace('!SYSTEMATIC!', job_syst)
    job = job.replace('!SUFFIX!', suf)
    job = job.replace('!BOOTSTRAP!', job_boot)
    job = job.replace('!INJECTION!', job_inj)
    if runranking:  job = job.replace('Systematics=', 'Ranking=')

    # Cleanup of ":"
    compiled = re.compile(r':{2,}')
    job = compiled.sub(':', job)
    job = job.rstrip(r':')
    job = job.lstrip(r':')


    return job+'\n'

def get_histo_name(
    jobline:     str,
    trex_config: str,
    cfg_regions: List[str],
    cfg_samples: List[str],
    cfg_jobname: str,
    excl_samp:   str,
    excluded:    List[str],
) -> List[str]:
    '''
    Method to extract the names of all the ROOT files that should be produced by
    the arguments constructed using this script along with the names of all nominal
    histograms expected within a given file.

    Args:
        jobline (str): The part of the trex-fitter job with flags
        trex_config (str):  The path to trex-fitter config
        cfg_regions (List[str]): The list of regions found in config file
        cfg_samples (List[str]): The list of samples found in config file
        cfg_jobname (str): The job name specified in the config
        excl_samp   (str): Comma-separated string of all samples (pyregex) that should be excluded when using othersamples
        excluded    (List[str]): A list of excluded objects (pyregex) that should be excluded from all jobs

    Returns:
        A list with elements being outfile_i/nominal_j for i and j indexing the output files
        expected and number of nominals expected, respectively
    '''

    x_excluded = lambda x, y: any(re.search('^'+excl+'$', x) is not None for excl in y)
    # ================================== Breakdown the job  =====================================
    regions = next((piece for piece in jobline.split(':') if "Regions=" in piece), '').split(',')
    samples = next((piece for piece in jobline.split(':') if "Samples=" in piece), '').split(',')
    excludes = next((piece for piece in jobline.split(':') if "Exclude=" in piece), '')
    # ================================== List of Samples  =====================================
    # If no samples flag was found, then we have to do some work to get names of nominal histograms
    if samples == ['']:
        # Then this either othersamples, or we are not splitting by sample
        if excludes == '':
            # Then we are not splitting by sample for sure (and no excluded stuff)
            samples = cfg_samples # All config samples should be in ROOT file
        elif excl_samp in excludes:
            excl_samp_lst_rgx = [es.replace('*','.*') for es in excl_samp.split(',')]
            # Then this is othersamples job (cannot use suffix in case sample removed)
            samples = [s for s in cfg_samples if not x_excluded(s, excl_samp_lst_rgx+excluded)]
        elif excl_samp not in excludes:
            # Then we are not splitting by sample, but some samples may be excluded
            samples = [s for s in cfg_samples if not x_excluded(s, excluded)]

    else: # if Samples are defined, they can be regex, so exact sample names are extracted from trex cfg
        samples_regex = [s.replace('*','.*').replace("Samples=","") for s in samples]
        samples = [cs for cs in cfg_samples if any(re.search(f'^{s}$', cs) is not None for s in samples_regex)]

    # ================================== Suffix used =====================================
    suffix = next((piece for piece in jobline.split(':') if "SaveSuffix=" in piece), None)
    if suffix is not None:  suffix = suffix.replace("SaveSuffix=","").strip()
    else:   suffix = ''

    # ================================== List of regions =====================================
    if regions == ['']:
        # Then we are not splitting by region
        if excludes == '':
            # Then we are not excluding any regions
            regions = cfg_regions # All config region should be in ROOT file
        else:
            # Then we are not splitting by region, but some regions may be excluded
            regions = [r for r in cfg_regions if not x_excluded(r, excluded)]
    else:  # if regions are defined, they can be regex, so exact region names are extracted from trex cfg
        regions_regex = [r.replace('*','.*').replace("Regions=","") for r in regions]
        regions = [cr for cr in cfg_regions if any(re.search(f'^{r}$', cr) is not None for r in regions_regex)]

    # ================================== Workout expected =====================================
    expected_outfiles = []
    # Loop over exact region names
    for region in regions:
        # Output files are defined by region, jobname and suffix
        histofile = f"{cfg_jobname}_{region}_histos{suffix}.root"
        # Loop over exact sample names
        for sample in samples:
            # Nominal histogram name is defined by region and sample
            nominal = f'{region}_{sample}'.lstrip('_')
            # Save the root file name followed by the nominal histo name
            # This produceds N = Number of nominals entries per output file
            expected_outfiles.append(f'{histofile}!!{nominal}!!trex-fitter n {trex_config} {jobline.strip()}')

    return expected_outfiles


def x_in_trex_xs(x: str, trex_xs: str) -> bool:
    '''
    Method to make sure the requested object exists in the trex-fitter config
    Args:
        x (str): The object name specified in job config
        trex_xs (str):  All the trex-fitter config options for this object type
    Returns:
        valid (bool)
    '''

    regex_x = x.replace("*",".*")
    valid = False
    for trex_x in trex_xs:
        if re.match(f'{regex_x}',f'{trex_x}') is not None:
            valid = True
            break
    return valid

def hupdate_file_maker(regions, job, outdir, histfolder, bootstrap):

    def getline(bootstrap_idx):
        lines = []
        for region in regions:
            folder = histfolder if '/Histograms/' in histfolder else histfolder+'Histograms/'
            if bootstrap_idx is not None:
                # Check nominal histogram file for region exists
                nom_histos_file = f"{bootstrap.nom_hfolder}/{bootstrap.nom_jobname}_{region}_histos.root"
                if not os.path.exists(nom_histos_file):
                    logger.error(f'Nominal histogram file {nom_histos_file} does not exist')
                line = f'hupdate.exe {folder}/{job}_{region}_histos_{bootstrap.param}_{bootstrap_idx}.root {folder}/{job}_{region}_histos_*_{bootstrap.param}_{i}.root & \n'
            else:
                line = f'hupdate.exe {folder}/{job}_{region}_histos.root {folder}/{job}_{region}_histos_*.root & \n'

            lines.append(line)

        return lines

    all_lines = []
    if bootstrap.n > 0:
        fname = f'{outdir}hupdate_Bootstrap.sh'
        for i in range(bootstrap.n):
            lines = getline(i)
            for line in lines:
                if i == bootstrap.n-1 and line == lines[-1]:
                    line = line.replace('&&','&')
                all_lines.append(line)
    else:
        fname = f'{outdir}hupdate.sh'
        lines = getline(None)
        for line in lines:
            all_lines.append(line)

    chunks = np.array_split(all_lines, np.ceil(len(all_lines)/100))
    for i, chunk in enumerate(chunks):
        chunk_fname = fname.replace('.sh',f'_{i}.sh')
        hupdate_file = open(chunk_fname, 'w')
        for line in chunk:
            hupdate_file.write(line)
        hupdate_file.close()

if __name__ == '__main__':  main()