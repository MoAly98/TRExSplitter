# TRExSplitter
Simple utility code to help parallelising trex-fitter jobs. This utlity is independent of which step of TRExFitter you are trying to parallelise, however, the command line arguments supported are:

- `Regions`
- `Samples`
- `Systematics`
- `SaveSuffix`
- `BootstrapIdx`
- `NPValuesFromFitResults`

So you can parallelise jobs that require the above command line options only. If you want to parallelise jobs that require a different command line option, please open and issue on this repository.

## How this works

The utility will start by traversing your `TRExFitter` configuration file, working out all the available samples, regions and systematics. Then based on how you request to split your jobs, the utlity will build an argument list with `N` jobs split in the way you request.

## Getting this utility

You should clone this repository to get the code:

```
git clone git@github.com:MoAly98/TRExSplitter.git
```

## Setting up the utility
You will need to have some python packages to run this code. You can install these in your favourite python environment with

```bash
pip install -r requirements.txt
```

Note that if you are on `lxplus`, you should run the above command in a virtual environment (`conda`/`venv`).

## Running the utility

Often one is looking to optimize the processing of `TRExFitter` commands, such as running the `n`-step or `r`-step. This repository can help
in preparign a list of jobs for the user that they can run using their favourite method later on.

To run the code, execute

```bash
python jobs/trex_job_creator.py  <your_creation_config>
```

where `<your_creation_config>` is a `JSON` file containing the settings that specify how you would like to split your jobs, with a lot of flexibility.

The above will run the `n` step. To run a different step, use:

```bash
python jobs/trex_job_creator.py <your_creation_config> --step <MyStep>
```

To suffix your output files, you should use the `-s` option.

### The Configuration File

The configuration file is `json` file where you can specify flexibly how to split your `trex-fitter` jobs by regions/samples/systematics. The setting names in the configuration file are not case-senstive.

The configuration file looks like:

```json
{
    "Regions": ...,
    "Samples": ...,
    "Systematics": ....,
    "General": {
        "TRExConfig": ...,
        "OutDir":....,
        ....
    }
}
```

where the fields shown in the schematic are the required fields in the config. More details are in [this section](#required-parts).

An example of a configuration file can be found in this repository (`creation_config.json`).


#### Required Parts

The following fields are needed to specify what you want to parallelise. You can split by all/any of the main components of a `TRExFitter` config (samples, regions and systematics). If you want to completely parallelise one component (e.g. only run one region per job), use the `Full Split` syntax. If you don't want to split a component at all (e.g. run all regions in each job),  use the `No Split` syntax.

| Settings| Description| Full Split | No Split  |
|---------|------------| ---------- | ----------|
| Regions | The list of regions to split by         | \* | \*\* |
| Samples | The list of samples to split by         | \* | \*\* |
| Systematics | The list of systematics to split by | \* | \*\* |

The required general settings are:

| Settings| Description|
|---------|------------|
| TRExConfig | The full path to the relevant `TRExFitter` configuration file |
| OutDir     | The full path to where the job list should be saved |

#### Optional Parts

The (optional) parallelisations settings:

| Settings| Description|
|---------|------------|
| Exclude | A list of Regions/Samples/Systematics that should be excluded from **ALL** jobs |


The optional general settings:

| Settings| Description|
|---------|------------|
| othersysts   | Boolean to determine if systematics not specified in the `Systematics` field should be combined into one job (except Excluded). The suffix of these jobs will contain `_othersysts`. |
| othersamples | Boolean to determine if samples not specified in the `Samples` field should be combined into one job (except Excluded). The suffix of these jobs will contain `_othersamples`. |
| SplitBySys   | List of samples to split by systematics. This is useful if you in general don't want to split by systematics, except for some samples.  |
| NotSplitBySys| List of samples that should **NOT** be split by systematics (e.g. Fakes) |
| MaxJobsPerChunk |  The maximum number of jobs to be dumped into 1 file. This is useful if you are going to produce many parallel jobs and want to run `N` jobs at a time.|
| RemoveFromSuff | [REGION, SAMPLE, SYSTEMATIC] which specifies if the suffix should not contain the name of the chosen component.|
| OnlyKeep | A filter for jobs based on the suffix. Use this if you only want to keep jobs with a certain suffix (e.g. split by systematics and only keep get other systematic jobs) |

The optional bootstrapping settings:
**Use these if you want to parallelise bootstrapping replica histgoram-creation or fits**

| Settings| Description|
|---------|------------|
| Bootstrap      | The number of bootstrap replicas you aim to produce |
| BootstrapParam | The systematic/sample being bootstrapped |

The optional injection settings:

**Use these if you want to parallelise injecting fit results into a new fit at scale, e.g. from bootstrapping**

| Settings| Description|
|---------|------------|
| InjectFromCfg  | `TRExFitter` Configuration file that produces fit results to be injected in your new fit |

### The Output

The code will produce 3 types of files:

1. Job list: `arguments_X.txt` files will be produced (1 if not chunking, multiple if chunking with `MaxJobsPerChunk`). This contains a line-by-line list of the `TRExFitter` jobs that you want to run.

2. Histogram Merging Script: `hupdate_X.sh` file will be produced which contain `hupdate.exe` executable calls to merge all histograms belonging to one region into the same file.

3. Histogram Summary: `histogram_files.txt` will be produced If you are running `n`/`b` step. They summarise which samples, regions and systematics have been processed by each of the jobs created. These can be used to debug `TRExFitter` histogramming outputs in external scripts.
