# Enforce explicit decision of output directory
if '--directory' not in sys.argv and '-d' not in sys.argv and 'assembly_workdir' not in config:
    print(sys.argv)
    raise KeyError('Assembly working directory not specified.')

# Get information on cores per node
try:
    CORES = int(config['cores'])
except KeyError:
    print('Cores not specified in config file.')
    raise
except TypeError:
    print('Cores cannot be interpreted as integer')
    raise
if CORES < 1:
    raise ValueError('Cores must be a positive integer')
config['cores'] = CORES

# Parse datalines
samples, experiments, runs = dict(), dict(), dict()
for line in config['data']:
    fields = line.split()

    if len(fields) not in (4, 5):
        raise ValueError('Data field {} needs 4 or 5 fields.'.format(line))

    run, experiment, sample, *paths = fields

    # No cross-category non-uniqueness in names
    if (run in experiments or run in samples or run in runs):
        raise ValueError('Run {} is not unique in data.'.format(run))
    if (experiment in runs or experiment in samples):
        raise ValueError('Experiment {} is also a run/sample'.format(experiment))
    if (sample in runs or sample in experiments):
        raise ValueError('Sample {} is also a run/experiment'.format(sample))

    # File exists
    for path in paths:
        if not os.path.isfile(path):
            raise FileNotFoundError(path)

    # Add them to dictionaries
    if sample in samples:
        samples[sample].append(experiment)
    else:
        samples[sample] = [experiment]
    if experiment in experiments:
        experiments[experiment].append(run)
    else:
        experiments[experiment] = [run]
    runs[run] = paths

# Not both SE and PE files in same experiment
for experiment, runlist in experiments.items():
    if len({len(runs[run]) for run in runlist}) != 1:
        raise ValueError('Experiment {} has both SE and PE runs.'.format(experiment))

# No experiment is in two samples (runs are ensured unique above)
if len(set().union(*samples.values())) != sum(len(v) for v in samples.values()):
    raise ValueError('No experiments can be in two different samples.')

# Create IS_PE for easy lookup which experiments are paired end
IS_PE = {exp: len(runs[runlist[0]]) == 2 for exp, runlist in experiments.items()}

# Finally add them to the config
config['samples'] = samples
config['experiments'] = experiments
config['runs'] = runs

# Get information on maximum number of mismatches
try:
    MM = int(config['adapterremoval_mm'])
except KeyError:
    print('Max mismatches (mm) not specified in config file.')
    raise
except ValueError:
    print('Max mismatches (mm) cannot be interpreted as integer')
    raise
if MM < 0:
    raise ValueError('Max mismatches (mm) must be a non-negative integer')

# Get information on DNA adapters, if None, do not specify parameter
try:
    ADAPTER1 = config["adapterremoval_adapter1"]
    ADAPTER2 = config["adapterremoval_adapter2"]
    ADAPTERLIST = config["adapterremoval_adapterfile"]
except KeyError:
    print('Following keys must exists in config file: '
          'adapterremoval_adapter1, '
          'adapterremoval_adapter2, '
          'adapterremoval_adapterfile')
    raise

if ADAPTER1 and ADAPTER2:
    if ADAPTERLIST:
        print('Specify either adapters or adapterlist')
        raise KeyError('Specify either adapters or adapterlist')

    ADAPTERARGS_PE = "--adapter1 {} --adapter2 {} ".format(pcr1, pcr2)
    ADAPTERARGS_SE = "--adapter1 {} ".format(pcr1)
    ADAPTERARGS_LIST = ""
    
else:
    if ADAPTERLIST:
        ADAPTERARGS_LIST = "--adapter-list {}".format(ADAPTERLIST)
    else:
        ADAPTERARGS_LIST = ""

    ADAPTERARGS_PE = ""
    ADAPTERARGS_SE = ""

# Get information on assembler kmers
try:
    ASSEMBLER_KMERS = config['assembler_kmers']
    kmers = [int(k) for k in ASSEMBLER_KMERS.split(',')]
except KeyError:
    print('Assembler kmers not specified in config file.')
    raise
except TypeError:
    raise ValueError('Assembler kmers cannot be interpreted as comma-seperated string of integers')
else:
    if not kmers == sorted(kmers) or not all((k%2 and k>10 and k<128) for k in kmers):
        raise ValueError('Assembler kmer sizes is not a sorted list of odd '
                         'positive ascending integers 10<k<128.')
    del kmers

# Get information on choice of assembler
try:
    ASSEMBLER = config['assembler'].lower()
except KeyError:
    print('Choice of assembler not specified in config file')
    raise
if ASSEMBLER not in ('megahit', 'spades'):
    raise ValueError('Assembler must be either Megahit or Spades')
    
# If megahit, get information on mink_1pass
if ASSEMBLER == 'megahit':
    try:
        KMIN_1PASS = config['megahit_mink_1pass']
    except KeyError:
        print('Megahit kmer 1pass not specified in config file.')
        raise
    if not isinstance(KMIN_1PASS, bool):
        raise TypeError('Megahit kmin 1pass must be boolean (true or false). '
        'Be sure to write it in lower case letters, JSON-style.')

# Else, assert it's SPAdes
else:
    assert ASSEMBLER == 'spades'

# Get information on programs to be used
names = ['AdapterRemoval', ASSEMBLER]
keys = tuple(name.lower() + '_path' for name in names)

for name, key in zip(names, keys):
    if key not in config:
        raise KeyError('{} path not specified in config file'.format(name))

    path = config[key]
    if not os.path.isfile(path):
        raise FileNotFoundError('{} is not an existing file'.format(path))

# Enforce Python version and only one PE experiments sample when using SPAdes
# (as of Spades 3.10, the exp restriction applies when using the --meta flag)
if ASSEMBLER == 'spades':
    # Check Python version
    allowed_versions = ((2,4), (2,5), (2,6), (2,7), (3,2), (3,3), (3,4), (3,5), (3,6))
    present_version = (sys.version_info.major, sys.version_info.minor)
    
    if present_version not in allowed_versions:
        print('Warning: Python version {}.{} not supported by SPAdes v 3.10. '
                      'Supported versions are 2.(4-7) and 3.(2-6). '
                      'Attempting to use python 3.5 directly from executable.'.format(
                      present_version[0], present_version[1]))
    
        # Attempt to use alternative python 3.5 executable:
        try:
            PYTHON_PATH = config['alternative_python_path']
        except KeyError:
            print('Alternative python 3.5 path needed but not specified in config file.')
            raise
        if not os.path.isfile(PYTHON_PATH):
            raise FileNotFoundError('Alternative python 3.5 path {} '
                                    'needed but not found.'.format(PYTHON_PATH))
    
    else: # present version is okay:
        PYTHON_PATH = sys.executable
    
    del allowed_versions, present_version
    
    # Check that all experiments are paired end
    if not all(IS_PE.values()):
        print('metaSPAdes v 3.10 does not support single-end libraries.')
        raise ValueError

################## CREATE OUTPUT LISTS REPRESENTING MERGING/TRIMMING STEP #########

MERGE_OUTPUT = list()
TRIM_OUTPUT = list()
for experiment in config['experiments'].keys():
    if IS_PE[experiment]:
        MERGE_OUTPUT.append('merged/{}.fw.fastq.gz'.format(experiment))
        MERGE_OUTPUT.append('merged/{}.rv.fastq.gz'.format(experiment))
        TRIM_OUTPUT.append('trim/{}.fw.fastq.gz'.format(experiment))
        TRIM_OUTPUT.append('trim/{}.rv.fastq.gz'.format(experiment))

    else:
        MERGE_OUTPUT.append('merged/{}.fastq.gz'.format(experiment))
        TRIM_OUTPUT.append('trim/{}.fastq.gz'.format(experiment))

print('Parsed assembly config file.')
