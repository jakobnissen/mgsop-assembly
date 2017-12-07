"Submission script for Jakob's pipeline. See readme for details."

import argparse
import sys
import os
import subprocess
import json

if sys.version_info < (3, 5):
    raise OSError('Requires Python v. 3.5 or later.')
    
# Invoke argparse
usage = """python assemble.py config group [-o outdir]
[-j maxjobs] [-c a=b x=y ... ] [-t targetfile] [--commit]"""

parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    usage=usage, add_help=False)

# Positional arguments
parser.add_argument('configfile', help='configuration file')
parser.add_argument('group', help='group for Computerome queue')

# Optional arguments
options = parser.add_argument_group('options')
options.add_argument('-h', action='help', help='show this help message and exit')
options.add_argument('-o', dest='outdir', help='Working directory')
options.add_argument('-n', dest='name', help='name of your analysis run.')
options.add_argument('-j', dest='maxjobs', type=int, default=20,
    help='maximal number of simultaneous jobs to run [20]')
options.add_argument('-t', dest='target', default='',
    help='target file [default set by snakefile]')
options.add_argument('-c', dest='config', nargs='+', metavar='k=v', default=[],
    help='key=val pairs for configuration')
options.add_argument('--commit', action='store_true',
    help='do not dry-run, execute snakefile')

# If no arguments, print help
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit()

args = parser.parse_args()

# Set outdirstring
if args.outdir is None:
    args.outdirstring = ''
else:
    args.outdirstring = '-d {}'.format(args.outdir)

# Set dryrun
if args.commit:
    args.dryrun = ''
else:
    args.dryrun = 'n'
    
# Set snakefile path
args.snakefile = os.path.join(sys.path[0], 'Snakefile')

# Test if config file exists:
if not os.path.isfile(args.configfile):
    raise FileNotFoundError(args.configfile)
args.configfilestring = '--configfile {}'.format(args.configfile)
    

validconfigs = {'cores', 'adapterremoval_mm', 'adapterremoval_adapter1',
                'adapterremoval_adapter2', 'adapterremoval_adapterfile',
                'assembler', 'assemble_contigs', 'assembler_kmers',
                'megahit_mink_1pass'}

# Set configstring
if len(args.config) != 0:      
    for config in args.config:
        key, equal, value = config.partition('=')
        if not equal or not value:
            raise ValueError('Config "{}" does not follow pattern k=v'.format(config))
        
        if key not in validconfigs:
            raise KeyError('Config key "{}" is not valid. '
                           'Choose among {}'.format(key, ', '.join(validconfigs)))
    
    args.configstring = '-C ' + ' '.join(args.config)
else:
    args.configstring = ''

# Fill in all fields execute command
template = """snakemake -s {0.snakefile} -c "qsub -l \
walltime={{params.walltime}},nodes=1:ppn={{params.ppn}},mem={{params.mem}} \
-A {0.group} -W group_list={0.group}" {0.configfilestring} -pk{0.dryrun}j \
{0.maxjobs} {0.outdirstring} {0.target} -w 60 {0.configstring}"""

command = template.format(args)

subprocess.run(command, shell=True, check=True)
