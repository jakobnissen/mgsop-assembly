# Load the general configuration file 
configfile: srcdir('config.yaml')

include: 'scripts/parseconfig.snakefile'

# Set working directory to that found in the passed configfile
# This is overwritten if output directory is given by command line
workdir: config['assembly_workdir']

# This helps with network latency.
shell.prefix('sleep 10; ')

rule all:
    input:
        # Merge runs to experiments
        #MERGE_OUTPUT
        # Trim using AdapterRemoval
        #TRIM_OUTPUT
        # Assemble
        expand('assembly/{sample}.fna', sample=config['samples'])

########################## MERGING OF RUNS ###########################
# Runs are defined as repititions of identical library setups.
# They can therefore be concatenated together to one file immediately.
# This is only done if multiple runs are found per experiment, whether
# that is the case is detected by the input functions to adapterremoval
######################################################################

def merge_runs_pe_input_fw(wildcards):
    return [config['runs'][run][0] for run in config['experiments'][wildcards.experiment]]
    
def merge_runs_pe_input_rv(wildcards):
    return [config['runs'][run][1] for run in config['experiments'][wildcards.experiment]]

rule merge_runs_pe:
    input:
        forward = merge_runs_pe_input_fw,
        reverse = merge_runs_pe_input_rv 
    output:
        forward = temp('merged/{experiment}.fw.fastq.gz'),
        reverse = temp('merged/{experiment}.rv.fastq.gz')
    params: 
        walltime = 3600,
        ppn = 1,
        mem = '1gb'
    shell:
        'cat {input.forward} > {output.forward}; '
        'cat {input.reverse} > {output.reverse}'

def merge_runs_se_input(wildcards):
    return [config['runs'][run][0] for run in config['experiments'][wildcards.experiment]]

rule merge_runs_se:
    input: merge_runs_se_input
    output: temp('merged/{experiment}.fastq.gz')
    params:
        walltime = 3600,
        ppn = 1,
        mem = '1gb'
    shell: 'cat {input} > {output}'

####################### ADAPTERREMOVAL ########################
# For function of input functions: See MERGING OF RUNS rule.
###############################################################

def adapterremoval_input_fw(wildcards):
    # If more than one run exist per experiment
    if len(config['experiments'][wildcards.experiment]) > 1:
        return rules.merge_runs_pe.output.forward
    
    else: # If only one run per experiment
        runname = config['experiments'][wildcards.experiment][0]
        return config['runs'][runname][0]
        
def adapterremoval_input_rv(wildcards):
    # If more than one run exist per experiment
    if len(config['experiments'][wildcards.experiment]) > 1:
        return rules.merge_runs_pe.output.reverse
    
    else: # If only one run per experiment
        runname = config['experiments'][wildcards.experiment][0]
        return config['runs'][runname][1]

def adapterremoval_input_se(wildcards):
    # If more than one run exist per experiment
    if len(config['experiments'][wildcards.experiment]) > 1:
        return rules.merge_runs_se.output
    
    else: # If only one run per experiment
        runname = config['experiments'][wildcards.experiment][0]        
        return config['runs'][runname][0]

rule adapterremoval_pe:
    input:
        forward = adapterremoval_input_fw,
        reverse = adapterremoval_input_rv
    output:
        forward = 'trim/{experiment}.fw.fastq.gz',
        reverse = 'trim/{experiment}.rv.fastq.gz'
    threads: 2
    log: 'log/trim/{experiment}.trim.log'
    params: 
        mm = MM,
        walltime = 864000,
        ppn = 2,
        mem = '5gb'
    shell:
        '{config[adapterremoval_path]} '
        # Input files
        '--file1 <(zcat {input.forward}) --file2 <(zcat {input.reverse}) '
        # Output files
        '--output1 >(gzip > {output.forward}) --output2 >(gzip > {output.reverse}) '
        '--singleton /dev/null --discarded /dev/null --settings /dev/null '
        "2> {log} "
        # Adapters
        '{ADAPTERARGS_PE} {ADAPTERARGS_LIST} '
        # Other parameters:
        '--minlength 30 --trimns --trimqualities --minquality 2 '
        '--qualitybase 33 --qualitymax 43 --mm {params.mm} '
        '--threads {threads}'
        
rule adapterremoval_se:
    input: adapterremoval_input_se
    output: 'trim/{experiment}.fastq.gz'
    threads: 2
    log: 'log/trim/{experiment}.trim.log'
    params: 
        mm = MM,
        walltime = 864000,
        ppn = 2,
        mem = '5gb'
    shell:
        '{config[adapterremoval_path]} '
        # Input files
        '--file1 <(zcat {input}) '
        # Output files
        '--output1 >(gzip > {output}) '
        '--discarded /dev/null --settings /dev/null '
        "2> {log} "
        # Adapters
        '{ADAPTERARGS_SE} {ADAPTERARGS_LIST} '
        # Other parameters:
        '--minlength 30 --trimns --trimqualities --minquality 2 '
        '--qualitybase 33 --qualitymax 43 --mm {params.mm} '
        '--threads {threads}'
        
########################### ASSEMBLY ###################################
# The assembly input functions determine whether each experiment is
# single or paired ended. The single assembly_input_command then
# gathers information from the input functions and presents it as
# a string of files for the assembler to use.
# This is necessary to allow any combination of PE and SE files.

# Spades does not support more than one run as of Spades 3.11.0
# This is enforced in assembly_config.snakefile, and we need not concern us
# with this limitation here.

# Getting contigs is a separate rule because when the assembler crashes,
# you want to be  sure it's not the mv command failing, then snakemake 
# deleting the directory. It also makes the log files more transparent.
########################################################################

    
if ASSEMBLER == 'megahit':
    # We sort to make sure the fw and rv are in sync
    def megahit_input_fw(wildcards):
        return sorted(['trim/{}.fw.fastq.gz'.format(exp)
                       for exp in config['samples'][wildcards.sample] if IS_PE[exp]])

    def megahit_input_rv(wildcards):
        return sorted(['trim/{}.rv.fastq.gz'.format(exp)
                       for exp in config['samples'][wildcards.sample] if IS_PE[exp]])

    def megahit_input_se(wildcards):
        return ['trim/{}.fastq.gz'.format(exp)
                       for exp in config['samples'][wildcards.sample] if not IS_PE[exp]]

    def megahit_input_command(wildcards, input):
        pe_files = '-1 {} -2 {}'.format(','.join(input.fw), ','.join(input.rv)) if input.fw else ''
        se_files = ' -r {}'.format(','.join(input.se)) if input.se else ''
        return '{}{}'.format(pe_files, se_files)

    rule megahit:
        input:
            fw = megahit_input_fw,
            rv = megahit_input_rv,
            se = megahit_input_se
        output: temp('assembly/{sample}')
        threads: config['cores']
        log: 'log/assembly/{sample}.log'
        params:
            infiles = megahit_input_command,
            out_prefix = '{sample}',
            kmers = ASSEMBLER_KMERS,
            kmin_1pass = '--kmin-1pass' if KMIN_1PASS else '', 
            memory = 99000000000,            
            walltime = 864000,
            ppn = config['cores'],
            mem = '100gb'
        shell:
            '{config[megahit_path]} {params.infiles} '
            '-t {threads} -o {output} -m {params.memory} '
            '--out-prefix {params.out_prefix} '
            '--k-list {params.kmers} {params.kmin_1pass} '
            '2> {log}'
            
    rule get_megahit_contigs:
        input: 'assembly/{sample}'
        output: 'assembly/{sample}.fna'
        threads: 1
        params:
            filename = 'assembly/{sample}/{sample}.contigs.fa',
            walltime = 864000,
            ppn = 1,
            mem = '1gb',
        shell: 'mv {params.filename} {output}'

else: # Assembler is Spades
    def spades_input(wildcards):
        if IS_PE[wildcards.experiment]:
            return ['trim/{}.fw.fastq.gz'.format(wildcards.experiment),
                    'trim/{}.fw.fastq.gz'.format(wildcards.experiment)]
        else:
            return 'trim/{}.fastq.gz'.format(wildcards.experiment)

    def spades_input_command(wildcards):
        if IS_PE[wildcards.experiment]:
            return '-1 trim/{}.fw.fastq.gz -2 trim/{}.rv.fastq.gz'.format(
                wildcards.experiment, wildcards.experiment)
        else:
            return '-s trim/{}.fastq.gz'.format(wildcards.experiment)
            
    def get_spades_contigs_input(wildcards):
        return expand('assembly/{exp}', exp=config['samples'][wildcards.sample])

    def spades_contigs_command(wildcards):
        experiments = config['samples'][wildcards.sample]
        dest = 'assembly/{}.fna'.format(wildcards.sample)
        
        # If only one exp, move, do not concatenate.
        if len(experiments) == 1:
            source = 'assembly/{}/contigs.fasta'.format(experiments[0]) 
            return 'mv {} {}'.format(source, dest)
            
        else:
            sources = expand('assembly/{exp}/contigs.fasta', exp=config['samples'][wildcards.sample])
            return 'cat {} > {}'.format(' '.join(sources), dest)
        
    def get_spades_logs_filenames(wildcards):
        return expand('assembly/{exp}/spades.log', exp=config['samples'][wildcards.sample])

    rule spades:
        input: spades_input
        output: temp('assembly/{experiment}')
        threads: config['cores']
        params:
            python_exec = PYTHON_PATH,
            infiles = spades_input_command,
            memory = 99, # Must be a little less than cluster-params.
            kmers = ASSEMBLER_KMERS,
            walltime = 864000,
            ppn = config['cores'],
            mem = '100gb',
        shell:
            '{params.python_exec} {config[spades_path]} --meta '
            '{params.infiles} -k {params.kmers} '
            '-t {threads} -m {params.memory} -o {output}'
            
    rule get_spades_contigs:
        input: get_spades_contigs_input
        output:
            contigs = 'assembly/{sample}.fna',
            logs = 'log/assembly/{sample}.log'
        threads: 1
        params:
            command = spades_contigs_command,
            lognames = get_spades_logs_filenames,
            walltime = 864000,
            ppn = 1,
            mem = '1gb',
        shell:
            'cat {params.lognames} > {output.logs}; {params.command}'

def cleanup(directory):
    # Move the numerous "snakejob" files that Computerome can generate.
    files = [f for f in os.listdir(directory) if f.startswith('snakejob.')]
    snakejobdir = os.path.join(directory, 'snakejobs')

    try:
        os.mkdir(snakejobdir)
    except FileExistsError:
        if os.path.isdir(snakejobdir):
            pass
        else:
            print('A non-directory file called "snakejobs" already exists. Cannot clean.')
            return

    for file in files:
        source = os.path.join(directory, file)
        destination = os.path.join(snakejobdir, file)

        os.rename(source, destination)        

onsuccess:
    cleanup(os.getcwd())

onerror:
    cleanup(os.getcwd())
