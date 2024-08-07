#
# Based on the work of Morgante, P; Peverati, R.; SoftwareX, 2018, *submitted*
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation version 3.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

import csv
import re
import os
import glob
from operator import itemgetter
from datetime import timedelta
#
#    --- ACCDB WORKFLOW FILE ---
#
#   This Snakemake file is heavily modified from the original ACCDB Snakemake file. Please refer to the README.
# 

configfile: "config.yaml"

wildcard_constraints:
     molecule='[-\w+=.]+',
     method='[-\w+=.]+',
     dft='[-\w+=.]+',
     db='[-\w+=.]+',

# Input/output directories
DATABASES_DIR = config['DATABASES_DIR']
GEOMETRIES_DIR = config['GEOMETRIES_DIR']
ENERGIES_DIR=config['ENERGIES_DIR']
OUTPUTS_DIR=config['OUTPUTS_DIR']

DATABASES = config['DATABASES']

METHODS = config['METHODS']
RULES = config['RULES']

# TEMPL1 = 'qcengine.tmpl'
TEMPL1 = config['TEMPL1']

# NUMBER OF PROCESSORS/CORES/THREADS
nproc = config['PROC']

# A method for dealing with the different units being thrown around everywhere
# Energy of 1 Hartree in the designated unit. For Psi4, this is 1.0 since it outputs in Hartrees.
# If your software outputs in eV, for example, this would be 27.211396641308 Hartrees/eV
# Find this in the first ROW of https://ryutok.github.io/EnergyConversionTable/
units = config['UNITS']

# -------------------------------------------------------------------------------------------------------------------- #
#                                                  SUPPORT FUNCTIONS                                                   #
# -------------------------------------------------------------------------------------------------------------------- #

def mean(l): return sum(l) / len(l)
def meanabs(l): return mean([*map(abs, l)])
def rms(l): return mean([*map(lambda n: n**2.0, l)])**0.5
def maxabs(l): return max(map(abs, l))

def scalef(f, scale):
    """Returns a function that calls `f` and multiplies its output by `scale`."""
    return lambda *args, **kwargs: f(*args, **kwargs) * scale

def rule_matches(match, **attribs):
    """
        Evaluate whether the provided `**attribs` match the values in the `match` dict.
        For every item in `match`, the corresponding value in `**attribs` must equal it, or be contained in it if the
        value in `match` is a list. See the example `config.yaml`.
    """
    for k, v in match.items():
        if not k in attribs:
            return False
        elif isinstance(v, list):
            if not attribs[k] in v: return False
        elif match[k] != v:
            return False
    return True
def evaluate_rules(rules, **attribs):
    """
        Evaluates the given set of rules and returns the first matching rule, or None if nothing matches. Matches are
        evaluated by `rule_matches`. See the example `config.yaml`.
    """
    for rule in rules:
        if rule_matches(rule['match'], **attribs):
            return rule
    return None

def read_database_eval(db, rules=RULES, db_dir=DATABASES_DIR, csv_name='DatasetEval.csv'):
    """
        Parse a `DatasetEval.csv` file into a list of entries. Each entry has:
        * `db` - Database name
        * `rules` - Set of rules to determine calculation methods for datapoints
        * `db_dir` - Database directory to search in. Defaults to `DATABASES_DIR`
        * `csv_name` - The name of the value CSV. Defaults to `'IndValues.csv'`
    """
    path = os.path.join(db_dir, db, csv_name)
    results_list = []
    with open(path, newline='') as csvfile:
        evalreader = csv.reader(csvfile, dialect='excel')
        
        def evalRow(row):
            obj = {
              'name': row[0], # Name of the datapoint
              'dataset': '_'.join(row[0].split('_')[:-1]), # Name of the source dataset. TODO: This logic is a bit questionable
              'refval': float(row[-1]), # Reference energy in undetermined units -- See below
              'deps': row[2::2], # Dependencies; Raw calculations that make up this datapoint
              'coeffs': [*map(float, row[1:-1:2])], # Coefficients to apply to `deps` to make the datapoint `name`
            }
            
            matched_rule = evaluate_rules(rules,
                DATABASE=db,
                DATASET=obj['dataset'],
                POINT=obj['name']
            )
            if matched_rule is None:
              raise Exception("Datapoint {} from database {} does not match any rules!", obj['name'], db)
            obj['method'] = matched_rule['method']
            
            _units = {**units, **matched_rule['units']} if 'units' in matched_rule else units
            obj['refval'] /= _units['INP_DATASETEVAL'] # Units coming IN -> Divide
            
            if len(obj['deps']) != len(obj['coeffs']):
                raise Exception("Malformed row in dataset eval file: " + str(test))
            return obj
        return [*map(evalRow, filter(lambda r: len(r) > 1, evalreader))]

def read_deviations(db, out_dir=OUTPUTS_DIR, csv_name='IndValues.csv'):
    """
        Parse an `IndValues.csv` file into a map of (dft, dataset) -> list of datapoint deviations from reference (DP
        minus reference). Units are not transformed; Unit values from `IndValues.csv` are preserved. Returns a tuple of
        a list of DFT names, a list of dataset names, and the above map.
        * `db` - Database name
        * `out_dir` - Output directory to search in. Defaults to `OUTPUTS_DIR`
        * `csv_name` - The name of the value CSV. Defaults to `'IndValues.csv'`
    """
    path = os.path.join(out_dir, db, csv_name)
    results_list = []
    with open(path, newline='') as csvfile:
        evalreader = csv.reader(csvfile, dialect='excel')
        
        dft_names = evalreader.__next__()[3:]
        
        dft_db_map = {}
        def point(dft, ds_name, point):
            if not (dft, ds_name) in dft_db_map: dft_db_map[(dft, ds_name)] = []
            dft_db_map[(dft, ds_name)].append(point)
        
        datasets = set()
        for row in evalreader:
            # dp_name = row[0] # datapoint name; we don't track it
            ds_name = row[1]
            ref = float(row[2])
            datasets.add(ds_name)
            
            i = 3
            for dft in dft_names:
                point(dft, ds_name, abs(float(row[i]) - ref))
                i += 1
        
        return (dft_names, [*datasets], dft_db_map)

def get_dep_set(points):
    """
        Returns a set of dependency calculations for a list of points returned by `read_database_eval`. Each point is a
        tuple of (method, molecule geometry)
    """
    s = set()
    for p in points:
        for d in p['deps']:
            s.add((p['method'], d))
    return s

def get_regex_result(regexp, *path):
    "Find the last instance of `regexp` in the `mol.out` file at `*path`."
    with open(os.path.join(*path, 'mol.out'), 'r') as qcengine_out:
        # we are only interested in the last occurrence
        # https://stackoverflow.com/a/54277955/7853604
        match = None
        for m in re.finditer(regexp, qcengine_out.read()): match = m.groupdict()
        return match

def get_full_energy(point, dft, energy_out_dir=ENERGIES_DIR):
    """
        Calculate the complete energy of `point` based on the provided output files. All of the molecular geometries
        that make up `point` are considered and are summed with their leading coeffficients.
    """
    engine_exp = config["REGEXP"]
    energy_exp = '(?P<energy>[-+]?\d+\.\d+)'
    total_exp = re.compile(engine_exp + energy_exp)
    
    s = 0
    for i in range(len(point['deps'])):
        re_res = get_regex_result(total_exp, energy_out_dir, dft, point['method'], point['deps'][i])
        s += point['coeffs'][i] * float(re_res['energy']) / units['INP_QCENGINE'] # Divide since units are coming in
    return s

def get_runtime_seconds(point, dft, energy_out_dir=ENERGIES_DIR):
    """
        Calculate the runtime of `point` based on the provided output files in seconds. The runtimes of each component
        molecule are considered.
    """
    total_exp = re.compile(config["TIME_REGEXP"])
    
    s = 0
    for i in range(len(point['deps'])):
        re_res = get_regex_result(total_exp, energy_out_dir, dft, point['method'], point['deps'][i])
        # Find total in seconds based on all available named groups; Note that they're all optional, but one needs to
        # be specified to get a non-zero answer of course
        days =     float(re_res['days'])    if 'days'    in re_res else 0
        hours =   (float(re_res['hours'])   if 'hours'   in re_res else 0) + 24 * days
        minutes = (float(re_res['minutes']) if 'minutes' in re_res else 0) + 60 * hours
        seconds = (float(re_res['seconds']) if 'seconds' in re_res else 0) + 60 * minutes
        seconds += (float(re_res['millis']) if 'millis'  in re_res else 0) / 1000.0
        s += seconds
    return s

def gen_output_rows(points, headers, handler):
    """
        Given a list of input `headers`, generate entries in a dict. Each entry is an array of mapped `points`, which
        are mapped by calling `handler` with (point, header) as arguments. Useful for, say, making a spreadsheet with
        several parameters.
    """
    gen_column = lambda header: [handler(point, header) for point in points]
    
    return dict(zip(headers, map(gen_column, headers)))

def write_csv(file, **kwargs):
    """
        Write a CSV file with headers matching the kwargs, and items matching the iterable arguments.
    """
    with open(file, 'w', newline='') as output_file:
        csv_data = csv.writer(output_file, dialect='excel')
        
        csv_data.writerow([*kwargs.keys()])
        for key in kwargs.keys(): kwargs[key] = kwargs[key].__iter__()
        
        n_keys = len(kwargs.keys())
        n_stop = 0
        while True:
            row = []
            for arg in kwargs.values():
                try: row.append(arg.__next__())
                except StopIteration: n_stop += 1
            if n_stop == n_keys: return
            if n_stop > 0: raise RuntimeError('Arguments to writeCsv have different lengths!')
            csv_data.writerow(row)

# -------------------------------------------------------------------------------------------------------------------- #
#                                                       RULES                                                          #
# -------------------------------------------------------------------------------------------------------------------- #

# No need to run simple calculations on a cluster; This would take longer than running locally
localrules: ALL, IND_VALUES, RUN_TIMES, QCENGINE_INPUT

rule ALL:
    input:
        expand('{out_DIR}/{db}/IndValues.csv', db=DATABASES.keys(), out_DIR=OUTPUTS_DIR),
        expand('{out_DIR}/{db}/RunTimes.csv', db=DATABASES.keys(), out_DIR=OUTPUTS_DIR),
        expand('{out_DIR}/{db}/MADValues.csv', db=DATABASES.keys(), out_DIR=OUTPUTS_DIR),
        expand('{out_DIR}/{db}/RMSDValues.csv', db=DATABASES.keys(), out_DIR=OUTPUTS_DIR),
        expand('{out_DIR}/{db}/MDValues.csv', db=DATABASES.keys(), out_DIR=OUTPUTS_DIR),
        expand('{out_DIR}/{db}/MAXValues.csv', db=DATABASES.keys(), out_DIR=OUTPUTS_DIR),
        expand('{out_DIR}/{db}/MINValues.csv', db=DATABASES.keys(), out_DIR=OUTPUTS_DIR),
        expand('{out_DIR}/{db}/AMAXValues.csv', db=DATABASES.keys(), out_DIR=OUTPUTS_DIR)

# All of the output files that need to be build for a complete database
all_mol_outs = unpack(
    lambda wildcards: expand(
        '{e_DIR}/{dft}/{m[0]}/{m[1]}/mol.out',
        m=get_dep_set(read_database_eval(wildcards.db)),
        dft=DATABASES[wildcards.db]['dfts'],
        e_DIR=ENERGIES_DIR
    )
)

stat_files = {
    'MAD_VALUES': ('MADValues.csv', meanabs),
    'RMSD_VALUES': ('RMSDValues.csv', rms),
    'MD_VALUES': ('MDValues.csv', mean),
    'MAX_VALUES': ('MAXValues.csv', max),
    'MIN_VALUES': ('MINValues.csv', min),
    'AMAX_VALUES': ('AMAXValues.csv', maxabs),
}
for rname in stat_files.keys():
    rule:
        name: rname
        output: OUTPUTS_DIR + '/{db}/' + stat_files[rname][0]
        input: OUTPUTS_DIR + '/{db}/IndValues.csv'
        run:
            DB=DATABASES[wildcards.db]
            (dfts, datasets, dev_map) = read_deviations(wildcards.db)
            # Turn our map to lists of deviations into a map to single statistical measurements (floats)
            # Use the function provided to aggregate the list; See stat_files
            # Access by `rule` here (NOT `rname`) to work around Snakemake bug
            dev_map = {k: stat_files[rule][1](v) for k, v in dev_map.items()}
            
            write_csv(
                output[0],
                Functional = dfts,
                **{ds: [dev_map[(dft, ds)] for dft in dfts] for ds in datasets}
            )

rule IND_VALUES:
    # Generate the `IndValues.csv` file in the given output directory for the given database.
    # `IndValues.csv` should replicate the results in `Database/{db}/IndValues.csv` to within ~10^(-6) Hartree at least.
    output: OUTPUTS_DIR + '/{db}/IndValues.csv'
    input: all_mol_outs
    run:
        DB=DATABASES[wildcards.db]
        eval_points = read_database_eval(wildcards.db)
        write_csv(
            output[0],
            RefNames        = map(lambda p: p["name"], eval_points),
            DatasetRefNames = map(lambda p: p["dataset"], eval_points),
            RefValues       = map(lambda p: p["refval"] * units["OUT_INDVALUES"], eval_points), # Units going OUT -> Multiply
            **gen_output_rows(eval_points, DB["dfts"], scalef(get_full_energy, units["OUT_INDVALUES"]))
        )
rule RUN_TIMES:
    # Generate the `RunTimes.csv` file in the given output directory for the given database.
    output: OUTPUTS_DIR + '/{db}/RunTimes.csv'
    input: all_mol_outs
    run:
        DB=DATABASES[wildcards.db]
        eval_points = read_database_eval(wildcards.db)
        write_csv(
            output[0],
            RefNames        = map(lambda p: p["name"], eval_points),
            DatasetRefNames = map(lambda p: p["dataset"], eval_points),
            **gen_output_rows(eval_points, DB["dfts"], get_runtime_seconds)
        )

rule QCENGINE_RUN:
    # Invoke the calculation program specified in your config file
    input:      '{out_DIR}/{dft}/{method}/{molecule}/mol.in'
    output:     '{out_DIR}/{dft}/{method}/{molecule}/mol.out'
    # Config for running remotely (ex, in a Slurm cluster)
    resources:
        cpus_per_task=4,
        mem="16G",
        runtime=240,
        nodes=1,
        tasks=1,
    run:
        shell(config["QCENGINE_CALL"])

rule QCENGINE_INPUT:
    input:      GEOMETRIES_DIR + '/{molecule}.xyz'
    output:     '{out_DIR}/{dft}/{method}/{molecule}/mol.in'
    run:
        # Make parent dirs if necessary
        shell('mkdir -p {out_DIR}/{dft}/{method}/{molecule}/'.format(**wildcards))
        
        # Build the input file based on the molecular geometry and template file
        for file_name in output:
            molecule = wildcards.molecule
            dft = wildcards.dft
            method = METHODS[wildcards.method]
            template = TEMPL1
            with open(GEOMETRIES_DIR + '/{molecule}.xyz'.format(**wildcards), 'r') as f:
                f.readline()                   # skip first line
                molecule_data = f.read()[:-1]  # skip last NL
            with open(file_name, 'w') as f:
                f.write(open(template).read().format(
                    dft=dft,
                    molecule=molecule,
                    molecule_data=molecule_data,
                    nproc=nproc,
                    params=method,
                ))

