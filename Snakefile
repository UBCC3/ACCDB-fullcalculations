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
     basis='[-\w]+',
     restriction='RHF|UHF',

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

# Energy of 1 Hartree in your output unit. For Psi4, this is 1.0 since it outputs in Hartrees.
# If your software outputs in eV, for example, this would be 27.211396641308 Hartrees/eV
# Find this in the first ROW of https://ryutok.github.io/EnergyConversionTable/
escale = config['ENERGY_SCALE']

hartree_to_kcal_mol = 627.509

# -------------------------------------------------------------------------------------------------------------------- #
#                                                  SUPPORT FUNCTIONS                                                   #
# -------------------------------------------------------------------------------------------------------------------- #

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
        * `name` - The name of the datapoint
        * `dataset` - The name of the dataset from which the point originates. Corresponds to `DatasetRefNames` in the
            `IndValues.csv` file for the particular dataset.
        * `refval` - The reference energy value in **HARTREES**
        * `deps` - A list of molecular geometries that make up this datapoint
        * `coeffs` - Coefficients for the energy value of each geometry in `deps`. Often, this is 1, -1, -1 for a diamer
            and two monomer geometries, respectively. These coefficients are used to calculate the final energy.
        * `method` - Based on the set of rules in `RULES`, a calculation method is determined.
    """
    path = os.path.join(db_dir, db, csv_name)
    results_list = []
    with open(path, newline='') as csvfile:
        evalreader = csv.reader(csvfile, dialect='excel')
        
        def evalRow(row):
            obj = {
              'name': row[0], # Name of the datapoint
              'dataset': '_'.join(row[0].split('_')[:-1]), # Name of the source dataset. TODO: This logic is a bit questionable
              'refval': float(row[-1]), # Reference energy in Hartrees
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
            
            if len(obj['deps']) != len(obj['coeffs']):
                raise Exception("Malformed row in dataset eval file: " + str(test))
            return obj
        return [*map(evalRow, filter(lambda r: len(r) > 1, evalreader))]

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
        s += point['coeffs'][i] * float(re_res['energy']) / escale
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

# -------------------------------------------------------------------------------------------------------------------- #
#                                                       RULES                                                          #
# -------------------------------------------------------------------------------------------------------------------- #

rule ALL:
    input:
        expand('{out_DIR}/{db}/IndValues.csv', db=DATABASES.keys(), out_DIR=OUTPUTS_DIR),
        expand('{out_DIR}/{db}/RunTimes.csv', db=DATABASES.keys(), out_DIR=OUTPUTS_DIR)

for result in ["IndValues.csv", "RunTimes.csv"]:
    rule:
# Generate the `IndValues.csv` and `RunTimes.csv` files in the given output directory for the given database.
# `IndValues.csv` should replicate the results in `Database/{db}/IndValues.csv` to within less than a percent.
        name: f"COMPILED_{result}"
        output: OUTPUTS_DIR + '/{db}/' + result
        input: unpack(lambda wildcards: expand('{e_DIR}/{dft}/{m[0]}/{m[1]}/mol.out', m=get_dep_set(read_database_eval(wildcards.db)), dft=DATABASES[wildcards.db]['dfts'], e_DIR=ENERGIES_DIR)), # If somebody knows how to get this to span more than one line (instead of being a very ugly one-liner without comments), please let me know!
        run:
            DB=DATABASES[wildcards.db]
            with open(output[0], 'w', newline='') as output_file:
                csv_data = csv.writer(output_file, dialect='excel')
                
                if result == 'IndValues.csv':
                    csv_data.writerow(['RefNames', 'DatasetRefNames', 'RefValues', *DB['dfts']])
                else:
                    csv_data.writerow(['RefNames', 'DatasetRefNames', *DB['dfts']])
                
                eval_points = read_database_eval(wildcards.db)
                for point in eval_points:
                    if result == 'IndValues.csv':
                        # Note that we write everything in kcal/mol because these studies keep switching units without
                        # ever actually LABELLING the units or specifying that a switch is happening. Oh well...
                        csv_data.writerow([
                            point["name"],
                            point["dataset"],
                            point["refval"] * hartree_to_kcal_mol,
                            *map(lambda dft: get_full_energy(point, dft) * hartree_to_kcal_mol, DB['dfts'])
                        ])
                    else:
                        csv_data.writerow([
                            point["name"],
                            point["dataset"],
                            *map(lambda dft: get_runtime_seconds(point, dft), DB['dfts'])
                        ])

rule QCENGINE_RUN:
#This rule calls the program.
    input:      '{out_DIR}/{dft}/{method}/{molecule}/mol.in'
    output:     '{out_DIR}/{dft}/{method}/{molecule}/mol.out'
    run:
        shell(config["QCENGINE_CALL"])

rule QCENGINE_INPUT:
#This rule writes the input file for the program of choice, based on the template given in the configfile
    input:      '{out_DIR}/{dft}/{method}/{molecule}/.keep'
    output:     '{out_DIR}/{dft}/{method}/{molecule}/mol.in'
    run:
        for file_name in output:
            molecule = wildcards.molecule
            dft = wildcards.dft
            method = METHODS[wildcards.method]
            template = TEMPL1
            with open('tolaunch/{molecule}.xyz'.format(molecule=wildcards.molecule), 'r') as f:
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

rule LAUNCH_DIR:
#This rule copies the xyz from the source folder to the 'tolaunch' directory.
      input:      GEOMETRIES_DIR + '/{molecule}.xyz'
      output:     'tolaunch/{molecule}.xyz'
      run:
                  shell('cp {input} {output}')

rule QCENGINE_DIRS:
#This rule builds the input and output files.
    input:      'tolaunch' + '/{molecule}.xyz'
    output:     '{out_DIR}/{dft}/{method}/{molecule}/.keep'
    shell:      'touch {output}'
