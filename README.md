 ```
      [.           [..       [..   [.....    [.. [..      
     [. ..      [..   [.. [..   [..[..   [.. [.    [..    
    [.  [..    [..       [..       [..    [..[.     [..   
   [..   [..   [..       [..       [..    [..[... [.      
  [...... [..  [..       [..       [..    [..[.     [..   
 [..       [..  [..   [.. [..   [..[..   [.. [.      [.   
[..         [..   [....     [....  [.....    [.... [..    
 ```
 
# Full DFT Reference Calculation Tool

This tool is based on [ACCDB v1.0](https://github.com/peverati/ACCDB/tree/master); Morgante, P; Peverati, R.; J Comput
Chem, 2018, 40 (6), 839-848. [LINK](https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.25761).

The original ACCDB could do single-point calculations of every datapoint in a dataset, but...
1. You had to globally set all options other than the method and basis set -- Hence, manual configuration was required
to replicate the results of the study.
2. Only single point energies were calculated; A separate program was required to calculate the final output values
based on the component single-point energies in `DatasetEval.csv`.
3. No support files for compute clusters
4. Didn't extract the runtime
5. Units didn't match the source units of the study
6. Calculation of multiple DFTs differing only in dispersion is very inefficient as it requires a full SCF cycle

Our goal is to evaluate new DFTs (or, at least, newly implemented in the software) against a reference database. Hence,
we needed the full calculation to be completed in a reproducible manner.

See `config.example.yaml` for an example config file (comments are included to explain).

Our modifications to ACCDB redesign the architecture to work in the following steps:
1. First, you request the databases that you would like to build (stored in `Databases/`) as well as the DFTs (formerly
*methods*) that you would like to include by placing them in the relevant `DATABASES` parameter.
2. This software then reads the appropriate `DatasetEval.csv` file in the database to figure out which molecules need to
be built.
3. The `RULES` are used to pick out the appropriate `METHODS` (basis set & grid, currently) to build the datapoint for
the database. Each single point energy is stored in `Wavefunctions/<DFT name>/<method>/<datapoint name>/mol.out`.
  a. In addition to the energy (in `mol.out`), the wavefunction is exported to `mol.wf`
4. Similarly, calculations are then done using those wavefunction files to get the dispersion energy for each single
point and stored in `Dispersion/<DFT name>/<dispersion>/<method>/<datapoint name>/mol.out`
5. Each row in `DatasetEval.csv` represents a weighed sum of component single-point calculations; For example, an
interaction energy might be calculated by 1 times the energy of molecules together, minus the energy of each molecule on
its own. The software then computes the weighed sum of the individual datapoints, and places the output in
`Outputs/<DB name>/IndValues.csv` in a way that should reproduce the `IndValues.csv` in the `Databases/` directory.
6. Also, `RunTimes.csv` file is generated in a similar format with run times in seconds.

## Modifications

I'm working on a version to **make dispersion calculations more efficient/modular**. You're looking at it. Click
[here](https://github.com/UBCC3/ACCDB-fullcalculations/tree/master) to go back to the master branch.

## How to Run

The easiest way to use this software is probably to use our [Psi4 docker containers](https://github.com/UBCC3/psi4-docker-utils).
A `docker-compose` file is provided, so simply run:
```
docker-compose run snakemake_psi4
```
Using Docker helps make sure we have a consistent and reproducible calculation environment, in addition to making it
easy to swap the runtime version of Psi4 for development purposes. However, using Docker isn't always possible (such as
on university computing clusters). The most universal way to run this software is with Conda. At the time of writing,
that looks like:
```
conda create -n accdb
conda activate accdb
conda install -c conda-forge -c bioconda psi4 snakemake
snakemake
```

If you would like to make use of the parallel computing features (via Slurm), you could run:
```
conda create -n accdb
conda activate accdb
conda install -c conda-forge -c bioconda psi4 snakemake snakemake-executor-plugin-slurm
. fix_paths.sh
screen -LS sn snakemake -j 64 --executor slurm --default-resources slurm_account=<your-slurm-act> --rerun-incomplete
```
...which creates a GNU screen running Snakemake with the Slurm executor for 64 parallel calculations. The GNU screen
ensures that you can log out without killing Snakemake and come back in a few days. Simply resume the screen with
`screen -r sn` (note that I made `sn` the name of the screen). When you're in the screen, press `Ctrl+A D` to exit. A
logfile is also generated in `screenlog.0`. You can use `tail -f screenlog.0` to follow along without actually entering
the screen.

It is also necessary to **source fix_paths.sh** in order to work around read-only filesystems on many compute clusters;
Snakemake puts its cache and temp files in your home directory by default. `fix_paths.sh` fixes that by placing them in
`.tmp` in your local directory. Note that:
1. The `.tmp` file will be placed in whichever directory you run the script from
2. You **MUST** source the file (via `. fix_paths.sh` or `source fix_paths.sh`), not simply run it; Otherwise, the
environment variable changes are not captured.

**Note that you will probably need/want to customize the Slurm options a bit more, possibly by modifying the Snakefile itself!**

## Using other software
Following the modular approach of ACCDB, software other than Psi4 can be used. Simply create a `.tmpl` file for the
software input, then change the `QCENGINE_CALL` command and the two `REGEXP`s. Examples of other software usage has not
been provided, sorry.

## Units, units, units...
A large problem that we encountered was the different units between databases, software, and even within the same
database. For example, `MGCDB84` places reference values in `DatasetEval.csv` in Hartrees, but final values in
`IndValues.csv` are stored in kcal/mol instead. These issues are mitigated by the `UNITS` field in the config file --
While the software works in Hartrees internally, every unit going in and out of this software is convered according to
the scaling factors set in `UNITS` below, which can be obtained from the first *row* of [this table](https://ryutok.github.io/EnergyConversionTable/).
```yaml
# A method for dealing with the different units being thrown around everywhere
# Energy of 1 Hartree in the designated unit. For Psi4, this is 1.0 since it outputs in Hartrees.
# If your software outputs in eV, for example, this would be 27.211396641308 Hartrees/eV
# Find this in the first ROW of https://ryutok.github.io/EnergyConversionTable/
UNITS:
  # Energy scaling factor for the output of your program
  INP_QCENGINE: 1.0 # Psi4 1.0 Hartree/Hartree
  # Energy scaling factor for the study's reference values stored in DatasetEval.csv
  INP_DATASETEVAL: null # Set per-database above.
  
  # Energy scaling factor for what YOU want IndValues.csv to be in
  # Be aware that this must be set correctly to directly reproduce the values of studies!
  OUT_INDVALUES: 627.509 # Units switch once again: 627.509 Hartrees/kcal/mol
```

The value for `INP_DATASETEVAL` can also be overridden in the `RULES` field like so:
```yaml
RULES:
  ...
  - match:
      DATABASE: MGCDB84
    method: def2-QZVPPD_99x590
    units: { INP_DATASETEVAL: 1.0 }
  ...
```
This allows you to set the units properly per-study.

## The current version of ACCDB is v1.0 and it includes the following databases:
- **MGCDB84**: Head-Gordon's database, version 2017 [4,986 data in 91 subsets].
- **GMTKN55**: Grimme's main group thermochemistry, kinetics and noncovalent interactions database [1,499 data in 55 subsets] https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/GMTKN/gmtkn55
- **Minnesota Database** (version 2015B, excluding geometries and solid state databases): Truhlar's database [471 data in 34 subsets] https://comp.chem.umn.edu/db2015/
- **DP284**: Head-Gordon's 2018 dipole and polarizabilities databases [284 data in 2 subsets].
- **Metals&EE**: Collection of databases containing transition metals [210 data in 8 subsets].
- **W4-17**: Martin's Weizman 2017 database [1,042 data in 5 subsets].

PLEASE READ CAREFULLY THE "**HOW TO CITE**" SECTION BELOW WHEN CITING ANY OF THE DATABASES COLLECTED HERE!

## How to cite:
This work is published under GNU-GPL license and credit must be given to the authors of the original material, as appropriate. Most of the databases are collected from free sources on the internet, and we act only as a central repository for them. (If you are the author or the owner of one of such databases and want your data to be excluded from this repository, please do contact us, and we will remove it immediately.)
Some databases and the automation workflows have been created by us, and are presented here for the first time. Such databases and software are subject to the same license and are freely available to use and modify, as long as proper credit is given to the ACCDB project (see below).
Everything in ACCDB is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 

**Citations for published scientific material:**

*For pre-existing databases or subsets:* 

You **must** cite the appropriate original reference for each database, as indicated within each folder, and you can cite the ACCDB project in addition, if you find this useful. 

*For the following new features presented here for the first time:* 

If you are using MN-RE or W4-17-RE you should cite the ACCDB project (see below), **and** the original paper by R.Bartlett et. al., as indicated in the README of each database.  If you are using the Porphyrins subset you should cite the ACCDB project (see below), **and** the original papers by K. Pierloot et. al. where the benchmarks where obtained, as indicated in the README of each database.  If you are using any of the automation scripts you should cite the ACCDB project (see below), **and** the snakemake original article.

*The reference for the ACCDB project is*: 

Morgante, P; Peverati, R.; J Comput Chem, 2018, 40 (6), 839-848. [LINK](https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.25761)

