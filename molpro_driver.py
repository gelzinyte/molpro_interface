#!/usr/bin/env python
# molpro driver - Alan Nichol
# Adapted from castep_driver by James Kermode

# This script can be used to interface QUIP to MOLPRO. Within QUIP, you
# should create a FilePot potential like this:
# 
# type(Potential) :: pot
# ...
# call initialise(pot, "FilePot command=/path/to/molpro_driver.py")
#
# The variable MOLPRO_TEMPLATE also needs to point to a file containing the
# molpro input file, with geom= as a placeholder (you don't need to specify the file,
# because QUIP will pass the correct file along and the driver will update the input file
# accordingly.
# For now, if you want forces to be calculated, the FORCE keyword has to be in your molpro input
# and in this driver you need to set extract_forces=True.
# Running eval with this driver looks like:
#   ./eval at_file=geom.xyz E F init_args={FilePot command=$QUIP_ROOT/quippy/scripts/molpro_driver.py}
# and the $PYTHONPATH envvar needs to be set as well.

# Support for CCSD(T)-F12 calculations needs work.

# arguments:
#                       paths should be either singly-quoted or {}'ed
#   template,
#   append_lines        extra lines to be appended to template; if needed
#   molpro              Command used to execute molpro, with a %s where seed name should go;
#   test_mode           If set to True, don't actually run MOLPRO
#   working_dir         Working directory for MOLPRO. Set this to a local scratch directory if network file performance is poor.
#   energy_from         defined in molpro.py: 'CCSD(T)-F12', 'CCSD(T)', 'MP2', 'DF-MP2', DF-RMP2', 'RKS', 'RHF', 'DF-RHF', 'HF', 'DF-HF'
#


# ----------------------------------------------------------------
# Dependancies 
# ----------------------------------------------------------------

# Python standard library modules
import sys, string, os, os.path, shutil, glob, operator, xml.dom.minidom, logging, subprocess

import numpy as np

# TODO set normally once done
import molpro
from molpro import MolproDatafile

from ase import Atoms
from ase.io import read, write

from util import parse_params

# Set up logging
log = logging.getLogger('molpro_driver')
format = logging.Formatter('%(name)s - %(levelname)-10s - %(asctime)s - %(message)s')
handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.INFO)
handler.setFormatter(format)
log.addHandler(handler)
log.propagate = False
log.level = logging.INFO


def die(message):
    "Print error message and abort"
    log.critical(message)
    os.chdir(orig_dir)
    sys.exit(1)


class ParamError(Exception):
    pass


# Save starting directory
orig_dir = os.getcwd()

if len(sys.argv) < 3:
    die('%s usage: <xyzfile> <outputfile> [KEY=VALUE]...' % sys.argv[0])

xyzfile = sys.argv[1]
outfile = sys.argv[2]
geom = "geom_plain.xyz"
log.info("output to %s", outfile)

args_str = ''
if len(sys.argv) > 3:
    args_str = ' '.join(sys.argv[3:])

calc_args_str = parse_params(args_str)  # sits in util module, turns key=val pairs into dict
# what if I made it so that you passed append_lines="hf;ccsd(t)-f12;angstrom etc etc separated by ';'"
# then calc_args_str['append_lines'].split(';') makes a list of lines to be written.

print("Using calc args: {}".format(calc_args_str))

stem = os.path.basename(xyzfile)
stem_split = os.path.splitext(stem)
if stem_split[1] == '.xyz':  # Remove extension
    stem = stem_split[0]
logfile = stem + ".log.xyz"

# remove old input file, if it's there
old_input_fname = os.path.join(stem, stem)
if os.path.isfile(old_input_fname):
    os.remove(old_input_fname)

# ----------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------


# Template used for input file
# need to allow option for this to be blank so long as append_lines was read
if 'MOLPRO_TEMPLATE' in os.environ:
    MOLPRO_TEMPLATE = os.environ['MOLPRO_TEMPLATE']
elif 'template' in calc_args_str:
    MOLPRO_TEMPLATE = calc_args_str['template']
    del calc_args_str['template']
else:
    MOLPRO_TEMPLATE = orig_dir + '/template'

# extra lines (if any) to be appended to the template file
lines_to_append = []
if 'append_lines' in calc_args_str:
    log.info(calc_args_str['append_lines'])
    lines_to_append = calc_args_str['append_lines'].split(' ')
    del calc_args_str['append_lines']
    log.info(str(lines_to_append))
# Command used to execute molpro, with a %s where seed name should go
if 'MOLPRO' in os.environ:
    MOLPRO = os.environ['MOLPRO']
elif 'molpro' in calc_args_str:
    MOLPRO = calc_args_str['molpro']
    del calc_args_str['molpro']
else:
    MOLPRO = '~/molpro %s'

# If there's no %s, put seed name at end of string
if MOLPRO.find('%s') == -1:
    MOLPRO = MOLPRO + ' %s'

# If set to True, don't actually run MOLPRO
TEST_MODE = False
if 'test_mode' in calc_args_str:
    TEST_MODE = calc_args_str['test_mode']
    del calc_args_str['test_mode']

# Working directory for MOLPRO. Set this to a local scratch
# directory if network file performance is poor.
WORKING_DIR = '.'
if 'working_dir' in calc_args_str:
    WORKING_DIR = calc_args_str['working_dir']
    del calc_args_str['working_dir']

ENERGY_FROM = None
if 'energy_from' in calc_args_str:
    ENERGY_FROM = calc_args_str['energy_from']
    log.info("the energy of the frame is from " + ENERGY_FROM)
    del calc_args_str['energy_from']

BATCH_READ = False
BATCH_QUEUE = False

# forces will only be calculated by molpro if you request them
# extract_forces just determines whether they are passed to Atoms object

# ----------------------------------------------------------------

if os.path.splitext(MOLPRO_TEMPLATE)[1] != '.xml':
    # Read template input file
    try:
        datafile = molpro.MolproDatafile(MOLPRO_TEMPLATE)
        # datafile.write()

    except IOError:
        die("Can't open input file %s" % MOLPRO_TEMPLATE)
    except ValueError as message:
        die(str(message))

    # need to add XML handling here

# Read extended XYZ input file containing cluster
cluster = read(xyzfile)  # just first conformation

# remove old output file, if it's there
if os.path.exists(outfile):
    os.remove(outfile)

work_dir = os.path.join(WORKING_DIR, stem)

# Make working directory if necessary
if not os.path.isdir(work_dir):
    os.mkdir(work_dir)
os.chdir(work_dir)
#shutil.copyfile(os.path.join(orig_dir, MOLPRO_TEMPLATE), os.path.join(orig_dir, work_dir, MOLPRO_TEMPLATE))
#file_from = os.path.join(orig_dir, 'methane.inp')
#file_to = work_dir
#log.info("file_from: %s" % file_from)
#log.info("file_to: %s" % file_to)
#shutil.copyfile(os.path.join(orig_dir, 'methane'), os.path.join(work_dir, 'methane'))

if not BATCH_READ:
    # Load up old cluster, if it's there
    old_cluster_fname = stem + '_old.xyz'
    if os.path.exists(old_cluster_fname):
        log.info('found old cluster in file %s' % old_cluster_fname)
        try:
            # only reads the first conformation
            old_cluster = read(old_cluster_fname)
        except IOError:
            die('error opening old cluster file %s' % old_cluster_fname)

        if (len(cluster) == len(old_cluster)):
            # log.info('RMS position difference is %.3f A' % rms_diff(cluster.positions, old_cluster.postions))
            log.info('RMS positions difference is !!rms_diff() not implemented!! A')
        else:
            log.warning('number mismatch with previous cluster')

    else:
        old_cluster = cluster

    # Read template into MolproDatafile object
    datafile = MolproDatafile(datafile=MOLPRO_TEMPLATE)

    # Update the reference to the geometry file, or create one if it's not there
    if 'GEOMETRY' not in datafile.keys() and 'GEOM' not in datafile.keys():
        temp = MolproDatafile()
        if 'MEMORY' in datafile.keys():
            temp['MEMORY'] = datafile['MEMORY']
        temp['GEOM'] = []
        temp['GEOM'].append('=' + geom)
        for key in datafile.keys():
            if key != 'MEMORY':
                temp[key] = datafile[key]
        datafile = temp.copy()
    else:
        try:
            datafile['GEOM'].append('=' + geom)
        except:
            datafile['GEOMETRY'].append('=' + geom)

    # Append given lines to datafile
    if len(lines_to_append) > 0:
        for line in lines_to_append:
            datafile.parse_line(line)

    # And write to the molpro input file
    # datafile.write(datafile=stem)
    # write an ordinary xyz file (i.e. just species and positions)
    # cluster.write(dest=geom, properties=['species', 'pos'])
    write(geom, cluster, 'xyz')  # not sure this is completely ok

    # check if FORCE keyword present
    extract_forces = False
    if 'FORCE' in datafile.keys():
        extract_forces = True

# now invoke MolPro
if not BATCH_READ and not BATCH_QUEUE:
    if not molpro.run_molpro(datafile, MOLPRO, stem, test_mode=TEST_MODE):
        log.error('molpro run failed')

log.info("Molpro run has finished")
#shutil.copyfile(os.path.join(work_dir, stem+'.out'), os.path.join(orig_dir, stem+'.out'))
#shutil.copyfile(os.path.join(work_dir, geom), os.path.join(orig_dir, geom))
#shutil.copyfile(os.path.join(work_dir, stem), os.path.join(orig_dir, stem+'produced_input'))

# parse the XML output for energy, forces
cluster = molpro.read_xml_output(
    stem + '.xml', energy_from=ENERGY_FROM, extract_forces=extract_forces,
    datafile=datafile, cluster=cluster)

cluster.set_cell([20, 20, 20])

oldxyzfile = stem + '_old.xyz'
# Save cluster for comparison with next time
cluster.write(oldxyzfile, format='xyz')
# Also append it to the log file
logfile = open(logfile, 'a')
oldxyzfile = open(oldxyzfile, 'r')
logfile.writelines(oldxyzfile.readlines())

# Finally change back to original working directory and write output file
os.chdir(orig_dir)
cluster.write(outfile, format='extxyz')
