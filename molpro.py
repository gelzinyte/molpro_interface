import os
import sys
import re
import subprocess
import shutil

import numpy as np
from lxml import etree as et
from pathlib import Path
from collections import OrderedDict
from collections import abc
import xml.dom.minidom as minidom

import ase
from ase.units import Hartree, Bohr
from ase import Atoms

from ase.calculators.calculator import FileIOCalculator
from ase.calculators.calculator import Calculator
from ase.calculators.calculator import ReadError
from ase.calculators.calculator import SCFError
from ase.calculators.calculator import Parameters
from ase.calculators.calculator import all_changes

# TODO add a function to use the calculator parallelly?
# TODO figure out when at.info or at.arrays entries are too long for Molpro
#      and deal with it
# TODO add option for geometry optimisation with Molpro
#


__all__ = ['Molpro']



class Molpro(FileIOCalculator):
    r"""
Molpro Interface Documentation

The module defines an interface to MOLPRO 2012 (?).
Based on Python 2 molpro driver for quippy by Alan Nichol, James Kermode
and others (who?). Adapted to ASE by Elena Gelžinytė.

=========================   ==================================================
Calculator arguments        Description
=========================   ==================================================
``directory``               Relative path where all Molpro-related outputs will
                            be placed. Will be created if doesn't exist.


``label``                   The prefix for .inp, .xyz, .out, .xml, etc files.

``molpro_command``          Command to run molpro. Can also be set via the bash
                            environment variable ``MOLPRO_COMMAND`. If none is
                            given or found, will default to ``molpro``.

``scratch_dir``             Scratch directory for Molpro.

``restart``                 TODO

``ignore_bad_restart_file`` TODO


=========================   ==================================================

=========================   ==================================================
Molpro arguments            Description
=========================   ==================================================

``task``                    Either ``energy_only`` or ``gradient`` for calcu-
                            lating energies only or graidnets as well with
                            a Molpro call. Mandatory.

``command``                 Electronic structure method to execute. Also used
                            to extract the appropriate value from the output
                            file. Mandatory.

``basis``                   Basis for the method. Mandatory.

``functional``              Functional, if appropriate.

``memory``                  Memory in 'amount, unit' e.g. '300, w'

``command_block``           Molpro command block in the form of
                            '{COMMAND, options \\n directives \\n data \\n}'

``maxit``                   Maximum number of SCF iterations. Molpro default
                            is 60.

``template_path``           A path to template to be used instead of suplying
                            keyword arguments. Must have a ``geomtyp=xyz``
                            and an unasigned ``geom=`` which is linked to the
                            relevant file by the calculator. ``command`` must
                            also be specified to indicate which energy to pick
                            from the Molpro output file.


=========================   ==================================================


=========================   ==================================================
Internal Setting            Description
=========================   ==================================================

``_prepare_input_only``     (``=False``) If set to ``True``, calculator will
                            stop before calling Molpro.

``_run_molpro``             (``=True``) Whether to call Molpro. Will try to
                            read output even if Molpro is not called.

``_overwrite_old_outputs``  (``=True``) Removes old output files if found
                            before calling Molpro. Otherwise old outputs
                            are appended with a number.

``discard_results_on_any_change``
                            (''=True'') Resets the calculator if any of te
                            parameters have changed.

``_rename_existing_dir``   (``=True``) (TODO not implemented) when using a new instance
                           of the calculator, this will move directories out of
                           the way that would be overwritten otherwise,
                           appending a date string.

=========================   ==================================================


End Molpro Interface Documentation
"""



    implemented_properties = ['energy', 'forces']

    default_parameters = dict(task='gradient',
                              command='rks',
                              basis='6-31G*',
                              functional='b3lyp')

    discard_results_on_any_change = True

    supported_commands = ['CCSD(T)-F12', 'CCSD(T)', 'MP2', 'DF-MP2',
            'DF-RMP2', 'RKS', 'UKS', 'RHF', 'DF-RHF', 'HF', 'DF-HF']

    supported_tasks = ['energy_only', 'gradient']



    def __init__(self, restart=None, ignore_bad_restart_file=False, label='molpro',
               atoms=None, molpro_command=None, directory='MOLPRO',
                scratch_dir=None, **kwargs):

        self.__name__ = 'Molpro'
        self.molpro_command = get_molpro_command(molpro_command)


        # TODO sort out how command=molpro_command works out in FileIOCalculator
        FileIOCalculator.__init__(self, restart=restart,
                ignore_bad_restart_file=ignore_bad_restart_file, label=label,
                atoms=atoms, command=molpro_command, directory=directory, **kwargs)


        self.scratch_dir = scratch_dir

        self._prepare_input_only = False
        self._run_molpro = True
        self._overwrite_old_outputs = True
        # self._rename_existing_dir
        # TODO add option to clean up all the files after a run

        # Check that necesary calculator keyword arguments are present and take required values
        p = self.parameters
        if 'task' not in p.keys() or 'command' not in p.keys() or\
                    'basis' not in p.keys():
            raise RuntimeError('Need to specify a task, command and basis at the least ')

        if p.task not in self.supported_tasks:
            raise RuntimeError(f'{p.task} not supported, "task" must be one of {self.supported_tasks}')

        if p.command.upper() not in self.supported_commands:
            raise RuntimeError(f'{p.command} not supported, "command" must be one of {self.supported_commands}')


        if self.scratch_dir is not None:
            # TODO deal with existing dir and create a renamed one
            if not os.path.exists(self.scratch_dir):
               os.makedirs(self.scratch_dir)



    def calculate(self, atoms=None, properties=['energy', 'forces'], system_changes=all_changes):

        # TODO pass properties depending of what task is chosen
        # TODO figure out how properties are passed around
        Calculator.calculate(self, atoms, properties, system_changes) # sets self.atoms to atoms and creates directory
        self.write_input(atoms, properties, system_changes)
        if not self._prepare_input_only:
            if self._run_molpro:
                self.run()
            self.read() # pass properties too?
            # self.push_old_state()

    def write_input(self, atoms, properties=None, system_changes=None):

        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        p = self.parameters
        p.write(self.label + '.ase')
        ase.io.write(self.label + '.xyz', atoms)  # TODO check why not self.atoms


        if 'template_path' not in p.keys():
            # template with empty 'geom=' was not given, so writing from parameters
            with open(self.input_fname, 'w') as f:
                if 'memory' in p.keys():
                    # assume p.memory is given as a string 'amount, units'
                    f.write(f'memory, {p.memory}\n')
                f.write(f'geomtyp=xyz\n')
                f.write(f'geom={self.label}.xyz\n')
                f.write(f'basis={p.basis}\n')
                if 'command_block' in p.keys():
                    if p.command not in p.command_block:
                        raise RuntimeError(f'set command and command in command block should mach')  # use p.command to read stuff out of output.
                    f.write(f'{p.command_block}')
                    if p.command_block[:-2] != '\n':
                        f.write('\n')
                else:
                    f.write(f'{p.command}, {p.functional}')
                    if 'maxit' in p.keys():
                        f.write(f', maxit={p.maxit}')
                    f.write(f'\n')
                if p.task == 'gradient':
                    f.write('force\n')
        else:
            # template was given, so have to add the xyz file
            with open(p.template_path, 'r') as f:
                new_template = ''
                for line in f:
                    if 'geom=' in line:
                        new_line = f'geom={self.label}.xyz\n'
                        new_template += new_line
                    else:
                        new_template += line
                if p.task == 'gradient' and 'force' not in new_template:
                    new_template += 'force\n'

            with open(self.input_fname, 'w') as f:
                f.write(new_template)



    def run(self):
        ''' Excecutes Molpro and prints any stdout and stderr '''

        if self._overwrite_old_outputs:
            if os.path.exists(self.output_out_fname):
                os.remove(self.output_out_fname)
            if os.path.exists(self.output_xml_fname):
                os.remove(self.output_xml_fname)


        line_to_execute = f'{self.molpro_command} {self.input_fname} -o {self.output_out_fname}'
        if self.scratch_dir is not None:
            line_to_execute += f' -d {self.scratch_dir}'
        stdout, stderr = shell_stdouterr(line_to_execute)

        if stdout:
            print(f'Molpro call stdout:\n{stdout}')

        if stderr:
            print(f'Molpro call stderr:\n{stderr}')



    def read(self):
        """Read atoms, parameters and calculated properties from output file."""

        # maybe introduce function to read open file or just generic file as with castep.
        # Deal with reading from .out not only .xml, also from a fileobject. Or
        # is that not necesary?

        # check if .xml output is there
        if not os.path.isfile(self.label + '.xml'):
            raise ReadError('xml output file does not exist')

        # check for errors in the output file
        catch_molpro_errors(self.label+'.xml', self.parameters.command)

        # read parameters
        self.parameters = Parameters.read(self.label + '.ase')


        # read atoms from .xml file
        energy_from = self.parameters.command
        if self.parameters.task == 'gradient':
            extract_forces = True
        elif self.parameters.task == 'energy_only':
            extract_forces = False
        else:
            extract_forces = None

        self.atoms = read_xml_output(self.label+'.xml', energy_from=energy_from, extract_forces=extract_forces)

        # read calculated properties
        self.results['energy'] = self.atoms.info['energy']
        if self.parameters.task == 'gradient':
            self.results['forces'] = self.atoms.arrays['forces']





    @property
    def name(self):
        return self.__name__

    @property
    def input_fname(self):
        return self.label + '.inp'

    @property
    def output_out_fname(self):
        return self.label + '.out'
    
    @property
    def output_xml_fname(self):
        return self.label + '.xml'


class MolproDatafile(OrderedDict):
    """Class to wrap a molpro datafile"""

    def __init__(self, datafile=None, xml=None, atoms=None):
        OrderedDict.__init__(self)
        if datafile is not None:
            self.read(datafile)
        elif xml is not None:
            self.read_xml(xml)
        elif atoms is not None:
            self.update_from_atoms(atoms)

    def copy(self):
        new = MolproDatafile()
        new.update(self)  # this is inherited form OrderedDict
        return new

    def parse_line(self, line, key=None):

        # split a line by the first non-alphanumeric character, then cats that character
        # to start of second item in list, if there is one.
        if key is not None:
            # If key is present, indicates a commmand block - use semicolons
            # to split up lines containing directives and such
            # TODO Is there another use of the 'key' argument?
            #      It doesn't look like it...
            if not self[key]:
                self[key].append(';' + line)
            else:
                self[key][0] += (';' + line)
        else:
            nonalpha = re.compile(r'[^a-zA-Z0-9()\-\}\{]')  # '\-' works fine
            separator = re.findall(nonalpha, line)

            if len(separator) > 0:
                separator = separator[0]
            else:
                separator = ","
            # Note the colon ':' has special meaning for defining labels
            if separator == ':' and line.endswith(separator):
                # in which case the colon should be kept
                fields = [line, ]
            else:
                fields = line.split(separator, 1)
            key = fields[0].upper()
            if key not in self.keys():
                self[key] = []
                if len(fields) > 1 and fields[1] != '':
                    self[key].append(separator + fields[1])
            else:
                if key == "BASIS":
                    print('WARNING: CHANGING BASIS DOES UNPREDICTABLE THINGS')
                n = 2
                testkey = key
                while testkey in self.keys():  # find a key like charge#3 if charge and charge#2 are already taken
                    testkey = key + "#%s" % n
                    n += 1
                key = testkey
                self[key] = []
                if len(fields) > 1:
                    self[key].append(separator + fields[1])

    def read(self, datafile):
        if isinstance(datafile, abc.Mapping):
            self.update(datafile)
            return

        # If given a file first make a list containing the lines
        if type(datafile) == type(''):
            data = open(datafile, 'r')
            datalines = data.readlines()
        elif type(datafile) == type([]):
            datalines = datafile

        sub_lines = []
        current_key = None
        for compound_line in datalines:
            # molpro allows multiple statements on single line separated by ';' , so split these up
            lines = compound_line.split(';')
            for sub_line in lines:
                sub_lines.append(sub_line.strip())

        for line in sub_lines:
            # Skip comments and blank lines, NB molpro files may begin with  the line'***'
            if line.startswith('#') or line.startswith('!') or line.startswith('*') or line == '':
                continue


            self.parse_line(line, current_key)


    def read_xml(self, xml_output):
        """Read the input file from molpro output. Input should be filename, file-like object or list of lines"""
        if type(xml_output) == type(''):
            xml_output = open(xml_output, 'r')
            xml_output = xml_output.readlines()
        elif hasattr(xml_output, 'read'):
            xml_output = xml_output.readlines()

        # Remove newlines from end of each line in molpro_output
        xml_output = [line.strip() for line in xml_output]

        # Find the echo of the datafile in the molpro output
        try:
            datafile_start = xml_output.index('--><job>')
        except ValueError:
            raise ValueError('Unable to find echo of datafile in molpro output')

        datafile_lines = []
        i = datafile_start + 2  # skip over the comment indicator produced by molpro

        end_of_input = re.compile("Variables initialized")
        datafile_ended = False
        while datafile_ended == False:
            line = xml_output[i]
            i = i + 1
            if re.search("Variables initialized", line):
                datafile_ended = True
            elif line.strip() == '':
                continue  # skip arbitrary blank lines
            else:
                datafile_lines.append(line)

        self.read(datafile_lines)

    def write(self, datafile=sys.stdout):
        "Write molpro input file. datafile can be a filename or an open file"
        "At the moment is not used, since I'm constructing the template manually."

        if type(datafile) == type(''):
            datafile = open(datafile, 'w')

        for key, value in self.items():
            # iteritems important here because order of lines matters
            # if have multiple instances of a command, say, 'hf' or 'charge'
            # the n occurrences of that keyword after the first will have #n appended
            # e.g. hf, hf#2, hf#3, etc.
            if re.search('#', key):
                shortkey = key.split('#')[0]
            else:
                shortkey = key
            if len(value) > 1:
                # TODO this appears to be designed for geometry specifications
                #      but it also writes procedures when command blocks are
                #      intended - rethink this.
                datafile.write(shortkey + '={\n')
                for line in value:
                    datafile.write(line + '\n')
                datafile.write('}\n')
            elif value and value[0].startswith(';'):
                cmdblock_lines = value[0].split(';')
                datafile.write('{' + shortkey)
                for line in cmdblock_lines:
                    datafile.write(line + '\n')
                datafile.write('}\n')
            elif value:
                datafile.write('%s%s\n' % (shortkey, value[0]))
            else:
                datafile.write(shortkey + '\n')



def read_xml_output( xmlfile, energy_from=None, extract_forces=False,
                    extract_dipole=False, datafile=None, atoms=None):
    '''Returns atoms with energies and forces attached'''
    # TODO: maybe move to ase.io.molpro.py?

    if datafile is None:
        datafile = MolproDatafile(xml=xmlfile)
        if 'FORCE' in datafile:
            extract_forces = True

    energy_names = OrderedDict()
    energy_names['CCSD(T)-F12'] = ["total energy"]
    energy_names['CCSD(T)'] = ["total energy"]
    energy_names['MP2'] = ["total energy"]
    energy_names['DF-MP2'] = ["total energy"]
    energy_names['DF-RMP2'] = ["energy"]
    energy_names['RKS'] = ["Energy"]
    energy_names['UKS'] = ["Energy"]
    energy_names['RHF'] = ["Energy"]
    energy_names['DF-RHF'] = ["Energy"]
    energy_names['HF'] = ["Energy"]
    energy_names['DF-HF'] = ["Energy"]
    # etc

    gradient_names = OrderedDict()
    gradient_names['CCSD(T)'] = [""]
    gradient_names['RKS'] = ['RKS GRADIENT']
    gradient_names['UKS'] = ['UKS GRADIENT']
    gradient_names['MP2'] = ['MP2 GRADIENT']

    all_methods = OrderedDict()
    all_methods['HF'] = ["RHF"]
    all_methods['DF-HF'] = ["RHF"]
    all_methods['RHF'] = ["RHF"]
    all_methods['DF-RHF'] = ["RHF"]
    all_methods['MP2'] = ["MP2"]
    all_methods['DF-MP2'] = ["MP2"]
    all_methods['DF-RMP2'] = ["DF-RMP2"]
    all_methods['RKS'] = ["RKS"]
    all_methods['UKS'] = ["UKS"]
    all_methods['CCSD(T)-F12'] = ["CCSD(T)-F12a", "CCSD(T)-F12b"]
    all_methods['CCSD(T)'] = ["CCSD(T)"]

    if energy_from is None:
        raise RuntimeError(
            "don't know which energy to extract, use keyword "
            "energy_from with options " + str(
                [all_methods[k] for k in iter(all_methods)]).replace('[',
                                                                     '').replace(
                ']', ''))
    else:
        energy_from = energy_from.upper()

    # loop through datafile to look for methods.
    calcs = []  # holds the keys for getting correct method,
                # energy_name, gradient_name
    dfile_keys_stripped = [key.replace('{', '') for key in
                           datafile.keys()]
    data_keys_upper = [key.upper() for key in dfile_keys_stripped]


    for key in all_methods.keys():
        if key in data_keys_upper:
            calcs.append(key)
    dom = minidom.parse(xmlfile)

    # get atomic positions from the datafile
    elements = []
    position_matrix = []
    cml = dom.documentElement.getElementsByTagName('cml:atomArray')
    for l in cml[0].childNodes:
        if l.nodeType == 1:
            element = l.attributes['elementType'].value
            elements.append(element)
            posx = l.attributes['x3'].value.encode('ascii', 'ignore')
            posy = l.attributes['y3'].value.encode('ascii', 'ignore')
            posz = l.attributes['z3'].value.encode('ascii', 'ignore')
            position_matrix.append(
                [float(posx), float(posy), float(posz)])

    # create atoms object
    if atoms is None:
        position_matrix = np.array(position_matrix)
        atoms = Atoms(elements, positions=position_matrix)

    # extract from xml energy values for each of the methods found in datafile

    energy_found = False
    props = dom.documentElement.getElementsByTagName('property')
    for prop in props:
        prop_name = prop.attributes['name'].value
        prop_method = prop.attributes['method'].value
        for calc in calcs:
            if prop_name in energy_names[calc] and prop_method in \
                    all_methods[calc]:
                energy_param_name = "_".join([prop_method, prop_name])
                energy_param_name = energy_param_name.replace(" ", "_")
                energy_param = prop.attributes['value'].value
                my_energy = energy_param_name
                i_en = 1
                while my_energy in atoms.arrays.keys():
                    i_en += 1
                    my_energy = '_'.join([energy_param_name, str(i_en)])
                atoms.info[my_energy] = float(energy_param) * Hartree
                if prop_method == energy_from:
                    atoms.info['energy'] = float(energy_param) * Hartree
                    energy_found = True
            elif extract_dipole and prop_name == 'Dipole moment':
                dipole_param_name = "_".join([prop_method, prop_name])
                dipole_param_name = dipole_param_name.replace(" ", "_")
                log.info("found dipole moment: " + dipole_param_name)
                dipole_param = prop.attributes['value']
                atoms.arrays[dipole_param_name] = dipole_param

    if not energy_found:
        raise RuntimeError(
            f"couldn't find energy from {energy_from} prop method : "
            f"{prop_method}")

    # read gradients if requested
    if extract_forces:
        if not 'forces' in atoms.arrays.keys():
            atoms.arrays['forces'] = np.zeros((1,
                                                 3))  # suspicious -
            # only 1x3 matrix? also 1d or 2d array?

        grads = dom.documentElement.getElementsByTagName('gradient')
        force_matrix = grads[0].childNodes[0].data.split('\n')
        force_matrix = [str(i).split() for i in force_matrix]
        for _ in force_matrix:
            try:
                force_matrix.remove([])
            except ValueError:
                break
        force_matrix = [[(-1.0 * Hartree / Bohr) * float(j) for j in i]
                        for i in force_matrix]

        atoms.arrays['forces'] = np.array(force_matrix)

        if len(grads) != 1:
            for k in range(1, len(grads)):
                my_force = 'forces%s' % str(k + 1)
                force_matrix = grads[k].childNodes[0].data.split('\n')
                force_matrix = [str(i).split() for i in force_matrix]
                for i in force_matrix:
                    try:
                        force_matrix.remove([])
                    except ValueError:
                        break
                force_matrix = [
                    [(-1.0 * Hartree / Bohr) * float(j) for j in i]
                    for i in force_matrix]
                atoms.arrays[my_force] = np.array(force_matrix)

    return atoms



def catch_molpro_errors(filename, command):

    check_SCF_maxing_out(filename, command)
    catch_errors_and_warnings(filename)

def catch_errors_and_warnings(filename):
    '''checks for any "?Error", "No convergece", etc in the molpro output file'''
    pass

def check_SCF_maxing_out(filename, command, maxit=60):
    '''TODO Actually could just look for "No Convergence" in the output'''
    # and read other stuff too
    root = et.parse(filename).getroot()
    job = root.find('{http://www.molpro.net/schema/molpro-output}job')
    jobsteps = job.findall(
        '{http://www.molpro.net/schema/molpro-output}jobstep')
    for jobstep in jobsteps:
        if command.upper() in jobstep.attrib['command']:
            for child in jobstep:
                text = child.text
                if text is not None and 'ITERATION' in text:
                    text_lines = text.splitlines()
                    for idx, line in enumerate(text_lines):
                        if 'MAX. NUMBER OF ITERATIONS' in line:
                            maxit = int(re.search(r'\d+', line).group())
                        if 'Final occupancy' in line:
                            last_scf_line = text_lines[idx - 2].strip()
                            match = re.search(r'^\d+',
                                              last_scf_line).group()
                            last_iteration = int(match)
                            if last_iteration == maxit:
                                raise SCFError(
                                    'Maxed out of SCF iterations')
                            elif last_iteration > maxit:
                                raise RuntimeError(
                                    'last_iteration > self._scf_maxit, '
                                    'something is wrong with the code')
                        if '?APPARENTLY NO CONVERGENCE, EXIT AFTER THREE FURTHER ITERATIONS' in line:
                            raise SCFError('No convergence in SCF')



def get_molpro_command(molpro_command=''):

    if molpro_command:
        return molpro_command
    elif 'MOLPRO_COMMAND' in os.environ:
        return os.environ['MOLPRO_COMMAND']
    else:
        return 'molpro'


def shell_stdouterr(raw_command, cwd=None):
    """Abstracts the standard call of the commandline, when
    we are only interested in the stdout and stderr
    """
    stdout, stderr = subprocess.Popen(raw_command,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE,
                                      universal_newlines=True,
                                      shell=True, cwd=cwd).communicate()
    return stdout.strip(), stderr.strip()





