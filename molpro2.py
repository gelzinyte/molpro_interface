"""The module defines an interface to MOLPRO 2012 (?).
    Based on Python 2 molpro driver for quippy by Alan Nichol, James Kermode
    and others (????). Adapted to ASE by Elena Gelžinytė.
    Loosely based on castep ase calculator and Tamas' Orca (based on NWChem?).
"""

# TODO list:
#   Raise appropriate calculator errors
#   Clean up imports to follow python conventions
#   can I leave scratch directory to its own devices if not specified?

import os
import pdb
import ase
from lxml import etree as et
from ase.calculators.calculator import FileIOCalculator, \
                                       CalculatorSetupError, Calculator, \
                                       ReadError, PropertyNotImplementedError,\
                                       SCFError, PropertyNotPresent

from ase.calculators.calculator import Parameters
from ase.calculators.calculator import all_changes
from pathlib import Path


__all__ = ['Molpro']



class Molpro(FileIOCalculator):
    r"""
Molpro Interface Documentation


Introduction
============

Something something


Getting Started
===============

something something 
something something


Running the Calculator
======================


Arguments:
==========

=========================  ===================================================
Keyword                    Description
=========================  ===================================================
``directory``              Relative path where all Molpro-related items will
                           be placed. Will be created if doesn't exist.
                           TODO: add behaviour with how to deal with existing
                           directory.

``label``                  The prefix for .inp, .xyz, .out, .xml, etc files.

``molpro_command``         Command to run molpro. Can also be set via the bash
                           environment variable ``MOLPRO_COMMAND`. If none is
                           given or found, will default to ``molpro``.

``check_molpro_version``   Should implement???

``keyword_tolerance``      Should implement?

=========================  ===================================================


Arguments:
==========

=========================  ===================================================
Keyword                    Description
=========================  ===================================================


=========================  ===================================================


Additional Settings
===================

=========================  ===================================================
Internal Setting           Description
=========================  ===================================================
``_molpro_command``        (``=molpro``): the actual shell command used to
                           call MOLPRO

``_prepare_input_only``     (``=False``): If set to True, the calculator will
                           create \*.inp and \*.xyz files but not start the
                           calculation itself.

``_rename_existing_dir``   (``=True``) : when using a new instance
                           of the calculator, this will move directories out of
                           the way that would be overwritten otherwise,
                           appending a date string.

``_set_atoms``             (``=False``) : setting this to True will overwrite
                           any atoms object previously attached to the
                           calculator when reading the output file.  By de-
                           fault, the read() function will only create a new
                           atoms object if none has been attached and other-
                           wise try to assign forces etc. based on the atom's
                           positions.  ``_set_atoms=True`` could be necessary
                           if one uses internal geometry optimization
                           (``blah``)
                           because then the positions get out of sync.
                           *Warning*: this option is generally not recommended
                           unless one knows one really needs it. There should
                           never be any need, if MOLPRO is used as a
                           single-point calculator.

=========================  ===================================================

Special features
================


``.keyword.clear()``
  ??maybe implement?

``.initialize()``
  Creates all needed input in the ``_directory``. This could then be copied to
  and ran in a place without ASE or Python.

``print(calc)``
  Prints a short summary of the calculator settings and atoms.


Notes/Issues
============

blah blah


End Molpro Interface Documentation
"""

    # keyword arguments to update: 'memory', check write_input

    implemented_properties = ['energy', 'forces']

    default_parameters = dict(task='gradient',
                              command='rks',
                              basis='6-31G*',
                              functional='b3lyp')

    discard_results_on_any_change = True # resets calculator if any of the parameters have changed

    # from what was in read_xml_output. Where should assert that the command is in it?
    supported_commands = ['CCSD(T)-F12', 'CCSD(T)', 'MP2', 'DF-MP2',
            'DF-RMP2', 'RKS', 'UKS', 'RHF', 'DF-RHF', 'HF', 'DF-HF']



    def __init__(self, restart=None, ignore_bad_restart_file=False, label='MOLPRO/molpro',
               atoms=None, molpro_command=None, **kwargs):

        # TODO add an option to read out a template file to read out all the parameters
        # TODO function to make a template file out of everything
        # TODO should I input and link .xyz file or coordinates directly?
        # TODO deal with the directories - scratch (what's the default?), working
        # dir to write in and out all the tmp files, an option to not delete that
        # and how to deal with duplicate files that molpro just moves over
        # Does creating and removing directories with every calculation slow it down?
        # Maybe add a cleanup function that once everything is done removes all the tmp dirs
        # TODO add a function to use the calculator parallely.
        # TODO have a problem if .xyz file that is passed in has a number of at.info stuff.
        # Set it up so that the molpro input has no info, but when it's read out, at.info from original at file is returned.
        # actually probably taken care of by ase itself.
        # TODO when to make functions part of Molpro class and when just include them in this .py file?

        self.__name__ = 'Molpro'
        self.molpro_command = get_molpro_command(molpro_command)

        # TODO sort out how command=molpro_command works out in FileIOCalculator
        FileIOCalculator.__init__(self, restart=restart,
                ignore_bad_restart_file=ignore_bad_restart_file, label=label,
                atoms=atoms, command=molpro_command, **kwargs)


        # TODO How to deal with (internal) parameters like these that have pre-set values?
        self._no_scf_iterations = None
        self.prepare_input_only = False
        self._rename_existing_dir = True
        self._scf_maxit = 60

        # TODO should be set by the fileio calculator?
        self.set(**kwargs)



    def calculate(self, atoms=None, properties=['energy'], systm_changes=all_changes):

        # TODO pass properties depending of what task is chosen
        Calculator.calculate(self, atoms, properties, system_changes) # sets self.atoms to atoms and creates directory
        self.write_input(self, atoms, properties, system_changes)
        if not self.prepare_input_only:
            self.run()
            self.read() # pass properties too?
            self.push_old_state()

    def write_input(self, atoms, properties=None, system_changes=None):
        # TODO make so that either read from .xml, .out, .inp or something else.
        # creates a directory
        # TODO FileIOCalculator.write_input creates a directory, so does Calculator.calculate. Why duplicate? or could write input without calcultion?
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        p = self.parameters
        p.write(self.label + '.ase')
        ase.io.write(self.label + '.xyz', atoms)  # TODO check why not self.atoms

        # Could/should make a list of accepted commands and if template is given
        # just read the command out of there. That's assuming only one command is present though
        if 'command' not in p.keys():
            raise CalculatorSetupError('must give the command to run with molpro')

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
                        raise Exception(f'set command and command in command block should mach')  # use p.command to read stuff out of output.
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

        # TODO sort out scratch_dir
        line_to_execute = f'{self.molpro_command} {self.input_fname} -o {self.output_fname}'
        stdout, stderr = shell_stdouterr(line_to_execute)

        if stdout:
            print(f'Molpro call stdout:\n{stdout}')

        if stderr:
            print(f'Molpro call stderr:\n{stderr}')



    def read(self):
        """Read atoms, parameters and calculated properties from output file.
        TODO: currently everything is stripped in preparation to do MOLPRO calculation.
        Since MOLPRO doesn't take xyz comments. Need to deal with stripping atoms
         in preparation and somehow read back/preserve (?) """

        # maybe introduce function to read open file or just generic file as with castep.
        # deal with reading from .out not only .xml, also from a fileobject

        # check if .xml output is there
        if not os.path.isfile(self.label + '.xml'):
            raise ReadError('xml output file does not exist')
        elif os.path.isfile(self.label + '.out'):
            raise RuntimeError('.xml is not present and reading from .out files not implemented yet')

        # check for errors in the output file
        self.catch_molpro_errors(self.label+'.xml')

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

        # TODO deal with all the info entries from atoms object passed to the calculator

        atoms = self.read_xml_output(self.label+'.xml', energy_from=energy_from, extract_forces=extract_forces)
        self.atoms = atoms

        # read calculated properties
        self.results['energy'] = atoms.info['energy']
        if 'gradient' in self.parameters.keys():
            self.results['forces'] = atoms.arrays['forces']




    @property
    def name(self):
        return self.__name__

    @property
    def input_fname(self):
        return self.label + '.inp'

    @property
    def output_fname(self):
        return self._label + '.out'


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
            nonalpha = re.compile('[^a-zA-Z0-9()\-\}\{]')  # '\-' works fine
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
        position_matrix = np.array(position_matrix).T
        if not 'ANGSTROM' in datafile.keys() and not 'angstrom' in \
                                                     datafile.keys():
            position_matrix = position_matrix * (1.0 / 0.529177249)
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
    extract_errors_and_warnings(filename)

def extract_errors_and_warnings(filename):
    '''checks for any "?Error", etc in the molpro output file'''
    pass

def check_SCF_maxing_out(filename, command):
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
                            self._scf_maxit = maxit
                        if 'Final occupancy' in line:
                            last_scf_line = text_lines[idx - 2].strip()
                            # ^ and search are probably doubling up
                            match = re.search(r'^\d+',
                                              last_scf_line).group()
                            last_iteration = int(match)
                            if last_iteration == self._scf_maxit:
                                raise SCFError(
                                    'Maxed out of SCF iterations')
                            elif last_iteration > self._scf_maxit:
                                raise RuntimeError(
                                    'last_iteration > self._scf_maxit, '
                                    'something is wrong with the code')
                        if '?APPARENTLY NO CONVERGENCE, EXIT AFTER THREE FURTHER ITERATIONS' in line:
                            # TODO make a test for this one
                            raise SCFError('No convergence in SCF')

                        # TODO Spot messages for not converging and
                        #  stopping after a few iterations
                    # only interested in this (TODO: unless more than one command), so can return (?).
                    return



def get_molpro_command(molpro_command=''):
    # TODO documment this somewhere
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





