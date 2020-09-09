# TODO add stuff like __all__, if __name__=='main', etc
import logging
import os
import os.path
import re
import sys
from collections import OrderedDict, abc
from copy import deepcopy

import numpy as np
import xml.dom.minidom as minidom
from ase import Atoms
from ase.calculators.calculator import Calculator
from ase.io import write
from ase.io.extxyz import key_val_dict_to_str
from ase.units import Hartree, Bohr


# import molpro_old
# from molpro_old import MolproDatafile

# TODO list
# * checkout ase calculator erros that might be appropriate

class Molpro(Calculator):
    """
    Arguments:

    Keyword                 Description
    ======================= ====================================================
    calc_args               dictionary of arguments transfered from
                            molpro_driver.py.

    directory               Relative path where all Molpro-related items will be
                            placed. Will be created if doesn't exist.
                            TODO: add option to move existing directory.

    label                   Prefix for .inp, .xyz, .out, .xml, etc files.

    calculation_always_required
                            TODO


    """

    # TODO implement virial, (stress, local_virial, local_energy, stresses, energies?)

    def __init__(self, calc_args=None, directory='MOLPRO', label='molpro', atoms=None, \
                 calculation_always_required=False, **kwargs):

        self.__name__ = 'Molpro'

        # TODO set to only energy being default, as that only may be necessary sometimes
        self._default_properties = ['energy', 'forces']
        self.calculation_always_required = calculation_always_required  # TODO why?, why here?

        Calculator.__init__(self, restart=None, ignore_bad_restart_file=False, label=None, atoms=atoms, **kwargs)

        # TODO initialise other class attributes (?) to None (??)
        self.implemented_properties = ['energy', 'free_energy', 'forces']
        self._directory = directory
        self._label = label

        if atoms is not None:
            atoms.set_calculator(self)

        if not isinstance(calc_args, dict):
            raise TypeError('Please pass calc_args as dictionary, for now')
        self.calc_args = deepcopy(calc_args)

    def calculate(self, atoms=None, properties=None, system_changes=None, copy_all_results=True):
        # TODO what's up with system_changes, etc

        # TODO do I actually need this
        # handling the property inputs
        if properties is None:
            properties = self.get_default_properties()
        else:
            #properties = list(set(self.get_default_properties() + properties))
            properties = list(set(['energy'] + properties))


        for prop in properties:
            if prop not in self.implemented_properties:
                raise RuntimeError("Don't know how to calculate property '%s'" % prop)

        Calculator.calculate(self, atoms, properties, system_changes)

        calc_args = deepcopy(self.calc_args)
        label = self._label

        # do something if atoms is None

        # over from molpro_driver.py
        # -----------------------------------------------------------------------------------------

        # Set up logging
        log = logging.getLogger('molpro_driver')
        format = logging.Formatter('%(name)s - %(levelname)-10s - %(asctime)s - %(message)s')
        handler = logging.StreamHandler(sys.stdout)
        handler.setLevel(logging.INFO)
        handler.setFormatter(format)
        log.addHandler(handler)
        log.propagate = False
        log.level = logging.INFO

        log.disabled = True  # For now; somehow messages pile up in the cycle
        log.info("Using calc args: {}".format(key_val_dict_to_str(calc_args)))

        orig_dir = os.getcwd()

        # ----------------------------------------------------------------
        # Parameters
        # ----------------------------------------------------------------

        # Template used for input file
        # need to allow option for this to be blank so long as append_lines was read
        MOLPRO_TEMPLATE = calc_args['template']
        del calc_args['template']

        # extra lines (if any) to be appended to the template file
        lines_to_append = []
        if 'append_lines' in calc_args.keys():
            log.info(calc_args['append_lines'])
            lines_to_append = calc_args['append_lines'].split(' ')
            del calc_args['append_lines']
            log.info(str(lines_to_append))
        # Command used to execute molpro, with a %s where seed name should go
        if 'MOLPRO' in os.environ:
            MOLPRO = os.environ['MOLPRO']
        elif 'molpro' in calc_args.keys():
            MOLPRO = calc_args['molpro']
            del calc_args['molpro']
        else:
            MOLPRO = '~/molpro %s'

        # If there's no %s, put seed name at end of string
        if MOLPRO.find('%s') == -1:
            MOLPRO = MOLPRO + ' %s'

        # If set to True, don't actually run MOLPRO
        TEST_MODE = False
        if 'test_mode' in calc_args.keys():
            TEST_MODE = calc_args['test_mode']
            del calc_args['test_mode']

        ENERGY_FROM = None
        if 'energy_from' in calc_args.keys():
            ENERGY_FROM = calc_args['energy_from']
            log.info("the energy of the frame is from " + ENERGY_FROM)
            del calc_args['energy_from']

        # forces will only be calculated by molpro if you request them
        # extract_forces just determines whether they are passed to Atoms object

        # ----------------------------------------------------------------

        if os.path.splitext(MOLPRO_TEMPLATE)[1] != '.xml':
            # Read template input file
            try:
                datafile = MolproDatafile(MOLPRO_TEMPLATE)
                # datafile.write()

            except IOError:
                die("Can't open input file %s" % MOLPRO_TEMPLATE, log, orig_dir)
            except ValueError as message:
                die(str(message))

        directory = self._directory
        # Make working directory if necessary
        if not os.path.isdir(directory):
            os.mkdir(directory)
        os.chdir(directory)

        geom = f'{label}_input.xyz'
        write(geom, self.atoms, 'xyz')

        # Read template into MolproDatafile object
        datafile = MolproDatafile(datafile=MOLPRO_TEMPLATE)

        if 'forces' in properties and 'FORCE' not in datafile.keys():
            lines_to_append +=['FORCE']

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

        # check if FORCE keyword present
        extract_forces = False
        if 'FORCE' in datafile.keys():
            extract_forces = True
        if extract_forces == True and 'forces' not in properties: 
            properties += ['forces']

        if not run_molpro(datafile, MOLPRO, label, test_mode=TEST_MODE):
            log.error('molpro run failed')

        self.atoms = read_xml_output(
            label + '.xml', energy_from=ENERGY_FROM, extract_forces=extract_forces,
            datafile=datafile, cluster=self.atoms)

        os.chdir(orig_dir)

        self.results['energy'] = self.atoms.info['energy']
        if 'forces' in properties:
            self.results['forces'] = np.copy(self.atoms.arrays['forces'])

        if isinstance(copy_all_results, bool) and copy_all_results:
            #print('copying results')
            self.atoms.info['energy'] = self.results['energy']
            if 'forces' in self.results.keys():
                self.atoms.info['forces'] = self.results['forces'] 
                #at.arrays['forces'] = self.results['forces'].copy()

    def get_default_properties(self):
        return self._default_properties[:]


def run_molpro(datafile, molpro, stem, test_mode=False):
    # Invoke molpro and return true if completed successfully

    log = logging.getLogger('molpro_driver')

    # write datafile
    datafile.write(stem)

    # check command line
    if not '%s' in molpro: molpro = molpro + ' %s'

    if test_mode:
        log.info('test mode: not running molpro')

    else:
        # Remove old output file
        try:
            os.remove(stem + '.out')
        except:
            pass

        # log.info('contents of cwdir')
        # print(os.listdir())
        # log.info('molpro command:')
        # log.info(molpro % stem)
        os.system(molpro % stem)
        # log.info('contents of cwdir')
        # print(os.listdir())

        #    error = subprocess()
    got_error = False

    return not got_error


def die(message, log, orig_dir):
    "Print error message and abort"
    log.critical(message)
    os.chdir(orig_dir)
    sys.exit(1)


def _check_arg(arg):
    """Checks if the argument is True bool or string meaning True"""

    true_strings = ('True', 'true', 'T', 't', '.true.', '.True.')

    if arg is None:
        return 'n'
    else:
        if isinstance(arg, bool):
            if arg:
                return 'y'
            else:
                return 'n'
        if isinstance(arg, str):
            if arg in true_strings:
                return 'y'
            else:
                return 'n'
        return 'add'


# ### over from old molpro.py

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
                if key == "BASIS":  # this warning should be passed through logging
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

            '''
            # Check if any brackets in line
            open_bracket = re.search('{', line)
            close_bracket = re.search('}', line)

            if open_bracket and close_bracket:
                # superfluous brackets, no actual multi-line command or data
                # TODO this is bad if the brackets actually delimit a command block...
                line = line.replace("{", "")
                line = line.replace("}", "")
                open_bracket = False
                close_bracket = False
                self.parse_line(line.strip(), current_key)

            # check if command block starting
            elif open_bracket:
                if current_key == None:
                    line = re.sub('{', "", line).strip()
                    if line != "":
                        self.parse_line(line, current_key)
                        current_key = list(self.keys())[-1]
                    else:
                        raise ValueError("Parse error in datafile: standalone open bracket")
                else:
                    raise ValueError("Parse error in datafile: nesting of curly {} brackets is illegal")

            # check if end of command block reached
            elif close_bracket:
                if current_key == None:
                    raise ValueError("Parse error in datafile:  check pairing of curly {} brackets")
                else:
                    line = re.sub('}', "", line).strip()
                    if line != "":
                        self.parse_line(line, current_key)
                        current_key = None
                    else:
                        current_key = None
                        continue
            # normal line - no brackets
            '''
            # else:
            #     self.parse_line(line, current_key)
            self.parse_line(line, current_key)

    def read_from_molpro_output(self, molpro_output):
        """Read the input file from molpro output. Input should be filename, file-like object or list of lines"""
        if type(molpro_output) == type(''):
            molpro_output = open(molpro_output, 'r')
            molpro_output = molpro_output.readlines()
        elif hasattr(molpro_output, 'read'):
            molpro_output = molpro_output.readlines()

        # Remove newlines from end of each line in molpro_output
        molpro_output = [line.strip() for line in molpro_output]

        # Find the echo of the datafile in the molpro output
        try:
            datafile_start = molpro_output.index('default implementation of scratch files=df')
        except ValueError:
            raise ValueError('Unable to find echo of datafile in molpro output')

        datafile_lines = []
        i = datafile_start + 2  # skip over the blank line produced by molpro

        # If geometry is specified as a filename molpro will insert a comment which we should remove
        include = re.compile("Including file")
        end_of_input = re.compile("Variables initialized")
        datafile_ended = False
        while datafile_ended == False:
            line = molpro_output[i]
            i = i + 1
            if re.search("Variables initialized", line):
                datafile_ended = True
            elif re.search("Including file", line):
                continue
            elif line.strip() == '':
                continue  # skip arbitrary blank lines
            else:
                datafile_lines.append(line)

        self.read(datafile_lines)

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


def read_xml_output(xmlfile, energy_from=None, extract_forces=False, extract_dipole=False, datafile=None, cluster=None):
    # parse an xml output file and return cluster with updated info
    # datafile tells which energies, forces to look for, cluster Atoms object which gets returned,
    # this is echoed in the xml file so can be left out
    # If extract_forces is not given and the FORCE keyword is found in datafile,
    # the default is to set extract_forces=True

    log = logging.getLogger('molpro_driver')

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
        log.critical("don't know which energy to extract, use keyword energy_from with options " + str(
            [all_methods[k] for k in iter(all_methods)]).replace('[', '').replace(']', ''))

    # loop through datafile to look for methods.
    calcs = []  # holds the keys for getting correct method, energy_name, gradient_name
    dfile_keys_stripped = [key.replace('{', '') for key in datafile.keys()]
    data_keys_upper = [key.upper() for key in dfile_keys_stripped]
    # data_keys_upper = [key.upper() for key in datafile.keys()]
    for key in all_methods.keys():
        if key in data_keys_upper:
            calcs.append(key)
    dom = minidom.parse(xmlfile)

    elements = []
    position_matrix = []
    cml = dom.documentElement.getElementsByTagName('cml:atomArray')

    for l in cml[0].childNodes:
        if l.nodeType == 1:
            # element = l.attributes['elementType'].value.encode('ascii','ignore')
            element = l.attributes['elementType'].value
            # elements.append(atomic_numbers[element])
            elements.append(element)
            posx = l.attributes['x3'].value.encode('ascii', 'ignore')
            posy = l.attributes['y3'].value.encode('ascii', 'ignore')
            posz = l.attributes['z3'].value.encode('ascii', 'ignore')
            position_matrix.append([float(posx), float(posy), float(posz)])
    # TODO fix this
    if cluster is None:
        # cluster = Atoms(n=len(elements))
        # cluster.set_atoms(elements)
        # position_matrix = np.array(position_matrix).T
        position_matrix = np.array(position_matrix)
        # if not 'ANGSTROM' in datafile.keys() and not 'angstrom' in datafile.keys():
        #     position_matrix = position_matrix * (1.0 / 0.529177249)
        # cluster.pos[:,:]=position_matrix
        # #note this leaves the lattice undefined

        cluster = Atoms(elements, positions=position_matrix)

    # now look for each of these energies in xml file
    energy_found = False
    props = dom.documentElement.getElementsByTagName('property')
    for prop in props:
        # prop_name = prop.attributes['name'].value.encode('ascii','ignore')
        # prop_method = prop.attributes['method'].value.encode('ascii','ignore')
        prop_name = prop.attributes['name'].value
        prop_method = prop.attributes['method'].value
        for calc in calcs:
            if prop_name in energy_names[calc] and prop_method in all_methods[calc]:
                energy_param_name = "_".join([prop_method, prop_name])
                energy_param_name = energy_param_name.replace(" ", "_")
                # log.info("found "+energy_param_name)
                # dated routines for finding monomer pairs, triplets in Topology module
                # energy_param=prop.attributes['value'].value.encode('ascii','ignore')
                energy_param = prop.attributes['value'].value
                my_energy = energy_param_name
                i_en = 1
                while my_energy in cluster.arrays.keys():
                    i_en += 1
                    my_energy = '_'.join([energy_param_name, str(i_en)])
                cluster.info[my_energy] = float(energy_param) * Hartree
                if prop_method == energy_from:
                    cluster.info['energy'] = float(energy_param) * Hartree
                    energy_found = True
            elif extract_dipole and prop_name == 'Dipole moment':
                dipole_param_name = "_".join([prop_method, prop_name])
                dipole_param_name = dipole_param_name.replace(" ", "_")
                log.info("found dipole moment: " + dipole_param_name)
                dipole_param = prop.attributes['value']
                cluster.arrays[dipole_param_name] = dipole_param

    if not energy_found:
        log.critical(f"couldn't find energy from {energy_from} prop method : {prop_method}")

    # read gradients if requested
    if extract_forces:
        if not 'forces' in cluster.arrays.keys():
            # cluster.add_property('force', 0.0, n_cols=3)
            cluster.arrays['forces'] = np.zeros((1, 3))  # suspicious - only 1x3 matrix? also 1d or 2d array?

        grads = dom.documentElement.getElementsByTagName('gradient')
        force_matrix = grads[0].childNodes[0].data.split('\n')
        force_matrix = [str(i).split() for i in force_matrix]
        for i in force_matrix:
            try:
                force_matrix.remove([])
            except ValueError:
                break
        force_matrix = [[(-1.0 * Hartree / Bohr) * float(j) for j in i]
                        for i in force_matrix]

        cluster.arrays['forces'] = np.array(force_matrix)
        # cluster.arrays['forces'] = cluster.arrays['force']

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
                force_matrix = [[(-1.0 * Hartree / Bohr) * float(j) for j in i]
                                for i in force_matrix]
                cluster.arrays[my_force] = np.array(force_matrix)

    return cluster
