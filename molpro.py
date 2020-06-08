# TODO add stuff like __all__, if __name__=='main', etc
import sys, os, os.path, logging
import numpy as np
import ase
from ase.io import read, write
from ase.calculators.calculator import Calculator
from ase.io.extxyz import key_val_dict_to_str

from util import parse_params
import molpro_old
from molpro_old import MolproDatafile


class Molpro(Calculator):

    implemented_properties = ['energy', 'free_energy', 'forces']
    # TODO implement virial, (stress, local_virial, local_energy, stresses, energies?)


    def __init__(self, calc_args=None, work_dir='MOLPRO', calculation_always_required=False, **kwargs):
        # TODO what other properties do I need??

        self.__name__ = 'Molpro'

        # TODO set to only energy being default, as that only may be necessary sometimes
        self._default_properties = ['energy', 'forces']
        self.calculation_always_required = calculation_always_required # TODO why?, why here?

        # TODO quippy has more stuff in this initiate
        Calculator.__init__(self, **kwargs)

        # TODO initialise class attributes (?) to None (??)
        self._work_dir = work_dir
        self.atoms = None
        # if isinstance(calc_args, dict):
        #     calc_args = key_val_dict_to_str(calc_args)
        # elif calc_args is None:
        #     calc_args = ""
        # self.calc_args = calc_args
        if not isinstance(calc_args, dict):
            raise TypeError('Please pass calc_args as dictionary, for now')
        self.calc_args = calc_args



    def calculate(self, atoms=None, properties=None, system_changes=None,
                  forces=None, virial=None, local_energy=None,
                  local_virial=None, vol_per_atom=None,
                  copy_all_results=True, calc_args=None, add_arrays=None,
                  add_info=None, **kwargs):
        # TODO what's up with system_changes, etc

        # TODO do I actually need this
        # handling the property inputs
        if properties is None:
            properties = self.get_default_properties()
        else:
            properties = list(set(self.get_default_properties() + properties))

        if len(properties) == 0:
            raise RuntimeError('Nothing to calculate')

        for prop in properties:
            if prop not in self.implemented_properties:
                raise RuntimeError("Don't know how to calculate property '%s'" % prop)

        Calculator.calculate(self, atoms, properties, system_changes)

        # TODO should check if calculation is required
        # args_str = self.calc_args
        # if calc_args is not None:
        #     if isinstance(calc_args, dict):
        #         calc_args = key_val_dict_to_str(calc_args)
        #     args_str += ' ' + calc_args
        # if kwargs is not None:
        #     args_str += ' ' + key_val_dict_to_str(kwargs)

        calc_args = self.calc_args

        # over from molpro_driver.py
        # -----------------------------------------------------------------------------------------

        geom = 'molpro_input.xyz'
        write(geom, self.atoms, 'xyz')

        # Set up logging
        log = logging.getLogger('molpro_driver')
        format = logging.Formatter('%(name)s - %(levelname)-10s - %(asctime)s - %(message)s')
        handler = logging.StreamHandler(sys.stdout)
        handler.setLevel(logging.INFO)
        handler.setFormatter(format)
        log.addHandler(handler)
        log.propagate = False
        log.level = logging.INFO


        orig_dir = os.getcwd()

        # if len(sys.argv) < 3:
        #     die('%s usage: <xyzfile> <outputfile> [KEY=VALUE]...' % sys.argv[0], log, orig_dir)
        #
        # xyzfile = sys.argv[1]
        # outfile = sys.argv[2]
        # geom = "geom_plain.xyz"
        # log.info("output to %s", outfile)

        # args_str = ''
        # if len(sys.argv) > 3:
        #     args_str = ' '.join(sys.argv[3:])

        # calc_args = parse_params(args_str)  # sits in util module, turns key=val pairs into dict
        # what if I made it so that you passed append_lines="hf;ccsd(t)-f12;angstrom etc etc separated by ';'"
        # then calc_args['append_lines'].split(';') makes a list of lines to be written.

        # log.info("Using calc args: {}".format(calc_args))

        # stem = os.path.basename(xyzfile)
        # stem_split = os.path.splitext(stem)
        # if stem_split[1] == '.xyz':  # Remove extension
        #     stem = stem_split[0]
        # logfile = stem + ".log.xyz"
        #
        # # remove old input file, if it's there
        # old_input_fname = os.path.join(stem, stem)
        # if os.path.isfile(old_input_fname):
        #     os.remove(old_input_fname)

        ## for above: have calc_args passed via initiate or calculate or both.

        # ----------------------------------------------------------------
        # Parameters
        # ----------------------------------------------------------------

        # Template used for input file
        # need to allow option for this to be blank so long as append_lines was read
        if 'MOLPRO_TEMPLATE' in os.environ:
            MOLPRO_TEMPLATE = os.environ['MOLPRO_TEMPLATE']
        elif 'template' in calc_args.keys():
            MOLPRO_TEMPLATE = calc_args['template']
            del calc_args['template']
        else:
            MOLPRO_TEMPLATE = orig_dir + '/template'

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

        # Working directory for MOLPRO. Set this to a local scratch
        # directory if network file performance is poor.
        WORKING_DIR = '.'
        if 'working_dir' in calc_args.keys():
            WORKING_DIR = calc_args['working_dir']
            del calc_args['working_dir']

        ENERGY_FROM = None
        if 'energy_from' in calc_args.keys():
            ENERGY_FROM = calc_args['energy_from']
            log.info("the energy of the frame is from " + ENERGY_FROM)
            del calc_args['energy_from']

        BATCH_READ = False
        BATCH_QUEUE = False

        # forces will only be calculated by molpro if you request them
        # extract_forces just determines whether they are passed to Atoms object

        # ----------------------------------------------------------------

        if os.path.splitext(MOLPRO_TEMPLATE)[1] != '.xml':
            # Read template input file
            try:
                datafile = molpro_old.MolproDatafile(MOLPRO_TEMPLATE)
                # datafile.write()

            except IOError:
                die("Can't open input file %s" % MOLPRO_TEMPLATE, log, orig_dir)
            except ValueError as message:
                die(str(message))

        work_dir = self._work_dir
        # Make working directory if necessary
        if not os.path.isdir(work_dir):
            os.mkdir(work_dir)
        os.chdir(work_dir)

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

        # check if FORCE keyword present
        extract_forces = False
        if 'FORCE' in datafile.keys():
            extract_forces = True

        stem = 'molpro_out'

        if not BATCH_READ and not BATCH_QUEUE:
            if not molpro_old.run_molpro(datafile, MOLPRO, stem, test_mode=TEST_MODE):
                log.error('molpro run failed')

        self.atoms = molpro_old.read_xml_output(
            stem + '.xml', energy_from=ENERGY_FROM, extract_forces=extract_forces,
            datafile=datafile, cluster=self.atoms)

        os.chdir(orig_dir)

        self.results['energy'] = self.atoms.info['energy']
        if 'force' in properties:
            self.results['forces'] = np.copy(atoms.arrays['forces'])


    def get_default_properties(self):
        return self._default_properties[:]


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