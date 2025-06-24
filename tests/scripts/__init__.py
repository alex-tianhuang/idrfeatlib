import os
import pathlib
import subprocess as sp

SCRIPTS_DIR = "{}/scripts/".format(pathlib.Path(os.path.dirname(os.path.realpath(__file__))).parent.parent)

assert os.path.exists(SCRIPTS_DIR)

def run_script(basename, argv):
    """Run a script from the `scripts` directory in this repo."""
    sp.run(["python", "{}/{}.py".format(SCRIPTS_DIR, basename)] + argv, check=True)