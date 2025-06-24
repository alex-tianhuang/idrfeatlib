import os
import pathlib
import subprocess as sp

SCRIPTS_DIR = "{}/scripts/".format(pathlib.Path(os.path.dirname(os.path.realpath(__file__))).parent.parent)

assert os.path.exists(SCRIPTS_DIR)

def run_script(basename, argv):
    sp.run(["python", "{}/{}".format(SCRIPTS_DIR, basename)] + argv, check=True)