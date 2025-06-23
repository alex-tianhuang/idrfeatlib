Welcome to ``idrfeatlib``, the beginnings of a library for feature-based analysis
of protein sequences. This library is meant to be a highly extensible bioinformatics
framework for analyzing sequence features, finding novel sequences with similar
features as some target sequence, and so on.

Usage
=====

.. _installation:

Installation
------------

To use this ``idrfeatlib`` package,
pip install the idrfeatlib package from github:

.. code-block:: console

    (.venv) $ git clone https://github.com/alex-tianhuang/idrfeatlib
    (.venv) $ pip install ./idrfeatlib

You probably also should install `tqdm <https://pypi.org/project/tqdm/>`_,
as it is not a dependency of the libraries, but most of my scripts use it.
It is a progress bar library.

.. code-block:: console

    (.venv) $ pip install tqdm

If you have ``tqdm``, you can check that the installation worked by running
a basic IDR feature calculation script on some example fasta @ `test.fa`:

.. code-block:: console

    (.venv) $ python ./idrfeatlib/scripts/featurize.py test.fa /dev/stdout

.. _Computing Metric:

Computing Metric
----------------
Follow this section to compute the mean and weights associated with either
(a) An input fasta file of IDRs or
(b) An input fasta file of proteins + a csv of IDR regions

(a)
Suppose I have an input fasta file of IDRs @ `IDR_FILE.fasta`.
I will run the following:

.. code-block:: console

    (.venv) $ python scripts/featurize.py IDR_FILE.fasta OUTPUT_FEATURES.csv
    (.venv) $ python scripts/compute-metric.py OUTPUT_FEATURES.csv OUTPUT_ORIGIN_WEIGHTS.csv

(b)
Suppose I have an input fasta file of proteins @ `PROTEINS.fasta` + csv of IDR regions @ `IDRS.csv`.
The csv of IDR regions looks like:

.. code-block::

    ProteinID,RegionID,Start,Stop
    DDX3X,N-IDR,0,142
    protein_x,region3,24,72
    ...

I can then run:

.. code-block:: console

    (.venv) $ python scripts/featurize.py PROTEINS.fasta --input-regions IDRS.csv OUTPUT_FEATURES.csv
    (.venv) $ python scripts/compute-metric.py OUTPUT_FEATURES.csv OUTPUT_ORIGIN_WEIGHTS.csv --input-labels ProteinID RegionID

Using custom features
---------------------
Follow this section to add your own features.

Suppose I have a script @ `FOO.py`:

.. code-block:: python

    def bar(sequence, secret_param):
        return len(sequence) + secret_param

To add this feature to the available features, begin by copying the
:doc:`native features <native-features>` into a custom file `MY_FEATURES.json`.
It should look something like this:

.. code-block:: json

    {
        "features": {
            "AA_A": {
                "compute": "percent_residue",
                "residue": "A"
            },
            "AA_B": "...",
            "..."
        },
        "..."
    }

In the ``features`` section of that json, add this:

.. code-block:: json

    {
        "features": {
            "my_feature": {
                "compute": "custom",
                "libpath": "FOO.py",
                "funcname": "bar",
                "kwargs": {
                    "secret_param": 42
                }
            },
            "AA_A": {
                "compute": "percent_residue",
                "residue": "A"
            },
            "AA_B": "...",
            "..."
        },
        "..."
    }

Now, in subsequent analyses, add ``--feature-file MY_FEATURES.json`` to the end of every command.
For example, featurize a fasta file like so:

.. code-block:: console
    
    (.venv) $ python scripts/featurize.py IDR_FILE.fasta OUTPUT_FEATURES.csv --feature-file MY_FEATURES.json

Feature mimic design
--------------------
To design sequences, you'll need a weights file (see :ref:`Computing Metric`).

Suppose you have such a file @ `WEIGHTS.csv`, which roughly looks like:

.. code-block::

    Label,Feature1,Feature2,...
    origin,0,1,...
    weights,0.25,0.78,...

Use this to design on a protein fasta file + input regions csv like so:

.. code-block:: console
    
    (.venv) $ python scripts/feature-mimic.py PROTEINS.fa --input-regions REGIONS.csv WEIGHTS.csv DESIGN_OUTPUT.csv

Some helpful options are:

``--n-random``

    Design this many replicates starting with random sequences.
    e.g. if you want to design 10 sequences per region:

.. code-block:: console
    
    (.venv) $ python scripts/feature-mimic.py PROTEINS.fa --input-regions REGIONS.csv WEIGHTS.csv DESIGN_OUTPUT.csv --n-random 10

``--num-processes`` or ``-np``

    Use this many processes.
    If this number is greater than 1 (not default), you will need to install ``pathos`` and ``tqdm_pathos``.
    ``pathos`` is a multiprocessing library that is superior to the standard ``multiprocessing`` library
    (uses ``dill`` for pickling, and therefore doesn't fail all the time)
    ``tqdm_pathos`` makes a shared-process progress bar.

    e.g. to use 8 processes:

.. code-block:: console
    
    (.venv) $ python scripts/feature-mimic.py PROTEINS.fa --input-regions REGIONS.csv WEIGHTS.csv DESIGN_OUTPUT.csv -np 8

``--seeds-file``

    Design sequences starting with seeded random sequences.
    Seeds are provided in a `SEEDS.csv` file, which has columns "ProteinID" (and "RegionID" if ``--input-regions`` is specified)
    and a "Seed" column with integer seeds. Can be multiple seeds for one protein/region.
