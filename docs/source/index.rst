.. idrfeatlib documentation master file, created by
   sphinx-quickstart on Mon Feb  3 18:52:17 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

idrfeatlib documentation
========================

Welcome to ``idrfeatlib``, the beginnings of a library for feature-based analysis
of protein sequences. This library is meant to be a highly extensible bioinformatics
framework for analyzing sequence features, finding novel sequences with similar
features as some target sequence, and so on.

This project is currently in a somewhat unstable in development phase right now,
so I would say use this code with the knowledge that even core functionalities
may be renamed or otherwise changed, at least until I think of a proper versioning
protocol.

This means: if you include this as a dynamically installed dependency (i.e. Google
colab `!git clone`), you should point it to a specific branch/version number and
not just the latest or main branch. 

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   usage
   native-features