"""
`cardioemulator` library.

.. rubric:: Classes

.. autosummary::
   :toctree: _autosummary
   :template: custom-class-template.rst
   :recursive:

   Emulator
   Emulator_parametric


.. rubric:: Functions

.. autosummary::
   :toctree: _autosummary
   :recursive:

   build_emulator
   build_emulator_parametric
   get_PV_relationship
   transform_ES_elastance
   transform_time_shift
   compute_convergence



"""

from ._emulator import Emulator, Emulator_parametric, get_PV_relationship
from ._build_emulator import build_emulator
from ._build_emulator_parametric import build_emulator_parametric
from ._transform_emulator import transform_ES_elastance, transform_time_shift
from ._compute_convergence import compute_convergence

from ._version import __version__