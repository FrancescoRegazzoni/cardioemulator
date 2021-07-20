
cardioemulator
==============

``cardioemulator`` is a Python library that allows to construct zero-dimensional emulators of the cardiac electromechanical function.
These emulators, surrogating the pressure-volume relationship of the cardiac chambers, are build through a data-driven automated algorithm on the basis of pressure-volume recordings.
Once constructed, these emulators allow to simulate the cardiac function at a very low computational cost.
For more details we refer to [#paper]_.

Documentation
-------------

.. toctree::
   :maxdepth: 1

   docs/installation
   Usage<_examples/example.ipynb>
   API reference <_autosummary/cardioemulator.html#http://>

.. autosummary::
   :toctree: _autosummary
   :template: custom-module-template.rst
   :recursive:

   cardioemulator


References
----------

.. [#paper] F. Regazzoni, A. Quarteroni "`Accelerating the convergence to a limit cycle in 3D cardiac electromechanical simulations through a data-driven 0D emulator <https://doi.org/10.1016/j.compbiomed.2021.104641>`_", *Computers in Biology and Medicine* 135 (2021) 104641

Citing
------

Please cite this library as: ::

   @article{regazzoni2021cardioemulator,
      author  = {Regazzoni, F. and Quarteroni, A.},
      title   = {Accelerating the convergence to a limit cycle in 3D cardiac
                 electromechanical simulations through a data-driven 0D emulator},
      journal = {Computers in Biology and Medicine},
      volume  = {135},
      pages   = {104641},
      year    = {2021},
      issn    = {0010-4825},
   }

Contacts
--------

Francesco Regazzoni, MOX - Politecnico di Milano (francesco.regazzoni@polimi.it)


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`