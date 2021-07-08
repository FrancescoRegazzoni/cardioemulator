.. _installation:

============
Installation
============

The library requires `Python 3 <https://www.python.org/>`_.
The dependencies are listed in ``\requirements.txt``. If you are using ``pip``, you can install them by ::

    $ pip install -r requirements.txt

To use the library, you can choose among the following options.

Option 1: Install the library by `setuptools <https://setuptools.readthedocs.io/>`_
------------------------------------------------------------------------------------------

From the main folder of this repository run: ::

    $ pip install .


Option 2: Add the library path to Python paths
------------------------------------------------------------------------------------------

Linux / macOS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Add to your ``.bashrc``: ::

    export PYTHONPATH="${PYTHONPATH}:/path/to/cardioemulator"

Windows
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Add the path of ``cardioemulator`` to the ``PYTHONPATH`` environment variable.
Alternatively, you can write the path of ``cardioemulator`` inside a ``pth`` file, as described `here <https://docs.python.org/3/using/windows.html#finding-modules>`_.

Option 3: Add the library path within your python script
------------------------------------------------------------------------------------------

Put the following lines at the beginning of your Python script: ::

    import sys
    sys.path.append("/path/to/cardioemulator")
    import cardioemulator
