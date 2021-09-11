import io
from setuptools import setup
from setuptools import find_packages

with open("requirements.txt", "r") as f:
    requirements = [line.strip() for line in f.readlines()]

exec(open('cardioemulator/_version.py').read())

setup(
    name = "cardioemulator",
    version = __version__,
    description = "Library implementing a 0D emulator of cardiac 3D electromechanical models.",
    author="Francesco Regazzoni",
    author_email = "francesco.regazzoni@polimi.it",
    license = "MIT",
    install_requires = requirements,
    packages = find_packages(),
    scripts = [
        'apps/cardioemulator_build',
        'apps/cardioemulator_build_parametric',
        'apps/cardioemulator_transform',
        ],
)