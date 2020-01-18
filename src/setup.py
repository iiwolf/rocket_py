from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize("rocket_driver_c.pyx")
)