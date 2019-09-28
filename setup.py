from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy


modules = Extension(
    name="gimpact",
    sources=["gimpact.pyx"],
    language="c++",
    library_dirs = [r"D:\Software Libraries\GIMPACT\build\x64\Release"],
    libraries=["GIMPACT"],
    include_dirs=[r"D:\Software Libraries\GIMPACT\include", numpy.get_include()],
)
setup(
    name="gimpact",
    ext_modules=cythonize([modules])
)

#python setup.py build_ext --inplace
