# import os
# from setuptools import setup
# from setuptools.extension import Extension
# from Cython.Build import cythonize
# import numpy


# # Package meta-data.
# NAME = 'gimpact'
# DESCRIPTION = 'An unofficial python extension for the GImpact collision library.'
# URL = 'https://github.com/StephenNneji'
# AUTHOR = 'Stephen Nneji'
# VERSION = '0.1.0'

# path=os.path.abspath(os.path.dirname(__file__))
# try:
#     with open(os.path.join(path, 'README.md'), encoding='utf-8') as f:
#         long_description = '\n' + f.read()
# except FileNotFoundError:
#     long_description = DESCRIPTION


# modules = Extension(
#     name="gimpact",
#     sources=["gimpact.pyx", "GIMPACT/src/gimpact.cpp", "GIMPACT/src/gim_contact.cpp", "GIMPACT/src/gim_math.cpp", 
#              "GIMPACT/src/gim_boxpruning.cpp", "GIMPACT/src/gim_memory.cpp","GIMPACT/src/gim_tri_tri_overlap.cpp",
#              "GIMPACT/src/gim_trimesh.cpp","GIMPACT/src/gim_trimesh_capsule_collision.cpp", "GIMPACT/src/gim_trimesh_ray_collision.cpp",
#              "GIMPACT/src/gim_trimesh_sphere_collision.cpp", "GIMPACT/src/gim_trimesh_trimesh_collision.cpp"],
#     language="c++",
#     include_dirs=["GIMPACT/include", numpy.get_include()],
# )

# setup(name=NAME,
#       version=VERSION,
#       description=DESCRIPTION,
#       long_description=long_description,
#       long_description_content_type="text/markdown",
#       url=URL,
#       author=AUTHOR,
#       license='BSD',
#       classifiers=[
#           'Development Status :: 5 - Production/Stable',
#           'Intended Audience :: Developers',
#           'Intended Audience :: Science/Research',
#           'Topic :: Games/Entertainment :: Simulation',
#           'Topic :: Scientific/Engineering :: Physics',
#           'License :: OSI Approved :: BSD License',
#           'Operating System :: Microsoft :: Windows',
#           'Operating System :: POSIX',
#           'Programming Language :: Python :: 2.7',
#           'Programming Language :: Python :: 3',
#           'Programming Language :: Python :: 3.4',
#           'Programming Language :: Python :: 3.5',
#           'Programming Language :: Python :: 3.6',
#           'Programming Language :: Python :: 3.7',
#           'Programming Language :: Python :: Implementation :: CPython',
#           'Programming Language :: Cython'
#       ],
#       keywords='GImpact,Trimesh,Collision detection,Cython',
#       ext_modules=cythonize([modules]),
#       setup_requires = ["cython>=0.29.13", "numpy>=1.16.4"],
#       install_requires = ["numpy>=1.16.4"],
# )
from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy


modules = Extension(
    name="gimpact",
    sources=["gimpact.pyx", "GIMPACT/src/gimpact.cpp", "GIMPACT/src/gim_contact.cpp", "GIMPACT/src/gim_math.cpp", 
             "GIMPACT/src/gim_boxpruning.cpp", "GIMPACT/src/gim_memory.cpp","GIMPACT/src/gim_tri_tri_overlap.cpp",
             "GIMPACT/src/gim_trimesh.cpp","GIMPACT/src/gim_trimesh_capsule_collision.cpp", "GIMPACT/src/gim_trimesh_ray_collision.cpp",
             "GIMPACT/src/gim_trimesh_sphere_collision.cpp", "GIMPACT/src/gim_trimesh_trimesh_collision.cpp"],
    language="c++",
    include_dirs=["GIMPACT/include", numpy.get_include()],
)
setup(
    name="gimpact",
    ext_modules=cythonize([modules])
)

#python setup.py build_ext --inplace