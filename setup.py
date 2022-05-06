import os
from setuptools import setup
from setuptools.command.build_ext import build_ext
from setuptools.extension import Extension

# Package meta-data.
NAME = 'gimpact'
DESCRIPTION = 'An unofficial python extension for the GImpact collision library.'
URL = 'https://github.com/StephenNneji'
AUTHOR = 'Stephen Nneji'
EMAIL = 'steve.nneji@gmail.com'
VERSION = '0.1.1'

path=os.path.abspath(os.path.dirname(__file__))
try:
    with open(os.path.join(path, 'README.md'), encoding='utf-8') as f:
        long_description = '\n' + f.read()
except FileNotFoundError:
    long_description = DESCRIPTION


class build_ext(build_ext):
    def build_extensions(self):
        import numpy
        try:
            from Cython.Build import cythonize
            BUILD_CYTHON = True
        except ImportError:
            BUILD_CYTHON = False
            
        source = 'gimpact.pyx' if BUILD_CYTHON else 'gimpact.cpp'

        modules = Extension(
            name='gimpact',
            sources=[source, 'GIMPACT/src/gimpact.cpp', 'GIMPACT/src/gim_contact.cpp', 'GIMPACT/src/gim_math.cpp', 
                    'GIMPACT/src/gim_boxpruning.cpp', 'GIMPACT/src/gim_memory.cpp','GIMPACT/src/gim_tri_tri_overlap.cpp',
                    'GIMPACT/src/gim_trimesh.cpp','GIMPACT/src/gim_trimesh_capsule_collision.cpp', 'GIMPACT/src/gim_trimesh_ray_collision.cpp',
                    'GIMPACT/src/gim_trimesh_sphere_collision.cpp', 'GIMPACT/src/gim_trimesh_trimesh_collision.cpp'],
            language='c++',
            include_dirs=['GIMPACT/include', numpy.get_include()],
        )

        self.distribution.ext_modules[:] = cythonize([modules]) if BUILD_CYTHON else [modules]

        self.swig_opts = None
        super().finalize_options()

        super().build_extensions()


setup(name=NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=long_description,
      long_description_content_type="text/markdown",
      url=URL,
      author=AUTHOR,
      author_email=EMAIL,
      license='BSD',
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'Topic :: Games/Entertainment :: Simulation',
          'Topic :: Scientific/Engineering :: Physics',
          'License :: OSI Approved :: BSD License',
          'Operating System :: Microsoft :: Windows',
          'Operating System :: POSIX',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Programming Language :: Python :: Implementation :: CPython',
          'Programming Language :: Cython'
      ],
      keywords='GImpact,Trimesh,Collision detection,Cython',
      cmdclass={"build_ext": build_ext},
      ext_modules=[Extension("", [])],
      test_suite='tests',
      install_requires = ['numpy>=1.16.4'],
      tests_require=['numpy>=1.16.4'],
)
