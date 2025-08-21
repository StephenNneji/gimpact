import os
import sys
from setuptools import setup
from setuptools.command.build_ext import build_ext
from setuptools.extension import Extension

# Package meta-data.
NAME = 'gimpact'
DESCRIPTION = 'An unofficial python extension for the GImpact collision library.'
URL = 'https://github.com/StephenNneji'
AUTHOR = 'Stephen Nneji'
EMAIL = 'steve.nneji@gmail.com'
VERSION = '1.0.2'

path=os.path.abspath(os.path.dirname(__file__))
try:
    with open(os.path.join(path, 'README.md'), encoding='utf-8') as f:
        long_description = '\n' + f.read()
except FileNotFoundError:
    long_description = DESCRIPTION


# check whether compiler supports a flag
def has_flag(compiler, flag_name):
    import tempfile

    from setuptools.errors import CompileError

    with tempfile.NamedTemporaryFile("w", suffix=".cpp") as f:
        f.write("int main (int argc, char **argv) { return 0; }")
        try:
            compiler.compile([f.name], extra_postargs=[flag_name])
        except CompileError:
            return False
    return True


class build_ext(build_ext):
    """A custom build extension for adding compiler-specific options."""

    c_opts = {
        "msvc": ["/O2", "/EHsc"],
        "unix": ["-O2", "-std=c++11"],
    }
    if sys.platform == "darwin":
        darwin_opts = ["-stdlib=libc++", "-mmacosx-version-min=10.9", "-ffp-contract=off"]
        c_opts["unix"] = [*darwin_opts, "-O2"]


    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])

        if ct == "unix":
            if "-Wstrict-prototypes" in self.compiler.compiler_so:
                self.compiler.compiler_so.remove("-Wstrict-prototypes")

            opts.append(f'-DVERSION_INFO="{self.distribution.get_version()}"')
            if has_flag(self.compiler, "-fvisibility=hidden"):
                opts.append("-fvisibility=hidden")
        
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
            extra_compile_args=opts,
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
      cmdclass={"build_ext": build_ext},
      ext_modules=[Extension("", [])],
      install_requires=[
          "numpy >= 1.18.5",
      ],
      test_suite='tests'
)
