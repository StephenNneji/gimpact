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
