GImpact-Python
==================
An unofficial python extension for the [GImpact collision library](http://gimpact.sourceforge.net/manual/gimpact_manual.html). This extension integrates directly with GImpact's C++ API using Cython.

Features
---------
* Create trimesh object from numpy array
* Mesh decimation using Sven Forstmann's C++ mesh simplification [code](https://github.com/sp4cerat/Fast-Quadric-Mesh-Simplification) 
* Axis Aligned Bounding Box (AABB)
* AABB set for box prunning
* Collision of triangle mesh with the following
    * triangle mesh
    * sphere
    * capsule
    * plane
    * ray
 * Supports "first contact" or "all contacts" modes

Build Wheel
-----------
``` shell
pip install -q build
python -m build
```


Installation
------------
Build requires numpy and cython (tested on Linux and Windows).
``` shell
pip install gimpact
```

Example Usage
-------------
AABB
```  python
import gimpact


aabb1 = gimpact.AABB(-1, 1, -1, 1, -1, 1)
aabb2 = gimpact.AABB(-1, 1, -1, 1, 1.5, 2)

print(aabb1.intersects(aabb2))
print(aabb1.intersection(aabb2))

aabb1.merge(aabb2)
print(aabb1)
```

Box Prunning
``` python
import gimpact

aabb_set = gimpact.AABBSet(10)
print(len(aabb_set))
print(aabb_set.global_bounds)
for aabb in aabb_set:
    aabb.bounds = (0., 0., 0., 0., 0., 0.)

for aabb in aabb_set:
    print(aabb)

print(aabb_set.global_bounds)
pairs = aabb_set.find_intersections(aabb_set)
print(pairs)

del aabb_set
print(aabb.bounds)
```

Collision
```  python
import gimpact
import numpy as np

contacts = gimpact.trimesh_trimesh_collision(trimesh1, trimesh2)
contacts = gimpact.trimesh_sphere_collision(trimesh1, [0., 0., 0.], 1, True)
contacts = gimpact.trimesh_capsule_collision(trimesh1, np.array([0., 0., 0.]), np.array([1., 0., 0.]), 1, True)
contacts = gimpact.trimesh_plane_collision(trimesh1, [0., 0., 1., 0.], True)
```