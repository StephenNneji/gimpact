import gimpact
import numpy as np


def read_binary_stl(filename):
    with open(filename, 'rb') as f:
        f.seek(80)
        face_count = np.frombuffer(f.read(4),
                                   dtype=np.int32)[0]

        record_dtype = np.dtype([
            ('normals', np.float32, (3,)),
            ('vertices', np.float32, (3, 3)),
            ('attr', '<i2', (1,)),
        ])
        data = np.fromfile(f, dtype=record_dtype)

    if face_count != data.size:
        raise IOError('stl data has incorrect size')

    points = data['vertices'].reshape(-1, 3)
    indices = np.arange(face_count * 3)

    return {'vertices': points, 'indices': indices}

# Testing AABB Set

aabb_set = gimpact.AABBSet(10)
print(len(aabb_set))
print(aabb_set.global_bound)
for i in range(len(aabb_set)):
    aabb_set[i] = (0., 0., 0., 0., 0., 0.)

for aabb in aabb_set:
    print(aabb)

print('\n')
print(aabb_set.global_bound)
pairs = aabb_set.find_intersections(aabb_set)
print(pairs)

# Testing Trimesh Class

mesh = read_binary_stl(r'D:\Experiments\RB1910414\M14\M14.stl')
trimesh = gimpact.TriMesh(mesh['vertices'], mesh['indices'])
v = mesh['vertices']
print(trimesh.getTriangleCount())
print(trimesh.getTriangle(0))
del mesh
#t = np.identity(4, np.float32)
#t[2, 3] = -1500
#trimesh.transform(t)
print(trimesh.getTriangle(0))
v[0, 0] = -1
print('\n', v[0])
print('\n', trimesh.getTriangle(0))

mesh = read_binary_stl(r'D:\Software Development\SScanSS-2\Code\SScanSS-2\instruments\engin-x\models\z_stage.stl')
trimesh1 = gimpact.TriMesh(mesh['vertices'], mesh['indices'])
#trimesh1.transform(t)
mesh = read_binary_stl(r'D:\Software Development\SScanSS-2\Code\SScanSS-2\instruments\engin-x\models\y_stage.stl')
trimesh2 = gimpact.TriMesh(mesh['vertices'], mesh['indices'])
#trimesh2.transform(t)

print(trimesh1.getTriangleCount())
for i in range(trimesh1.getTriangleCount()):
    print(trimesh1.getTriangle(i))

print(trimesh2.getTriangleCount())
for i in range(trimesh2.getTriangleCount()):
    print(trimesh2.getTriangle(i))

contacts = gimpact.trimesh_trimesh_collision(trimesh1, trimesh2)
contacts = gimpact.trimesh_sphere_collision(trimesh1, np.array([0., 0., 0.], np.float32), 1)
contacts = gimpact.trimesh_capsule_collision(trimesh1, np.array([-1., 0., 0.], np.float32), np.array([1., 0., 0.], np.float32), 1)
contacts = gimpact.trimesh_plane_collision(trimesh1, np.array([0., 0., 0., 0.], np.float32))
for c in contacts:
   print(*c)

contact = gimpact.trimesh_ray_collision(trimesh2, np.array([0., 0., 0.], np.float32), np.array([-1., 0., 0.], np.float32), 1000)
print('\n', contact)


contact = gimpact.trimesh_ray_closest_collision(trimesh2, np.array([0., 0., 0.], np.float32), np.array([-1., 0., 0.], np.float32), 1000)
print('\n', contact)

t = trimesh.clone()
print(t.getTriangleCount())