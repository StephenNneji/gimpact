from cpython.ref cimport PyObject
from gimpact cimport *
import numpy as np


cdef object __create_uninitialized__ = object()

cdef class AABBSet:
    """Represents a collection of axis-aligned bounding boxes (AABBs)
    
    :param count: number of bounding box
    :type count: int
    """

    cdef GIM_AABB_SET _aabb_set
    
    def __cinit__(self, int count):
        gim_aabbset_alloc(&self._aabb_set, count)	

    def __dealloc__(self):
        gim_aabbset_destroy(&self._aabb_set)

    def __getitem__(self, index):
        if index < 0 or index >= self._aabb_set.m_count:
            raise IndexError('AABBSet index out of range')

        return AABBSet.aabb3f_to_tuple(self._aabb_set.m_boxes[index])

    def __setitem__(self, index, bounds):
        if index < 0 or index >= self._aabb_set.m_count:
            raise IndexError('AABBSet index out of range')

        self._aabb_set.m_boxes[index].minX = bounds[0]
        self._aabb_set.m_boxes[index].maxX = bounds[1]
        self._aabb_set.m_boxes[index].minY = bounds[2]
        self._aabb_set.m_boxes[index].maxY = bounds[3]
        self._aabb_set.m_boxes[index].minZ = bounds[4]
        self._aabb_set.m_boxes[index].maxZ = bounds[5]
        gim_aabbset_update(&self._aabb_set)

    @property
    def global_bound(self):
        """Gets an AABB bounding the entire collection

        :return: bounds of AABB (min X, max X, min Y, max Y, min Z, max Z)
        :rtype: Tuple[float, float, float, float, float, float]
        """
        
        gim_aabbset_update(&self._aabb_set) 
        return AABBSet.aabb3f_to_tuple(self._aabb_set.m_global_bound)

    def find_intersections(self, aabb_set):
        """Finds all intersections between this AABBs of this 
        AABBSet and those of another.

        :param aabb_set: the other AABBSet
        :type aabb_set: AABBSet
        :return: A list of intersecting index pairs
        :rtype: List[Tuple[int]]
        """
        
        cdef GDYNAMIC_ARRAY gim_pairs
        GIM_CREATE_PAIR_SET(gim_pairs)

        cdef object id = aabb_set._id()
        
        gim_aabbset_bipartite_intersections(&self._aabb_set, <GIM_AABB_SET*>PyLong_AsVoidPtr(id), &gim_pairs)
        pairs = []
        
        cdef int count = gim_pairs.m_size
        cdef GIM_PAIR* data = <GIM_PAIR*>gim_pairs.m_pdata 
        
        for pair in data[:count]:
            pairs.append((pair.m_index1, pair.m_index2))
        
        GIM_DYNARRAY_DESTROY(gim_pairs)

        return pairs

    def __len__(self):
        return self._aabb_set.m_count

    @staticmethod
    cdef aabb3f_to_tuple(aabb3f aabb):
        return (aabb.minX, aabb.maxX, aabb.minY, aabb.maxY, aabb.minZ, aabb.maxZ)
    
    def _id(self):
        return PyLong_FromVoidPtr(<void*>&self._aabb_set)


cdef class TriMesh:
    """Represents a triangle mesh
    
    :param vertices: array of 3D vertices
    :type vertices: Array[float]
    :param indices: array of indices
    :type indices: Array[int]
    :param copy: indicates that data should be copied
    :type copy: bool
    """

    cdef GIM_TRIMESH _trimesh
    cdef GBUFFER_MANAGER_DATA buffer_managers[G_BUFFER_MANAGER__MAX]
    
    def __cinit__(self, float[:, ::1] vertices, int[::1] indices, copy=True, flag=None):
        if flag is __create_uninitialized__:
            return
        
        if vertices is None or indices is None:
            raise ValueError('vertices and indices cannot be None.')
        
        cdef char _copy = 1 if copy else 0
        gim_init_buffer_managers(self.buffer_managers)

        if np.amin(indices) < 0 or np.amax(indices) >= vertices.shape[0]:
            raise ValueError("Vertex index out of range")

        gim_trimesh_create_from_data(self.buffer_managers, &self._trimesh, <vec3f*>&vertices[0, 0],
                                     vertices.shape[0], _copy, <GUINT32*>&indices[0], indices.size, _copy, 1)
        gim_trimesh_update(&self._trimesh)

    def __dealloc__(self):
        gim_trimesh_destroy(&self._trimesh)
        gim_terminate_buffer_managers(self.buffer_managers)

    def transform(self, float[:, ::1] matrix not None):
        if matrix.shape[0] != 4 or matrix.shape[1] != 4:
            raise ValueError("Transformation matrix should have dimension (4, 4)")

        gim_trimesh_set_tranform(&self._trimesh, <mat4f>&matrix[0, 0])
        gim_trimesh_update(&self._trimesh)
    
    def clone(self):
        """Returns a  clone of the trimesh.
        
        :return: trimesh clone
        :rtype: Trimesh
        """
        cdef GIM_TRIMESH dest_trimesh
        cdef GBUFFER_MANAGER_DATA buffer_managers[G_BUFFER_MANAGER__MAX]
        gim_init_buffer_managers(buffer_managers)

        gim_trimesh_copy(&self._trimesh, buffer_managers, &dest_trimesh, 0, 1)
        trimesh = TriMesh(None, None, flag=__create_uninitialized__)
        trimesh._trimesh = dest_trimesh
        trimesh.buffer_managers = buffer_managers
        
        return trimesh
    
    def getTriangleCount(self):
        """Returns the number of triangles in the TriMesh.
        
        :return: the number of triangles in the TriMesh
        :rtype: int
        """

        return gim_trimesh_get_triangle_count(&self._trimesh)

    def getTriangle(self, int idx):
        """returns vertices of the triangle with specified index

        :param idx: triangle index
        :type idx: int
        :return: vertices of triangle with specified index
        :rtype: Tuple[Tuple[float, float, float]]
        """
        cdef vec3f v0 = [0., 0., 0.]
        cdef vec3f v1 = [0., 0., 0.]
        cdef vec3f v2 = [0., 0., 0.]
        
        gim_trimesh_locks_work_data(&self._trimesh)
        gim_trimesh_get_triangle_vertices(&self._trimesh, idx, v0, v1, v2)
        gim_trimesh_unlocks_work_data(&self._trimesh)

        return ((v0[0], v0[1], v0[2]),
                (v1[0], v1[1], v1[2]),
                (v2[0], v2[1], v2[2]))
    
    def _id(self):
        return PyLong_FromVoidPtr(<void*>&self._trimesh)


def trimesh_trimesh_collision(trimesh1, trimesh2):
    """Determines contacts of a trimesh-trimesh collision
    
    :param trimesh1: first triangle mesh
    :type trimesh1: Trimesh
    :param trimesh2: second triangle mesh
    :type trimesh2: Trimesh
    :return: a list of Contacts
    :rtype: List[Contact]
    """
    
    cdef GDYNAMIC_ARRAY gim_contacts
    GIM_CREATE_CONTACT_LIST(gim_contacts)

    cdef object id1 = trimesh1._id()
    cdef object id2 = trimesh2._id()
    gim_trimesh_trimesh_collision(<GIM_TRIMESH*>PyLong_AsVoidPtr(id1), <GIM_TRIMESH*>PyLong_AsVoidPtr(id2), &gim_contacts)
       
    return extract_contact_data(gim_contacts)


def trimesh_sphere_collision(trimesh, float[::1] center, float radius):
    """Determines contacts of a trimesh-sphere collision
    
    :param trimesh: triangle mesh
    :type trimesh: Trimesh
    :param center: center of the sphere
    :type center: Array[float]
    :param radius: radius of the sphere
    :type radius: float
    :return: a list of Contacts
    :rtype: List[Contact]
    """
    cdef GDYNAMIC_ARRAY gim_contacts
    GIM_CREATE_CONTACT_LIST(gim_contacts)

    cdef object id = trimesh._id()
    gim_trimesh_sphere_collision(<GIM_TRIMESH*>PyLong_AsVoidPtr(id), <vec3f>&center[0], radius, &gim_contacts)
       
    return extract_contact_data(gim_contacts)


def trimesh_plane_collision(trimesh, float[::1] plane):
    """Determines contacts of a trimesh-plane collision
    
    :param trimesh: triangle mesh
    :type trimesh: Trimesh
    :param plane: plane parameters (a, b, c, d) in form of ax + by + cz + d = 0
    :type plane: Array[float]
    :return: a list of tuples containing point and penetration depth
    :rtype: List[Tuple[Array[float], float]]
    """
    
    cdef GDYNAMIC_ARRAY gim_contacts
    GIM_CREATE_TRIMESHPLANE_CONTACTS(gim_contacts)

    cdef object id = trimesh._id()
    
    gim_trimesh_plane_collision(<GIM_TRIMESH*>PyLong_AsVoidPtr(id), <vec4f>&plane[0],  &gim_contacts)
    contacts = []
    
    cdef int count = gim_contacts.m_size
    cdef vec4f* _contacts_data = <vec4f*>gim_contacts.m_pdata 
    
    for c in _contacts_data[:count]:
        t = np.array(<vec4f>c)
        contacts.append((t[0:3], t[3]))
    
    GIM_DYNARRAY_DESTROY(gim_contacts)

    return contacts


def trimesh_capsule_collision(trimesh, float[::1] point1, float[::1] point2, float radius):
    """Determines contacts of a trimesh-capsule collision
    
    :param trimesh: triangle mesh
    :type trimesh: Trimesh
    :param point1: first end-point of the capsule
    :type point1: Array[float]
    :param point2: second end-point of the capsule
    :type point2: Array[float]
    :param radius: radius of the sphere
    :type radius: float
    :return: a list of Contacts
    :rtype: List[Contact]
    """
    cdef GDYNAMIC_ARRAY gim_contacts
    GIM_CREATE_CONTACT_LIST(gim_contacts)

    cdef object id = trimesh._id()
    cdef GIM_CAPSULE_DATA capsule
    capsule.m_radius = radius
    capsule.m_point1 = <vec3f>&point1[0]
    capsule.m_point2 = <vec3f>&point2[0]
    
    gim_trimesh_capsule_collision(<GIM_TRIMESH*>PyLong_AsVoidPtr(id), &capsule, &gim_contacts)
       
    return extract_contact_data(gim_contacts)


def trimesh_ray_collision(trimesh, float[::1] origin, float[::1] direction, float tmax):
    """Determines contact of a trimesh-ray collision
    
    :param trimesh: triangle mesh
    :type trimesh: Trimesh
    :param origin: origin point of ray
    :type origin: Array[float]
    :param direction: direction vector of ray
    :type direction: Array[float]
    :param tmax: max distance param for ray.
    :type tmax: float
    :return: random contact if ray collides else None
    :rtype: Union[Contact, None]
    """

    cdef GIM_TRIANGLE_RAY_CONTACT_DATA gim_contact
    intersect = gim_trimesh_ray_collision(<GIM_TRIMESH*>PyLong_AsVoidPtr(trimesh._id()),
                                          <vec3f>&origin[0], <vec3f>&direction[0], tmax, &gim_contact)

    if intersect:
        pos= np.array(gim_contact.m_point)
        normal = np.array(gim_contact.m_normal)
        return Contact(pos, normal, gim_contact.tparam, gim_contact.m_face_id, None)


def trimesh_ray_closest_collision(trimesh, float[::1] origin, float[::1] direction, float tmax):
    """Determines closest contact of a trimesh-ray collision
    
    :param trimesh: triangle mesh
    :type trimesh: Trimesh
    :param origin: origin point of ray
    :type origin: Array[float]
    :param direction: direction vector of ray
    :type direction: Array[float]
    :param tmax: max distance param for ray.
    :type tmax: float
    :return: closest contact if ray collides else None
    :rtype: Union[Contact, None]
    """
    
    cdef GIM_TRIANGLE_RAY_CONTACT_DATA gim_contact
    intersect = gim_trimesh_ray_closest_collision(<GIM_TRIMESH*>PyLong_AsVoidPtr(trimesh._id()),
                                                  <vec3f>&origin[0], <vec3f>&direction[0], tmax, &gim_contact)

    if intersect:
        pos= np.array(gim_contact.m_point)
        normal = np.array(gim_contact.m_normal)
        return Contact(pos, normal, gim_contact.tparam, gim_contact.m_face_id, None)


cdef extract_contact_data(GDYNAMIC_ARRAY  gim_contacts):
    contacts = []
    
    cdef int count = gim_contacts.m_size
    cdef GIM_CONTACT* _contacts_data = <GIM_CONTACT*>gim_contacts.m_pdata 
    
    for c in _contacts_data[:count]:
        pos= np.array(c.m_point)
        normal = np.array(c.m_normal)
        contacts.append(Contact(pos, normal, c.m_depth, c.m_feature1, c.m_feature2))
    
    GIM_DYNARRAY_DESTROY(gim_contacts) 

    return contacts


class Contact:
    """Represents a collision contact
    
    :param point: intersection point
    :type point: Array[float]
    :param normal: friction direction
    :type normal: Array[float]
    :param depth: peneration depth
    :type depth: float
    :param feature1: index of colliding triangle in first mesh.
    :type feature1: int
    :param feature2: index of colliding triangle in second mesh.
    :type feature2: int
    """
    def __init__(self, point, normal, depth, feature1, feature2):
        self.point = point
        self.normal = normal
        self.depth = depth
        self.feature1 = feature1
        self.feature2 = feature2
    
    def __str__(self):
        return (f'Contact:\n\tpoint: {self.point}\n\tnormal: {self.normal}\n\tdepth: {self.depth}'
                f'\n\tfeature1: {self.feature1}, \n\tfeature2: {self.feature2}')

    def __repr__(self):
        return str(self)


def terminate():
    """terminate()"""

    gimpact_terminate()


def initialize():
    """init()"""

    gimpact_init()


initialize()
