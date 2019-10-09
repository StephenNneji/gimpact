from cpython.ref cimport PyObject
from gimpact cimport *
import numpy as np
cimport numpy as np

cdef object __create_uninitialized__ = object()


cdef class AABBSet:
    """Represents a collection of axis-aligned bounding boxes (AABBs)
    
    :param count: number of bounding box
    :type count: int
    """

    cdef GIM_AABB_SET _aabb_set
    
    def __cinit__(self, int count, flag=None):
        if flag is __create_uninitialized__:
            return

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
    def global_bounds(self):
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
    cdef object _aabb_set
    
    def __cinit__(self, vertices, indices, flag=None):
        if flag is __create_uninitialized__:
            return
        
        if vertices is None or indices is None:
            raise ValueError('vertices and indices cannot be None.')     
        
        cdef float[:, ::1] _vertices = np.array(vertices, dtype=np.float32, copy=False, order='C',)
        cdef int[::1] _indices = np.array(indices, dtype=np.int32, copy=False, order='C',)
        
        if np.amin(indices) < 0 or np.amax(indices) >= vertices.shape[0]:
            raise ValueError("Vertex index out of range")
        
        gim_init_buffer_managers(self.buffer_managers)
        gim_trimesh_create_from_data(self.buffer_managers, &self._trimesh, <vec3f*>&_vertices[0, 0],
                                     _vertices.shape[0], 1, <GUINT32*>&_indices[0], _indices.size, 1, 1)
        gim_trimesh_update(&self._trimesh)

        self._aabb_set = None

    def __dealloc__(self):
        gim_trimesh_destroy(&self._trimesh)
        gim_terminate_buffer_managers(self.buffer_managers)

    def transform(self, matrix not None):
        cdef float[:, ::1] _matrix = np.array(matrix, dtype=np.float32, copy=False, order='C',)
        if _matrix.shape[0] != 4 or _matrix.shape[1] != 4:
            raise ValueError("Transformation _matrix should have dimension (4, 4)")

        gim_trimesh_set_tranform(&self._trimesh, <mat4f>&_matrix[0, 0])
        gim_trimesh_update(&self._trimesh)

    @property
    def aabb_set(self):
        if self._aabb_set is None:
            self._trimesh.m_aabbset.m_shared = 1
            aabb_set = AABBSet(0, flag=__create_uninitialized__)
            aabb_set._aabb_set = self._trimesh.m_aabbset
            self._aabb_set = aabb_set

        return self._aabb_set

    @property
    def bounds(self):
        """Gets an AABB bounding the trimesh

        :return: bounds of AABB (min X, max X, min Y, max Y, min Z, max Z)
        :rtype: Tuple[float, float, float, float, float, float]
        """
        return self.aabb_set.global_bounds

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
    
    @property
    def triangle_count(self):
        """Returns the number of triangles in the TriMesh.
        
        :return: the number of triangles in the TriMesh
        :rtype: int
        """

        return gim_trimesh_get_triangle_count(&self._trimesh)

    def triangle(self, int idx):
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
    
    def decimate(self, int target_count):
        cdef GBUFFER_ARRAY* vertex_buffer = &self._trimesh.m_source_vertex_buffer
        cdef GBUFFER_ARRAY* index_buffer = &self._trimesh.m_tri_index_buffer
        gim_buffer_array_lock(vertex_buffer, 1)
        gim_buffer_array_lock(index_buffer, 1)
        vertices = np.asarray(<float[:vertex_buffer.m_element_count, :3]><float*>vertex_buffer.m_buffer_data)
        indices = np.asarray(<int[:index_buffer.m_element_count]><int*>index_buffer.m_buffer_data)
        tmp = np.unique(vertices[indices], return_inverse=True, axis=0)
        cdef float[:, ::1] _vertices = tmp[0]
        cdef int[::1] _indices = tmp[1].astype(np.int32)  # np.unique returns int64
        read_mesh(&_vertices[0, 0], _vertices.shape[0], &_indices[0], indices.size)
        gim_buffer_array_unlock(vertex_buffer)
        gim_buffer_array_unlock(index_buffer)
        
        simplify_mesh(target_count, 7, False)
        
        cdef vector[float] v
        cdef vector[int] i
        write_mesh(v, i)

        return TriMesh(np.array(v, np.float32).reshape(-1, 3), np.array(i, np.int32))

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


def trimesh_sphere_collision(trimesh, center, float radius):
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
    cdef vec3f _center = [center[0], center[1], center[2]]
    gim_trimesh_sphere_collision(<GIM_TRIMESH*>PyLong_AsVoidPtr(id), _center, radius, &gim_contacts)
       
    return extract_contact_data(gim_contacts)


def trimesh_plane_collision(trimesh, plane):
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
    cdef vec4f _plane = [plane[0], plane[1], plane[2], plane[3]]
    gim_trimesh_plane_collision(<GIM_TRIMESH*>PyLong_AsVoidPtr(id), _plane,  &gim_contacts)
    contacts = []
    
    cdef int count = gim_contacts.m_size
    cdef vec4f* _contacts_data = <vec4f*>gim_contacts.m_pdata 
    
    for c in _contacts_data[:count]:
        t = np.array(<vec4f>c)
        contacts.append((t[0:3], t[3]))
    
    GIM_DYNARRAY_DESTROY(gim_contacts)

    return contacts


def trimesh_capsule_collision(trimesh, point1, point2, float radius):
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
    capsule.m_point1 = [point1[0], point1[1], point1[2]]
    capsule.m_point2 = [point2[0], point2[1], point2[2]]
    
    gim_trimesh_capsule_collision(<GIM_TRIMESH*>PyLong_AsVoidPtr(id), &capsule, &gim_contacts)
       
    return extract_contact_data(gim_contacts)


def trimesh_ray_collision(trimesh, origin, direction, float tmax):
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
    cdef vec3f _origin = [origin[0], origin[1], origin[2]]
    cdef vec3f _direction = [direction[0], direction[1], direction[2]]
    intersect = gim_trimesh_ray_collision(<GIM_TRIMESH*>PyLong_AsVoidPtr(trimesh._id()),
                                          _origin, _direction, tmax, &gim_contact)

    if intersect:
        pos= np.array(gim_contact.m_point)
        normal = np.array(gim_contact.m_normal)
        return Contact(pos, normal, gim_contact.tparam, gim_contact.m_face_id, None)


def trimesh_ray_closest_collision(trimesh, origin, direction, float tmax):
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
    cdef vec3f _origin = [origin[0], origin[1], origin[2]]
    cdef vec3f _direction = [direction[0], direction[1], direction[2]]
    intersect = gim_trimesh_ray_closest_collision(<GIM_TRIMESH*>PyLong_AsVoidPtr(trimesh._id()),
                                                  _origin, _direction, tmax, &gim_contact)

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


def decimate(vertices, indices, int target_count):
    
    cdef float[:, ::1] _vertices = np.array(vertices, dtype=np.float32, copy=False, order='C',)
    cdef int[::1] _indices = np.array(indices, dtype=np.int32, copy=False, order='C',)
        
    if np.amin(indices) < 0 or np.amax(indices) >= vertices.shape[0]:
        raise ValueError("Vertex index out of range")
    
    cdef vector[float] v
    cdef vector[int] i
    
    read_mesh(&_vertices[0, 0], _vertices.shape[0], &_indices[0], indices.size)
    simplify_mesh(target_count, 7, False)
    write_mesh(v, i)
    
    return np.array(v, np.float32).reshape(-1, 3), np.array(i, np.int32)


cdef terminate():
    """terminate()"""

    gimpact_terminate()


cdef initialize():
    """init()"""

    gimpact_init()


initialize()
