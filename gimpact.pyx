from cpython.ref cimport PyObject
from gimpact cimport *
import numpy as np


cdef object __create_uninitialized__ = object()


cdef class TriMesh:
    """TriMesh object."""

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
    
    def copy(self):
        cdef GIM_TRIMESH dest_trimesh
        cdef GBUFFER_MANAGER_DATA buffer_managers[G_BUFFER_MANAGER__MAX]
        gim_init_buffer_managers(buffer_managers)

        gim_trimesh_copy(&self._trimesh, buffer_managers, &dest_trimesh, 0, 1)
        trimesh = TriMesh(None, None, flag=__create_uninitialized__)
        trimesh._trimesh = dest_trimesh
        trimesh.buffer_managers = buffer_managers
        
        return trimesh

    
    def getTriangleCount(self):
        """getTriangleCount() -> n

        Returns the number of triangles in the TriMesh."""

        return gim_trimesh_get_triangle_count(&self._trimesh)

    def getTriangle(self, int idx):
        """getTriangle(idx) -> (v0, v1, v2)

        @param idx: Triangle index
        @type idx: int
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
    cdef GDYNAMIC_ARRAY gim_contacts
    GIM_CREATE_CONTACT_LIST(gim_contacts)

    cdef object id1 = trimesh1._id()
    cdef object id2 = trimesh2._id()
    gim_trimesh_trimesh_collision(<GIM_TRIMESH*>PyLong_AsVoidPtr(id1), <GIM_TRIMESH*>PyLong_AsVoidPtr(id2), &gim_contacts)
       
    return extract_contact_data(gim_contacts)


def trimesh_sphere_collision(trimesh, float[::1] center, float radius):
    cdef GDYNAMIC_ARRAY gim_contacts
    GIM_CREATE_CONTACT_LIST(gim_contacts)

    cdef object id = trimesh._id()
    gim_trimesh_sphere_collision(<GIM_TRIMESH*>PyLong_AsVoidPtr(id), <vec3f>&center[0], radius, &gim_contacts)
       
    return extract_contact_data(gim_contacts)


def trimesh_plane_collision(trimesh, float[::1] plane):
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
    cdef GIM_TRIANGLE_RAY_CONTACT_DATA gim_contact
    intersect = gim_trimesh_ray_collision(<GIM_TRIMESH*>PyLong_AsVoidPtr(trimesh._id()),
                                          <vec3f>&origin[0], <vec3f>&direction[0], tmax, &gim_contact)

    if intersect:
        pos= np.array(gim_contact.m_point)
        normal = np.array(gim_contact.m_normal)
        return Contact(pos, normal, gim_contact.tparam, gim_contact.m_face_id, None)


def trimesh_ray_closest_collision(trimesh, float[::1] origin, float[::1] direction, float tmax):
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
