cdef extern from "stdlib.h":

    void* malloc(long)
    void free(void*)


cdef extern from "stdio.h":
    int printf(char*, ...)


cdef extern from "Python.h":
    object PyLong_FromVoidPtr(void *)
    void *PyLong_AsVoidPtr(object)


cdef extern from "GIMPACT/gimpact.h":
    void gimpact_init()
    void gimpact_terminate()


cdef extern from "GIMPACT/gim_math.h":
    ctypedef unsigned int GUINT32
    ctypedef float GREAL
    ctypedef GREAL vec3f[3]
    ctypedef GREAL vec4f[4]
    ctypedef GREAL mat4f[4][4]


cdef extern from "GIMPACT/gim_contact.h":
    cdef struct GIM_CONTACT:
        vec3f m_point
        vec3f m_normal
        GREAL m_depth
        void * m_handle1
        void * m_handle2
        GUINT32 m_feature1
        GUINT32 m_feature2

    cdef struct GIM_TRIANGLE_RAY_CONTACT_DATA:
        GREAL u
        GREAL v
        GREAL tparam
        GUINT32 m_face_id
        vec3f m_point
        vec3f m_normal

    void GIM_CREATE_CONTACT_LIST(GDYNAMIC_ARRAY contact_array)


cdef extern from "GIMPACT/gim_memory.h":
    cpdef enum:
        G_BUFFER_MANAGER_SYSTEM,
        G_BUFFER_MANAGER_SHARED,
        G_BUFFER_MANAGER__MAX   
    
    cdef struct GBUFFER_MANAGER_DATA:
        pass
    
    cdef struct GDYNAMIC_ARRAY:
        char * m_pdata
        GUINT32 m_size
        GUINT32 m_reserve_size
    
    void * gim_alloc(size_t size)

    void gim_init_buffer_managers(GBUFFER_MANAGER_DATA buffer_managers[])
    
    void gim_terminate_buffer_managers(GBUFFER_MANAGER_DATA buffer_managers[])

    void GIM_DYNARRAY_DESTROY(GDYNAMIC_ARRAY & array_data)


cdef extern from "GIMPACT/gim_boxpruning.h":
    cdef struct aabb3f:
        GREAL 	minX
        GREAL 	maxX
        GREAL 	minY
        GREAL 	maxY
        GREAL 	minZ
        GREAL 	maxZ

    struct GIM_RSORT_TOKEN:
        pass

    cdef struct GIM_AABB_SET:
        GUINT32 m_count
        aabb3f m_global_bound
        aabb3f * m_boxes
        GUINT32 * m_maxcoords
        GIM_RSORT_TOKEN * m_sorted_mincoords;
        char m_shared;
    
    cdef struct GIM_PAIR:
        GUINT32 m_index1
        GUINT32 m_index2

    void gim_aabbset_alloc(GIM_AABB_SET * aabbset, GUINT32 count)
    
    void gim_aabbset_destroy(GIM_AABB_SET * aabbset)

    void gim_aabbset_update(GIM_AABB_SET * aabbset) 	

    void gim_aabbset_bipartite_intersections(GIM_AABB_SET * aabbset1, GIM_AABB_SET * aabbset2, GDYNAMIC_ARRAY * collision_pairs)
    
    void GIM_CREATE_PAIR_SET(GDYNAMIC_ARRAY pair_array)


cdef extern from "GIMPACT/gim_tri_capsule_collision.h":
    cdef struct GIM_CAPSULE_DATA:
        GREAL m_radius
        vec3f m_point1
        vec3f m_point2


cdef extern from "GIMPACT/gim_trimesh.h":   
    cdef struct GIM_TRIMESH:
        pass

    void GIM_CREATE_TRIMESHPLANE_CONTACTS(GDYNAMIC_ARRAY contact_array)
    
    void gim_trimesh_create_from_data(GBUFFER_MANAGER_DATA buffer_managers[],
	            GIM_TRIMESH * trimesh, vec3f * vertex_array, GUINT32 vertex_count, char copy_vertices, 
	            GUINT32 * triindex_array, GUINT32 index_count,char copy_indices,char transformed_reply)
    
    void gim_trimesh_destroy(GIM_TRIMESH * trimesh)

    GUINT32 gim_trimesh_get_triangle_count(GIM_TRIMESH * trimesh)
    
    void gim_trimesh_set_tranform(GIM_TRIMESH * trimesh, mat4f transform)

    void gim_trimesh_get_triangle_vertices(GIM_TRIMESH * trimesh, GUINT32 triangle_index, vec3f v1,vec3f v2,vec3f v3)

    void gim_trimesh_locks_work_data(GIM_TRIMESH * trimesh)

    void gim_trimesh_unlocks_work_data(GIM_TRIMESH * trimesh)

    void gim_trimesh_copy(GIM_TRIMESH * source_trimesh, GBUFFER_MANAGER_DATA dest_buffer_managers[], GIM_TRIMESH * dest_trimesh, 
                          char copy_by_reference, char transformed_reply)	 		

    void gim_trimesh_update(GIM_TRIMESH * trimesh)

    void gim_trimesh_trimesh_collision(GIM_TRIMESH * trimesh1, GIM_TRIMESH * trimesh2, GDYNAMIC_ARRAY * contacts)
    
    void gim_trimesh_sphere_collision(GIM_TRIMESH * trimesh,vec3f center,GREAL radius, GDYNAMIC_ARRAY * contacts)
    
    void gim_trimesh_capsule_collision(GIM_TRIMESH * trimesh, GIM_CAPSULE_DATA * capsule, GDYNAMIC_ARRAY * contacts)

    void gim_trimesh_plane_collision(GIM_TRIMESH * trimesh, vec4f plane, GDYNAMIC_ARRAY * contacts)

    int gim_trimesh_ray_collision(GIM_TRIMESH * trimesh,vec3f origin,vec3f dir, GREAL tmax, GIM_TRIANGLE_RAY_CONTACT_DATA * contact)
    
    int gim_trimesh_ray_closest_collision(GIM_TRIMESH * trimesh,vec3f origin,vec3f dir, GREAL tmax, GIM_TRIANGLE_RAY_CONTACT_DATA * contact)
