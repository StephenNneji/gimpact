import unittest
import numpy as np
import gimpact


def create_cubiod(width=10.0, height=10.0, depth=10.0):

    # centre cube at the origin
    w = width / 2
    h = height / 2
    d = depth / 2
    
    ftl = [-w, h, d] # front top left vertex
    ftr = [w, h, d] # front top right vertex
    fbl = [-w, -h, d] # front bottom left vertex
    fbr = [w, -h, d] # front bottom right vertex
    btl = [-w, h, -d] # back top left vertex
    btr = [w, h, -d] # back top right vertex
    bbl = [-w, -h, -d] # back bottom left vertex
    bbr = [w, -h, -d] # back bottom right vertex
    
    vertices = [fbl, fbr, ftl, ftr, ftl, ftr, 
                btl, btr, btl, btr, bbl, bbr, 
                bbl, bbr, fbl, fbr, fbr,  bbr, 
                ftr, btr, bbl, fbl, btl, ftl]

    indices =  np.array([0, 1, 2, 2, 1, 3, 
                4, 5, 6, 6, 5, 7, 
                8, 9, 10, 10, 9, 11, 
                12, 13, 14, 14, 13, 15, 
                16, 17, 18, 18, 17, 19, 
                20, 21, 22, 22, 21, 23])
    
    return {'vertices': np.array(vertices), 'indices': indices}


class TestAABB(unittest.TestCase):
    def testAABBExtents(self):
        aabb = gimpact.AABB(-1, 1, -1, 1, -1, 1)
        self.assertAlmostEqual(aabb.min_x, -1, 2)
        self.assertAlmostEqual(aabb.max_x, 1, 2)
        self.assertAlmostEqual(aabb.min_y, -1, 2)
        self.assertAlmostEqual(aabb.max_y, 1, 2)
        self.assertAlmostEqual(aabb.min_z, -1, 2)
        self.assertAlmostEqual(aabb.max_z, 1, 2)

        aabb.min_x = -0.5
        aabb.max_x = 0.5
        aabb.min_y = -0.5
        aabb.max_y = 0.5
        aabb.min_z = -0.5
        aabb.max_z = 0.5
        np.testing.assert_array_almost_equal(aabb.bounds, (-0.5, 0.5, -0.5, 0.5, -0.5, 0.5), decimal=2)
        aabb.bounds = (1, -2, 3, -4, -5, 6)
        np.testing.assert_array_almost_equal(aabb.bounds, (1, -2, 3, -4, -5, 6), decimal=2)
        self.assertAlmostEqual(aabb.min_x, 1, 2)
        self.assertAlmostEqual(aabb.max_x, -2, 2)
        self.assertAlmostEqual(aabb.min_y, 3, 2)
        self.assertAlmostEqual(aabb.max_y, -4, 2)
        self.assertAlmostEqual(aabb.min_z, -5, 2)
        self.assertAlmostEqual(aabb.max_z, 6, 2)

    def testAABBMerge(self):
        aabb1 = gimpact.AABB(-1, 1, -1, 1, -1, 1)
        aabb2 = gimpact.AABB(-1, 1, -1, 1, 1.5, 2)
        aabb1.merge(aabb2)
        np.testing.assert_array_almost_equal(aabb1.bounds, (-1, 1, -1, 1, -1, 2), decimal=2)

    def testAABBIntersection(self):
        aabb1 = gimpact.AABB(-1, 1, -1, 1, -1, 1)
        aabb2 = gimpact.AABB(-1, 1, -1, 1, 1.5, 2)
        self.assertFalse(aabb1.intersects(aabb2))
        
        aabb1 = gimpact.AABB(-1, 1, -1, 0.5, -1, 2.5)
        aabb2 = gimpact.AABB(-1, 1, -1, 1, -1.5, 2)
        self.assertTrue(aabb1.intersects(aabb2))
        np.testing.assert_array_almost_equal(aabb1.intersection(aabb2).bounds, (-1, 1, -1, 0.5, -1, 2), decimal=2)
    
    def testClone(self):
        aabb1 = gimpact.AABB(-1, 1, -1, 1, -1, 1)
        aabb2 = aabb1.clone()
        np.testing.assert_array_almost_equal(aabb1.bounds, aabb2.bounds, decimal=2)
        aabb2.max_z = 5
        self.assertNotAlmostEqual(aabb1.max_z, aabb2.max_z)


class TestAABBSet(unittest.TestCase):
    def testBounds(self):
        aabb_set = gimpact.AABBSet(10)
        self.assertEqual(len(aabb_set), 10)

        for aabb in aabb_set:
            aabb.bounds = (-9., 0., 0., 0., 0., 10.)
        aabb_set[0] = (0., 6., 0., 0., -5., 0.)
        np.testing.assert_array_almost_equal(aabb_set.global_bounds.bounds, (-9., 6., 0., 0., -5., 10.), decimal=2)

        bounds = aabb.bounds
        del aabb_set
        np.testing.assert_array_almost_equal(aabb.bounds, bounds, decimal=2)

    def testIntersection(self):
        aabb_set = gimpact.AABBSet(2)
        aabb_set[0] = (-1, 1, -1, 1, -1, 1)
        aabb_set[1] = gimpact.AABB(-1, 1, -1, 1, 1.5, 2)

        np.testing.assert_array_almost_equal(aabb_set.global_bounds.bounds, (-1.,  1., -1.,  1., -1.,  2.), decimal=2)

        expected = [(0, 0), (1, 1)]
        pairs = aabb_set.find_intersections(aabb_set)
        self.assertEqual(len(pairs), 2)
        for pair in pairs:
            self.assertIn(pair, expected)
        
        aabb_set2 = gimpact.AABBSet(2)
        aabb_set2[0] = gimpact.AABB(-1, 1, -1, 1, 1.5, 2)
        aabb_set2[1] = (-1, 1, -1, 1, -1, 1)
        expected = [(0, 1), (1, 0)]
        pairs = aabb_set2.find_intersections(aabb_set)
        self.assertEqual(len(pairs), 2)
        for pair in pairs:
            self.assertIn(pair, expected)


class TestTrimesh(unittest.TestCase):
    def setUp(self):
        mesh = create_cubiod()
        self.vertices = mesh['vertices']
        self.indices = mesh['indices']
    
    def testDecimation(self):
        trimesh = gimpact.TriMesh(self.vertices, self.indices)
        trimesh2 = trimesh.decimate(8)
        self.assertEqual(trimesh2.triangle_count, 8)

    def testTriangle(self):
        trimesh = gimpact.TriMesh(self.vertices, self.indices)
        self.assertEqual(trimesh.triangle_count, self.indices.size//3)
        
        for i in range(trimesh.triangle_count):
            np.testing.assert_array_almost_equal(trimesh.triangle(i), 
                                                 self.vertices[self.indices[i*3:i*3+3]], decimal=2)

        trimesh2 = trimesh.clone()
        for i in range(trimesh.triangle_count):
            a = np.array(trimesh.triangle(i))
            b = np.array(trimesh2.triangle(i))
            np.testing.assert_array_almost_equal(a, b, decimal=2)

    def testTransform(self):
        trimesh = gimpact.TriMesh(self.vertices, self.indices)
        t = np.identity(4)
        t[0:3, 3] = [10, 10, 10]
        trimesh.transform(t)
        v = self.vertices[self.indices]
        v = v @ t[0:3, 0:3].transpose() + t[0:3, 3] 
        for i in range(trimesh.triangle_count):
            np.testing.assert_array_almost_equal(trimesh.triangle(i), v[i*3:i*3+3], decimal=2)

        trimesh2 = trimesh.clone()
        t = np.identity(4)
        t[0:2, 0:2] = [[0, -1], [1, 0]]
        trimesh2.transform(t)
        v = self.vertices[self.indices]
        v = v @ t[0:3, 0:3].transpose() 
        for i in range(trimesh2.triangle_count):
            np.testing.assert_array_almost_equal(trimesh2.triangle(i), v[i*3:i*3+3], decimal=2)

    def testAABBSet(self):
        trimesh = gimpact.TriMesh(self.vertices, self.indices)
        aabb_set = trimesh.aabb_set
        np.testing.assert_array_almost_equal(aabb_set.global_bounds.bounds, trimesh.bounds.bounds, decimal=2)

        t = trimesh.clone()
        aabb_set = t.aabb_set
        self.assertEqual(trimesh.triangle_count, t.triangle_count)
        del t
        self.assertEqual(len(trimesh.aabb_set), len(aabb_set))
    
    def testTrimeshCollision(self):
        trimesh1 = gimpact.TriMesh(self.vertices, self.indices)
        trimesh2 = trimesh1.clone()
    
        contacts = gimpact.trimesh_trimesh_collision(trimesh1, trimesh2, True)
        self.assertEqual(len(contacts), 1)

        t = np.identity(4)
        t[2, 3] = 11
        trimesh2.transform(t)
        contacts = gimpact.trimesh_trimesh_collision(trimesh1, trimesh2)
        self.assertEqual(len(contacts), 0)

        t[:3, 3] = [10, 10, 10]
        trimesh2.transform(t)
        contacts = gimpact.trimesh_trimesh_collision(trimesh1, trimesh2)
        np.testing.assert_array_almost_equal(contacts[0].point, [5, 5, 5], decimal=2)
        np.testing.assert_array_almost_equal(contacts[0].normal, [0, 0, -1], decimal=2)
        self.assertAlmostEqual(contacts[0].depth, 0, 2)
        self.assertEqual(contacts[0].feature1, 1)
        self.assertEqual(contacts[0].feature2, 4)

    def testCapsuleCollision(self):
        trimesh = gimpact.TriMesh(self.vertices, self.indices)
        
        radius = 1
        start = [-1., 0., 5.]
        stop = [1., 0., 5.]
        contacts = gimpact.trimesh_capsule_collision(trimesh, start, stop, radius, True)
        self.assertEqual(len(contacts), 1)
        
        contacts = gimpact.trimesh_capsule_collision(trimesh, start, stop, radius)
        self.assertEqual(len(contacts), 3)
        np.testing.assert_array_almost_equal(contacts[0].point, [-1, 0, 4], decimal=2)
        np.testing.assert_array_almost_equal(contacts[0].normal, [0, 0, -1], decimal=2)
        self.assertAlmostEqual(contacts[0].depth, 1, 2)
        self.assertEqual(contacts[0].feature1, 0)

        np.testing.assert_array_almost_equal(contacts[1].point, [0, 0, 4], decimal=2)
        np.testing.assert_array_almost_equal(contacts[1].normal, [0, 0, -1], decimal=2)
        self.assertAlmostEqual(contacts[1].depth, 1, 2)
        self.assertEqual(contacts[1].feature1, 1)
 
        np.testing.assert_array_almost_equal(contacts[2].point, [1, 0, 4], decimal=2)
        np.testing.assert_array_almost_equal(contacts[2].normal, [0, 0, -1], decimal=2)
        self.assertAlmostEqual(contacts[2].depth, 1, 2)
        self.assertEqual(contacts[2].feature1, 1)

    def testPlaneCollision(self):
        trimesh = gimpact.TriMesh(self.vertices, self.indices)

        contacts = gimpact.trimesh_plane_collision(trimesh, np.array([0., 0., 1., 0.]), True)
        self.assertEqual(len(contacts), 1)
        contacts = gimpact.trimesh_plane_collision(trimesh, np.array([0., 0., 1., 0.]))
        self.assertEqual(len(contacts), self.vertices.shape[0]/2)
        for contact, depth in contacts:
            self.assertAlmostEqual(contact[2], -5, 2)
            self.assertAlmostEqual(depth, 5, 2)

        contacts = gimpact.trimesh_plane_collision(trimesh, np.array([0., 0., 1., -5.01]))
        self.assertEqual(len(contacts), 0)
        contacts = gimpact.trimesh_plane_collision(trimesh, np.array([0., 0., -1., 5.01]))
        self.assertEqual(len(contacts), self.vertices.shape[0])     

    def testRayCollision(self):
        trimesh = gimpact.TriMesh(self.vertices, self.indices)
        start = [0., 0., 0.]
        direction = [-1., 0., 0.]
        depth = 100

        contact = gimpact.trimesh_ray_collision(trimesh, start, direction, depth)
        # No contact because ray originates from inside the object
        self.assertIsNone(contact)
        
        contact = gimpact.trimesh_ray_collision(trimesh, [0., 5.1, 0.], direction, depth)
        self.assertIsNone(contact)
        
        contact = gimpact.trimesh_ray_collision(trimesh, [2.5, -6.0, 0.], [0., 1., 0.], depth)
        np.testing.assert_array_almost_equal(contact.point, [2.5, -5, 0], decimal=2)
        np.testing.assert_array_almost_equal(contact.normal, [0, 1, 0], decimal=2)
        self.assertAlmostEqual(contact.depth, 1, 2)
        self.assertEqual(contact.feature1, 7)

        contact = gimpact.trimesh_ray_closest_collision(trimesh, np.array([5., 0., 5.]), np.array([0., 0., -1.]), 100)
        np.testing.assert_array_almost_equal(contact.point, [5, 0, 5], decimal=2)
        np.testing.assert_array_almost_equal(contact.normal, [0, 0, -1], decimal=2)
        self.assertAlmostEqual(contact.depth, 0, 2)
        self.assertEqual(contact.feature1, 1)

    def testSphereCollision(self):
        trimesh = gimpact.TriMesh(self.vertices, self.indices)

        radius = 1.0   
        center = [5., 5., 5.]    
        contacts = gimpact.trimesh_sphere_collision(trimesh, center, radius, True)
        self.assertEqual(len(contacts), 1)
        
        contacts = gimpact.trimesh_sphere_collision(trimesh, center, radius)
        self.assertEqual(len(contacts), 3)
        np.testing.assert_array_almost_equal(contacts[0].point, [5, 5, 4], decimal=2)
        np.testing.assert_array_almost_equal(contacts[0].normal, [0, 0, -1], decimal=2)
        self.assertAlmostEqual(contacts[0].depth, 1, 2)
        self.assertEqual(contacts[0].feature1, 1)

        np.testing.assert_array_almost_equal(contacts[1].point, [5, 4, 5], decimal=2)
        np.testing.assert_array_almost_equal(contacts[1].normal, [0, -1, 0], decimal=2)
        self.assertAlmostEqual(contacts[1].depth, 1, 2)
        self.assertEqual(contacts[1].feature1, 3)
 
        np.testing.assert_array_almost_equal(contacts[2].point, [4, 5, 5], decimal=2)
        np.testing.assert_array_almost_equal(contacts[2].normal, [-1, 0, 0], decimal=2)
        self.assertAlmostEqual(contacts[2].depth, 1, 2)
        self.assertEqual(contacts[2].feature1, 8)


if __name__ == '__main__':
    unittest.main()
