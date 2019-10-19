import unittest
import gimpact

class TestAABB(unittest.TestCase):
    def testAABBSet(self):
        aabb_set = gimpact.AABBSet(10)
        print(len(aabb_set))
        print(aabb_set.global_bounds)
        for aabb in aabb_set:
            aabb.bounds = (0., 0., 0., 0., 0., 0.)

        for aabb in aabb_set:
            print(aabb)

        print('\n')
        print(aabb_set.global_bounds)
        pairs = aabb_set.find_intersections(aabb_set)
        print(pairs)


class TestAABBSet(unittest.TestCase):
    def setUp(self):
        pass

    def testAABBSet(self):
        pass


class TestTrimesh(unittest.TestCase):
    def setUp(self):
        pass

    def testAABBSet(self):
        pass


if __name__ == '__main__':
    unittest.main()