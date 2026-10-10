"""Analytic geometry checks only; no native solver or hydrodynamic evidence."""
import math
import sys
from pathlib import Path
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from sphere_geometry import disk_rectangle, sphere_box, segment_sphere


class SphereGeometryTests(unittest.TestCase):
    def test_disk_full_quarter_and_strip(self):
        self.assertAlmostEqual(disk_rectangle(1,[-2,-2],[2,2]),math.pi)
        self.assertAlmostEqual(disk_rectangle(1,[0,0],[1,1]),math.pi/4)
        self.assertAlmostEqual(disk_rectangle(1,[-.5,-2],[.5,2]),
                               math.sqrt(.75)+math.pi/3)

    def test_full_sphere_and_octant(self):
        whole = sphere_box([0,0,0],1,[-1,-1,-1],[1,1,1])
        octant = sphere_box([0,0,0],1,[0,0,0],[1,1,1])
        self.assertAlmostEqual(whole['volume'],4*math.pi/3,places=7)
        self.assertAlmostEqual(octant['volume'],math.pi/6,places=7)
        self.assertEqual(whole['normal_area'],[0,0,0])
        for x in octant['normal_area']:
            self.assertAlmostEqual(x,math.pi/4)

    def test_classification_and_scaling(self):
        inside = sphere_box([0,0,0],2,[-.1]*3,[.1]*3)
        outside = sphere_box([0,0,0],1,[1]*3,[2]*3)
        self.assertEqual(inside['classification'],'inside')
        self.assertEqual(inside['solid_fraction'],1)
        self.assertEqual(outside['volume'],0)
        reference = sphere_box([0,0,0],1,[-.5,-.3,-.8],[.8,.9,.2])
        scaled = sphere_box([3,3,3],2,[2,2.4,1.4],[4.6,4.8,3.4])
        self.assertAlmostEqual(scaled['volume']/8,reference['volume'],places=9)
        for a,b in zip(scaled['solid_face_areas'],reference['solid_face_areas']):
            self.assertAlmostEqual(a/4,b,places=12)

    def test_subdivision_and_shared_face(self):
        whole = sphere_box([.13,-.07,.11],1,[-1.2]*3,[1.2]*3)
        left = sphere_box([.13,-.07,.11],1,[-1.2]*3,[.2,1.2,1.2])
        right = sphere_box([.13,-.07,.11],1,[.2,-1.2,-1.2],[1.2]*3)
        self.assertAlmostEqual(left['volume']+right['volume'],whole['volume'],places=7)
        self.assertEqual(left['solid_face_areas'][1],right['solid_face_areas'][0])
        for a,b,c in zip(left['normal_area'],right['normal_area'],whole['normal_area']):
            self.assertAlmostEqual(a+b,c,places=12)

    def test_spherical_cap_and_axis_permutation(self):
        cap = sphere_box([0]*3,1,[.3,-1,-1],[1,1,1])
        self.assertAlmostEqual(cap['volume'],math.pi*(2/3-.3+.3**3/3),places=7)
        rotated = sphere_box([0]*3,1,[-1,.3,-1],[1,1,1])
        self.assertAlmostEqual(cap['volume'],rotated['volume'],places=7)

    def test_crossing_tangent_and_inside_segment(self):
        hits = segment_sphere([-2,0,0],[2,0,0],[0]*3,1)
        self.assertEqual([h['t'] for h in hits],[.25,.75])
        self.assertEqual([h['normal'] for h in hits],[[-1,0,0],[1,0,0]])
        self.assertEqual(len(segment_sphere([-2,1,0],[2,1,0],[0]*3,1)),1)
        self.assertEqual(segment_sphere([-.2,0,0],[.2,0,0],[0]*3,1),[])

    def test_invalid_geometry(self):
        for radius in (0,-1,float('nan')):
            with self.assertRaises(ValueError):
                sphere_box([0]*3,radius,[-1]*3,[1]*3)
        with self.assertRaises(ValueError):
            segment_sphere([0]*3,[0]*3,[0]*3,1)


if __name__ == '__main__':
    unittest.main()
