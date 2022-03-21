// 
// Decompiled by Procyon v0.5.36
// 

package wblut.geom;

import java.util.ArrayList;

import wblut.external.poly2Tri.Triangulation;

public class WB_Poly2Tri {
	public static int[] triangulatePolygon(final WB_Polygon poly) {
		final int noc = poly.getNumberOfContours();
		final int[] nopc = poly.getNumberOfPointsPerContour();
		final double[][] points = new double[poly.getNumberOfPoints()][2];
		final WB_OrthoProject op = new WB_OrthoProject(poly.getPlane());
		final WB_Point p = new WB_Point();
		for (int i = 0; i < poly.getNumberOfPoints(); ++i) {
			op.mapPoint3D(poly.getPoint(i), p);
			points[i][0] = p.xd();
			points[i][1] = p.yd();
		}
		final ArrayList<ArrayList<Integer>> triangulation = Triangulation.triangulate(noc, nopc, points);
		final int[] triangles = new int[triangulation.size() * 3];
		int id = 0;
		for (final ArrayList<Integer> triangle : triangulation) {
			triangles[id++] = triangle.get(0);
			triangles[id++] = triangle.get(1);
			triangles[id++] = triangle.get(2);
		}
		return triangles;
	}
}
