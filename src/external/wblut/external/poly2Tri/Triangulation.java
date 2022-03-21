// 
// Decompiled by Procyon v0.5.36
// 

package wblut.external.poly2Tri;

import java.util.ArrayList;

public class Triangulation {
	public static boolean debug;
	public static String debugFileName;

	static {
		Triangulation.debug = false;
		Triangulation.debugFileName = "polygon_triangulation_log.txt";
	}

	public static ArrayList<ArrayList<Integer>> triangulate(final int numContures, final int[] numVerticesInContures,
			final double[][] vertices) {
		final Polygon p = new Polygon(numContures, numVerticesInContures, vertices);
		if (Triangulation.debug) {
			p.setDebugFile(Triangulation.debugFileName);
			p.setDebugOption(Triangulation.debug);
		} else {
			p.setDebugOption(false);
		}
		if (p.triangulation()) {
			return p.triangles();
		}
		return null;
	}
}
