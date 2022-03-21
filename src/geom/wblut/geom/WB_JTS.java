package wblut.geom;

import java.awt.Shape;
import java.awt.geom.AffineTransform;
import java.awt.geom.PathIterator;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.eclipse.collections.impl.list.mutable.FastList;
import org.eclipse.collections.impl.map.mutable.UnifiedMap;
import org.locationtech.jts.algorithm.CGAlgorithms;
import org.locationtech.jts.algorithm.ConvexHull;
import org.locationtech.jts.densify.Densifier;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.CoordinateList;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryCollection;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.MultiLineString;
import org.locationtech.jts.geom.MultiPolygon;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.PrecisionModel;
import org.locationtech.jts.operation.buffer.BufferOp;
import org.locationtech.jts.operation.buffer.BufferParameters;
import org.locationtech.jts.simplify.TopologyPreservingSimplifier;
import org.locationtech.jts.triangulate.ConformingDelaunayTriangulationBuilder;
import org.locationtech.jts.triangulate.DelaunayTriangulationBuilder;
import org.locationtech.jts.triangulate.quadedge.QuadEdgeSubdivision;

import wblut.core.WB_ProgressReporter;
import wblut.hemesh.HEC_FromFacelist;
import wblut.hemesh.HE_Face;
import wblut.hemesh.HE_Halfedge;
import wblut.hemesh.HE_Mesh;
import wblut.hemesh.HE_MeshOp;
import wblut.hemesh.HE_Path;
import wblut.hemesh.HE_Vertex;
import wblut.math.WB_Epsilon;

public class WB_JTS {
	private static final GeometryFactory JTSgf;
	private static final WB_Map2D XY;

	static {
		JTSgf = new GeometryFactory(new PrecisionModel(WB_Epsilon.SCALE));
		XY = new WB_GeometryFactory().createEmbeddedPlane();
	}

	private WB_JTS() {
	}

	static Polygon toJTSPolygon2D(final WB_Polygon poly) {
		final int[] npc = poly.getNumberOfPointsPerContour();
		Coordinate[] coords = new Coordinate[npc[0] + 1];
		int i;
		for (i = 0, i = 0; i < npc[0]; ++i) {
			coords[i] = toJTSCoordinate2D(poly.getPoint(i), i);
		}
		coords[i] = toJTSCoordinate2D(poly.getPoint(0), 0);
		final LinearRing shell = WB_JTS.JTSgf.createLinearRing(coords);
		final LinearRing[] holes = new LinearRing[poly.getNumberOfHoles()];
		int index = poly.getNumberOfShellPoints();
		for (i = 0; i < poly.getNumberOfHoles(); ++i) {
			coords = new Coordinate[npc[i + 1] + 1];
			coords[npc[i + 1]] = toJTSCoordinate2D(poly.getPoint(index), index);
			for (int j = 0; j < npc[i + 1]; ++j) {
				coords[j] = toJTSCoordinate2D(poly.getPoint(index), index);
				++index;
			}
			holes[i] = WB_JTS.JTSgf.createLinearRing(coords);
		}
		return WB_JTS.JTSgf.createPolygon(shell, holes);
	}

	static Geometry toJTSMultiPolygon2D(final List<WB_Polygon> polys) {
		final Polygon[] JTSpolys = new Polygon[polys.size()];
		for (int j = 0; j < polys.size(); ++j) {
			final WB_Polygon poly = polys.get(j);
			final int[] npc = poly.getNumberOfPointsPerContour();
			Coordinate[] coords = new Coordinate[npc[0] + 1];
			int i;
			for (i = 0, i = 0; i < npc[0]; ++i) {
				coords[i] = toJTSCoordinate2D(poly.getPoint(i), i);
			}
			coords[i] = toJTSCoordinate2D(poly.getPoint(0), 0);
			final LinearRing shell = WB_JTS.JTSgf.createLinearRing(coords);
			final LinearRing[] holes = new LinearRing[poly.getNumberOfHoles()];
			int index = poly.getNumberOfShellPoints();
			for (i = 0; i < poly.getNumberOfHoles(); ++i) {
				coords = new Coordinate[npc[i + 1] + 1];
				coords[npc[i + 1]] = toJTSCoordinate2D(poly.getPoint(index), index);
				for (int k = 0; k < npc[i + 1]; ++k) {
					coords[k] = toJTSCoordinate2D(poly.getPoint(index), index);
					++index;
				}
				holes[i] = WB_JTS.JTSgf.createLinearRing(coords);
			}
			JTSpolys[j] = WB_JTS.JTSgf.createPolygon(shell, holes);
		}
		return WB_JTS.JTSgf.createMultiPolygon(JTSpolys).buffer(0.0);
	}

	static Coordinate toJTSCoordinate2D(final WB_Coord point, final int i) {
		return new Coordinate(point.xd(), point.yd(), i);
	}

	static WB_Point createPoint2D(final Coordinate coord) {
		return new WB_Point(coord.x, coord.y);
	}

	static WB_Polygon createPolygonFromJTSPolygon2D(final Polygon JTSpoly) {
		final LineString shell = JTSpoly.getExteriorRing();
		Coordinate[] coords = shell.getCoordinates();
		final WB_Coord[] points = new WB_Coord[coords.length - 1];
		for (int i = 0; i < coords.length - 1; ++i) {
			points[i] = createPoint2D(coords[i]);
		}
		final int numholes = JTSpoly.getNumInteriorRing();
		if (numholes > 0) {
			final WB_Coord[][] holecoords = new WB_Coord[numholes][];
			for (int j = 0; j < numholes; ++j) {
				final LineString hole = JTSpoly.getInteriorRingN(j);
				coords = hole.getCoordinates();
				holecoords[j] = new WB_Coord[coords.length - 1];
				for (int k = 0; k < coords.length - 1; ++k) {
					holecoords[j][k] = createPoint2D(coords[k]);
				}
			}
			return new WB_GeometryFactory2D().createPolygonWithHoles(points, holecoords);
		}
		return new WB_GeometryFactory2D().createSimplePolygon(points);
	}

	static List<WB_Polygon> createPolygonsFromJTSGeometry2D(final Geometry geometry) {
		final List<WB_Polygon> polygons = new FastList<WB_Polygon>();
		for (int i = 0; i < geometry.getNumGeometries(); ++i) {
			final Geometry geo = geometry.getGeometryN(i);
			if (!geo.isEmpty()) {
				if (geo.getGeometryType().equals("Polygon")) {
					polygons.add(createPolygonFromJTSPolygon2D((Polygon) geo));
				} else if (geo.getGeometryType().equals("MultiPolygon")) {
					for (int j = 0; j < geo.getNumGeometries(); ++j) {
						final Geometry ggeo = geo.getGeometryN(j);
						polygons.add(createPolygonFromJTSPolygon2D((Polygon) ggeo));
					}
				} else if (geo.getGeometryType().equals("GeometryCollection")) {
					for (int j = 0; j < geo.getNumGeometries(); ++j) {
						final Geometry ggeo = geo.getGeometryN(j);
						polygons.addAll(createPolygonsFromJTSGeometry2D(ggeo));
					}
				}
			}
		}
		return polygons;
	}

	static List<WB_Polygon> createBufferedPolygons2D(final WB_Polygon poly, final double d) {
		final Polygon JTSpoly = toJTSPolygon2D(poly);
		final Geometry result = BufferOp.bufferOp(JTSpoly, d);
		return createPolygonsFromJTSGeometry2D(result);
	}

	static List<WB_Polygon> createBufferedPolygons2D(final Collection<? extends WB_Polygon> poly, final double d) {
		final Polygon[] allPoly = new Polygon[poly.size()];
		int i = 0;
		for (final WB_Polygon pol : poly) {
			allPoly[i++] = toJTSPolygon2D(pol);
		}
		final MultiPolygon collPoly = WB_JTS.JTSgf.createMultiPolygon(allPoly);
		final Geometry result = BufferOp.bufferOp(collPoly, d);
		return createPolygonsFromJTSGeometry2D(result);
	}

	static List<WB_Polygon> createBufferedPolygons2D(final WB_Polygon poly, final double d, final int n) {
		final Polygon JTSpoly = toJTSPolygon2D(poly);
		final BufferParameters parameters = new BufferParameters(n, 1, (n == 0) ? 2 : 1, 5.0);
		final Geometry result = BufferOp.bufferOp(JTSpoly, d, parameters);
		return createPolygonsFromJTSGeometry2D(result);
	}

	static List<WB_Polygon> createBufferedPolygonsStraight2D(final WB_Polygon poly, final double d) {
		final Polygon JTSpoly = toJTSPolygon2D(poly);
		final BufferParameters parameters = new BufferParameters(0, 1, 2, 5.0);
		final Geometry result = BufferOp.bufferOp(JTSpoly, d, parameters);
		return createPolygonsFromJTSGeometry2D(result);
	}

	static List<WB_Polygon> createBufferedPolygons2D(final Collection<? extends WB_Polygon> poly, final double d, final int n) {
		final Polygon[] allPoly = new Polygon[poly.size()];
		int i = 0;
		for (final WB_Polygon pol : poly) {
			allPoly[i++] = toJTSPolygon2D(pol);
		}
		final MultiPolygon collPoly = WB_JTS.JTSgf.createMultiPolygon(allPoly);
		final BufferParameters parameters = new BufferParameters(n, 1, (n == 0) ? 2 : 1, 5.0);
		final Geometry result = BufferOp.bufferOp(collPoly, d, parameters);
		return createPolygonsFromJTSGeometry2D(result);
	}

	static List<WB_Polygon> createBufferedPolygonsStraight2D(final Collection<? extends WB_Polygon> poly, final double d) {
		final Polygon[] allPoly = new Polygon[poly.size()];
		int i = 0;
		for (final WB_Polygon pol : poly) {
			allPoly[i++] = toJTSPolygon2D(pol);
		}
		final MultiPolygon collPoly = WB_JTS.JTSgf.createMultiPolygon(allPoly);
		final BufferParameters parameters = new BufferParameters(0, 1, 2, 5.0);
		final Geometry result = BufferOp.bufferOp(collPoly, d, parameters);
		return createPolygonsFromJTSGeometry2D(result);
	}

	static List<WB_Polygon> createBoundaryPolygons2D(final WB_Polygon poly) {
		final Polygon JTSpoly = toJTSPolygon2D(poly);
		final LineString result = JTSpoly.getExteriorRing();
		return createPolygonsFromJTSGeometry2D(WB_JTS.JTSgf.createPolygon(result.getCoordinates()));
	}

	static List<WB_Polygon> createRibbonPolygons2D(final WB_Polygon poly, final double d) {
		final Polygon JTSpoly = toJTSPolygon2D(poly);
		final Geometry clean = BufferOp.bufferOp(JTSpoly, 0.0);
		final Geometry outer = BufferOp.bufferOp(clean, d * 0.5);
		final Geometry inner = BufferOp.bufferOp(clean, -d * 0.5);
		final Geometry result = outer.difference(inner);
		return createPolygonsFromJTSGeometry2D(result);
	}

	static List<WB_Polygon> createRibbonPolygons2D(final Collection<? extends WB_Polygon> poly, final double d) {
		final Polygon[] allPoly = new Polygon[poly.size()];
		int i = 0;
		for (final WB_Polygon pol : poly) {
			allPoly[i++] = toJTSPolygon2D(pol);
		}
		final MultiPolygon collPoly = WB_JTS.JTSgf.createMultiPolygon(allPoly);
		final Geometry clean = BufferOp.bufferOp(collPoly, 0.0);
		final Geometry outer = BufferOp.bufferOp(clean, d * 0.5);
		final Geometry inner = BufferOp.bufferOp(clean, -d * 0.5);
		final Geometry result = outer.difference(inner);
		return createPolygonsFromJTSGeometry2D(result);
	}

	static List<WB_Polygon> createRibbonPolygons2D(final WB_Polygon poly, final double o, final double i) {
		final Polygon JTSpoly = toJTSPolygon2D(poly);
		final Geometry clean = BufferOp.bufferOp(JTSpoly, 0.0);
		final Geometry outer = BufferOp.bufferOp(clean, o);
		final Geometry inner = BufferOp.bufferOp(clean, -i);
		final Geometry result = outer.difference(inner);
		return createPolygonsFromJTSGeometry2D(result);
	}

	static List<WB_Polygon> createRibbonPolygons2D(final Collection<? extends WB_Polygon> poly, final double o, final double i) {
		final Polygon[] allPoly = new Polygon[poly.size()];
		int j = 0;
		for (final WB_Polygon pol : poly) {
			allPoly[j++] = toJTSPolygon2D(pol);
		}
		final MultiPolygon collPoly = WB_JTS.JTSgf.createMultiPolygon(allPoly);
		final Geometry clean = BufferOp.bufferOp(collPoly, 0.0);
		final Geometry outer = BufferOp.bufferOp(clean, o);
		final Geometry inner = BufferOp.bufferOp(clean, -i);
		final Geometry result = outer.difference(inner);
		return createPolygonsFromJTSGeometry2D(result);
	}

	static List<WB_Polygon> createSimplifiedPolygon2D(final WB_Polygon poly, final double tol) {
		final Polygon JTSpoly = toJTSPolygon2D(poly);
		final Geometry result = TopologyPreservingSimplifier.simplify(JTSpoly, tol);
		return createPolygonsFromJTSGeometry2D(result);
	}

	static List<WB_Polygon> createDensifiedPolygon2D(final WB_Polygon poly, final double max) {
		final Polygon JTSpoly = toJTSPolygon2D(poly);
		final Geometry result = Densifier.densify(JTSpoly, max);
		return createPolygonsFromJTSGeometry2D(result);
	}

	static List<WB_Polygon> unionPolygons2D(final WB_Polygon poly1, final WB_Polygon poly2) {
		final Polygon JTSpoly1 = toJTSPolygon2D(poly1);
		final Polygon JTSpoly2 = toJTSPolygon2D(poly2);
		final Geometry result = JTSpoly1.union(JTSpoly2);
		return createPolygonsFromJTSGeometry2D(result);
	}

	static List<WB_Polygon> unionPolygons2D(final WB_Polygon poly1, final Collection<? extends WB_Polygon> poly2) {
		final Polygon JTSpoly1 = toJTSPolygon2D(poly1);
		final Polygon[] allPoly2 = new Polygon[poly2.size()];
		int i = 0;
		for (final WB_Polygon poly3 : poly2) {
			allPoly2[i++] = toJTSPolygon2D(poly3);
		}
		final MultiPolygon collPoly2 = WB_JTS.JTSgf.createMultiPolygon(allPoly2);
		final Geometry result = JTSpoly1.union(collPoly2);
		return createPolygonsFromJTSGeometry2D(result);
	}

	static List<WB_Polygon> unionPolygons2D(final Collection<? extends WB_Polygon> poly1, final Collection<? extends WB_Polygon> poly2) {
		final Polygon[] allPoly1 = new Polygon[poly1.size()];
		int i = 0;
		for (final WB_Polygon poly3 : poly1) {
			allPoly1[i++] = toJTSPolygon2D(poly3);
		}
		final MultiPolygon collPoly1 = WB_JTS.JTSgf.createMultiPolygon(allPoly1);
		final Polygon[] allPoly2 = new Polygon[poly2.size()];
		i = 0;
		for (final WB_Polygon poly4 : poly2) {
			allPoly2[i++] = toJTSPolygon2D(poly4);
		}
		final MultiPolygon collPoly2 = WB_JTS.JTSgf.createMultiPolygon(allPoly2);
		final Geometry result = collPoly1.union(collPoly2);
		return createPolygonsFromJTSGeometry2D(result);
	}

	static List<WB_Polygon> subtractPolygons2D(final WB_Polygon poly1, final WB_Polygon poly2) {
		final Polygon JTSpoly1 = toJTSPolygon2D(poly1);
		final Polygon JTSpoly2 = toJTSPolygon2D(poly2);
		final Geometry result = JTSpoly1.difference(JTSpoly2);
		return createPolygonsFromJTSGeometry2D(result);
	}

	static List<WB_Polygon> subtractPolygons2D(final WB_Polygon poly1, final Collection<? extends WB_Polygon> poly2) {
		final Polygon JTSpoly1 = toJTSPolygon2D(poly1);
		final Polygon[] allPoly2 = new Polygon[poly2.size()];
		int i = 0;
		for (final WB_Polygon poly3 : poly2) {
			allPoly2[i++] = toJTSPolygon2D(poly3);
		}
		final MultiPolygon collPoly2 = WB_JTS.JTSgf.createMultiPolygon(allPoly2);
		final Geometry result = JTSpoly1.difference(collPoly2);
		return createPolygonsFromJTSGeometry2D(result);
	}

	static List<WB_Polygon> subtractPolygons2D(final Collection<? extends WB_Polygon> poly1, final Collection<? extends WB_Polygon> poly2) {
		final Polygon[] allPoly1 = new Polygon[poly1.size()];
		int i = 0;
		for (final WB_Polygon poly3 : poly1) {
			allPoly1[i++] = toJTSPolygon2D(poly3);
		}
		final MultiPolygon collPoly1 = WB_JTS.JTSgf.createMultiPolygon(allPoly1);
		final Polygon[] allPoly2 = new Polygon[poly2.size()];
		i = 0;
		for (final WB_Polygon poly4 : poly2) {
			allPoly2[i++] = toJTSPolygon2D(poly4);
		}
		final MultiPolygon collPoly2 = WB_JTS.JTSgf.createMultiPolygon(allPoly2);
		final Geometry result = collPoly1.difference(collPoly2);
		return createPolygonsFromJTSGeometry2D(result);
	}

	static List<WB_Polygon> subtractPolygons2D(final Collection<? extends WB_Polygon> poly1, final WB_Polygon poly2) {
		final Polygon[] allPoly1 = new Polygon[poly1.size()];
		int i = 0;
		for (final WB_Polygon poly3 : poly1) {
			allPoly1[i++] = toJTSPolygon2D(poly3);
		}
		final MultiPolygon collPoly1 = WB_JTS.JTSgf.createMultiPolygon(allPoly1);
		final Polygon JTSpoly2 = toJTSPolygon2D(poly2);
		final Geometry result = collPoly1.difference(JTSpoly2);
		return createPolygonsFromJTSGeometry2D(result);
	}

	static List<WB_Polygon> intersectPolygons2D(final WB_Polygon poly1, final WB_Polygon poly2) {
		final Polygon JTSpoly1 = toJTSPolygon2D(poly1);
		final Polygon JTSpoly2 = toJTSPolygon2D(poly2);
		final Geometry result = JTSpoly1.intersection(JTSpoly2);
		return createPolygonsFromJTSGeometry2D(result);
	}

	static List<WB_Polygon> intersectPolygons2D(final WB_Polygon poly1, final Collection<? extends WB_Polygon> poly2) {
		final Polygon JTSpoly1 = toJTSPolygon2D(poly1);
		final Polygon[] allPoly2 = new Polygon[poly2.size()];
		int i = 0;
		for (final WB_Polygon poly3 : poly2) {
			allPoly2[i++] = toJTSPolygon2D(poly3);
		}
		final MultiPolygon collPoly2 = WB_JTS.JTSgf.createMultiPolygon(allPoly2);
		final Geometry result = JTSpoly1.intersection(collPoly2);
		return createPolygonsFromJTSGeometry2D(result);
	}

	static List<WB_Polygon> intersectPolygons2D(final Collection<? extends WB_Polygon> poly1,
			final Collection<? extends WB_Polygon> poly2) {
		final Polygon[] allPoly1 = new Polygon[poly1.size()];
		int i = 0;
		for (final WB_Polygon poly3 : poly1) {
			allPoly1[i++] = toJTSPolygon2D(poly3);
		}
		final MultiPolygon collPoly1 = WB_JTS.JTSgf.createMultiPolygon(allPoly1);
		final Polygon[] allPoly2 = new Polygon[poly2.size()];
		i = 0;
		for (final WB_Polygon poly4 : poly2) {
			allPoly2[i++] = toJTSPolygon2D(poly4);
		}
		final MultiPolygon collPoly2 = WB_JTS.JTSgf.createMultiPolygon(allPoly2);
		final Geometry result = collPoly1.intersection(collPoly2);
		return createPolygonsFromJTSGeometry2D(result);
	}

	static List<WB_Polygon> symDifferencePolygons2D(final WB_Polygon poly1, final WB_Polygon poly2) {
		final Polygon JTSpoly1 = toJTSPolygon2D(poly1);
		final Polygon JTSpoly2 = toJTSPolygon2D(poly2);
		final Geometry result = JTSpoly1.symDifference(JTSpoly2);
		return createPolygonsFromJTSGeometry2D(result);
	}

	static List<WB_Polygon> symDifferencePolygons2D(final WB_Polygon poly1, final Collection<? extends WB_Polygon> poly2) {
		final Polygon JTSpoly1 = toJTSPolygon2D(poly1);
		final Polygon[] allPoly2 = new Polygon[poly2.size()];
		int i = 0;
		for (final WB_Polygon poly3 : poly2) {
			allPoly2[i++] = toJTSPolygon2D(poly3);
		}
		final MultiPolygon collPoly2 = WB_JTS.JTSgf.createMultiPolygon(allPoly2);
		final Geometry result = JTSpoly1.symDifference(collPoly2);
		return createPolygonsFromJTSGeometry2D(result);
	}

	static List<WB_Polygon> symDifferencePolygons2D(final Collection<? extends WB_Polygon> poly1,
			final Collection<? extends WB_Polygon> poly2) {
		final Polygon[] allPoly1 = new Polygon[poly1.size()];
		int i = 0;
		for (final WB_Polygon poly3 : poly1) {
			allPoly1[i++] = toJTSPolygon2D(poly3);
		}
		final MultiPolygon collPoly1 = WB_JTS.JTSgf.createMultiPolygon(allPoly1);
		final Polygon[] allPoly2 = new Polygon[poly2.size()];
		i = 0;
		for (final WB_Polygon poly4 : poly2) {
			allPoly2[i++] = toJTSPolygon2D(poly4);
		}
		final MultiPolygon collPoly2 = WB_JTS.JTSgf.createMultiPolygon(allPoly2);
		final Geometry result = collPoly1.symDifference(collPoly2);
		return createPolygonsFromJTSGeometry2D(result);
	}

	static List<WB_Polygon> constrainPolygons2D(final WB_Polygon poly, final WB_Polygon container) {
		final Polygon JTSpoly1 = toJTSPolygon2D(poly);
		final Polygon JTSpoly2 = toJTSPolygon2D(container);
		final Geometry result = JTSpoly1.intersection(JTSpoly2);
		return createPolygonsFromJTSGeometry2D(result);
	}

	static List<WB_Polygon> constrainPolygons2D(final WB_Polygon[] polygons, final WB_Polygon container) {
		final List<WB_Polygon> polys = new FastList<WB_Polygon>();
		for (final WB_Polygon poly : polygons) {
			final Polygon JTSpoly1 = toJTSPolygon2D(poly);
			final Polygon JTSpoly2 = toJTSPolygon2D(container);
			final Geometry result = JTSpoly1.intersection(JTSpoly2);
			polys.addAll(createPolygonsFromJTSGeometry2D(result));
		}
		return polys;
	}

	static List<WB_Polygon> constrainPolygons2D(final List<WB_Polygon> polygons, final WB_Polygon container) {
		final List<WB_Polygon> polys = new FastList<WB_Polygon>();
		for (final WB_Polygon poly : polygons) {
			final Polygon JTSpoly1 = toJTSPolygon2D(poly);
			final Polygon JTSpoly2 = toJTSPolygon2D(container);
			final Geometry result = JTSpoly1.intersection(JTSpoly2);
			if (!result.isEmpty()) {
				polys.addAll(createPolygonsFromJTSGeometry2D(result));
			}
		}
		return polys;
	}

	static WB_Polygon createPolygonConvexHull2D(final WB_Polygon poly) {
		final Polygon JTSpoly = toJTSPolygon2D(poly);
		final Geometry result = new ConvexHull(JTSpoly).getConvexHull();
		if (result.getGeometryType().equals("Polygon")) {
			return createPolygonFromJTSPolygon2D((Polygon) result);
		}
		return null;
	}

	static WB_Polygon createPolygonFromJTSPolygon(final Polygon JTSpoly, final WB_Map2D map) {
		final LineString shell = JTSpoly.getExteriorRing();
		Coordinate[] coords = shell.getCoordinates();
		final WB_Point[] points = new WB_Point[coords.length - 1];
		for (int i = 0; i < coords.length - 1; ++i) {
			map.unmapPoint3D((points[i] = createPoint2D(coords[i])), points[i]);
		}
		final int numholes = JTSpoly.getNumInteriorRing();
		if (numholes > 0) {
			final WB_Point[][] holecoords = new WB_Point[numholes][];
			for (int j = 0; j < numholes; ++j) {
				final LineString hole = JTSpoly.getInteriorRingN(j);
				coords = hole.getCoordinates();
				holecoords[j] = new WB_Point[coords.length - 1];
				for (int k = 0; k < coords.length - 1; ++k) {
					map.unmapPoint3D((holecoords[j][k] = createPoint2D(coords[k])), holecoords[j][k]);
				}
			}
			return new WB_GeometryFactory2D().createPolygonWithHoles(points, holecoords);
		}
		return new WB_GeometryFactory2D().createSimplePolygon(points);
	}

	static Polygon toJTSPolygon(final WB_Polygon poly, final WB_Map2D map) {
		final int[] npc = poly.getNumberOfPointsPerContour();
		Coordinate[] coords = new Coordinate[npc[0] + 1];
		int i;
		for (i = 0, i = 0; i < npc[0]; ++i) {
			coords[i] = toJTSCoordinate(poly.getPoint(i), i, map);
		}
		coords[i] = toJTSCoordinate(poly.getPoint(0), 0, map);
		final LinearRing shell = WB_JTS.JTSgf.createLinearRing(coords);
		final LinearRing[] holes = new LinearRing[poly.getNumberOfHoles()];
		int index = poly.getNumberOfShellPoints();
		for (i = 0; i < poly.getNumberOfHoles(); ++i) {
			coords = new Coordinate[npc[i + 1] + 1];
			coords[npc[i + 1]] = toJTSCoordinate(poly.getPoint(index), index, map);
			for (int j = 0; j < npc[i + 1]; ++j) {
				coords[j] = toJTSCoordinate(poly.getPoint(index), index, map);
				++index;
			}
			holes[i] = WB_JTS.JTSgf.createLinearRing(coords);
		}
		return WB_JTS.JTSgf.createPolygon(shell, holes);
	}

	static Geometry toJTSMultiPolygon(final List<WB_Polygon> polys, final WB_Map2D map) {
		final Polygon[] JTSpolys = new Polygon[polys.size()];
		for (int j = 0; j < polys.size(); ++j) {
			final WB_Polygon poly = polys.get(j);
			final int[] npc = poly.getNumberOfPointsPerContour();
			Coordinate[] coords = new Coordinate[npc[0] + 1];
			int i;
			for (i = 0, i = 0; i < npc[0]; ++i) {
				coords[i] = toJTSCoordinate(poly.getPoint(i), i, map);
			}
			coords[i] = toJTSCoordinate2D(poly.getPoint(0), 0);
			final LinearRing shell = WB_JTS.JTSgf.createLinearRing(coords);
			final LinearRing[] holes = new LinearRing[poly.getNumberOfHoles()];
			int index = poly.getNumberOfShellPoints();
			for (i = 0; i < poly.getNumberOfHoles(); ++i) {
				coords = new Coordinate[npc[i + 1] + 1];
				coords[npc[i + 1]] = toJTSCoordinate(poly.getPoint(index), index, map);
				for (int k = 0; k < npc[i + 1]; ++k) {
					coords[k] = toJTSCoordinate(poly.getPoint(index), index, map);
					++index;
				}
				holes[i] = WB_JTS.JTSgf.createLinearRing(coords);
			}
			JTSpolys[j] = WB_JTS.JTSgf.createPolygon(shell, holes);
		}
		return WB_JTS.JTSgf.createMultiPolygon(JTSpolys).buffer(0.0);
	}

	static Coordinate toJTSCoordinate(final WB_Coord point, final int i, final WB_Map2D map) {
		final WB_Point mp = new WB_Point();
		map.mapPoint3D(point, mp);
		return new Coordinate(mp.xd(), mp.yd(), i);
	}

	private static List<WB_Polygon> createPolygonsFromJTSGeometry(final Geometry geometry, final WB_Map2D map) {
		final List<WB_Polygon> polygons = new FastList<WB_Polygon>();
		for (int i = 0; i < geometry.getNumGeometries(); ++i) {
			final Geometry geo = geometry.getGeometryN(i);
			if (!geo.isEmpty()) {
				if (geo.getGeometryType().equals("Polygon")) {
					polygons.add(createPolygonFromJTSPolygon((Polygon) geo, map));
				} else if (geo.getGeometryType().equals("MultiPolygon")) {
					for (int j = 0; j < geo.getNumGeometries(); ++j) {
						final Geometry ggeo = geo.getGeometryN(j);
						polygons.add(createPolygonFromJTSPolygon((Polygon) ggeo, map));
					}
				} else if (geo.getGeometryType().equals("GeometryCollection")) {
					for (int j = 0; j < geo.getNumGeometries(); ++j) {
						final Geometry ggeo = geo.getGeometryN(j);
						polygons.addAll(createPolygonsFromJTSGeometry(ggeo, map));
					}
				}
			}
		}
		return polygons;
	}

	static WB_Polygon createPolygonConvexHull(final WB_Polygon poly) {
		final WB_Map2D map = new WB_PlanarMap(poly.getPlane(0.0));
		final Polygon JTSpoly = toJTSPolygon(poly, map);
		final Geometry result = new ConvexHull(JTSpoly).getConvexHull();
		if (result.getGeometryType().equals("Polygon")) {
			createPolygonFromJTSPolygon((Polygon) result, map);
		}
		return null;
	}

	static List<WB_Polygon> createBufferedPolygons(final WB_Polygon poly, final double d) {
		final WB_Map2D map = new WB_PlanarMap(poly.getPlane(0.0));
		final Polygon JTSpoly = toJTSPolygon(poly, map);
		final Geometry result = BufferOp.bufferOp(JTSpoly, d);
		return createPolygonsFromJTSGeometry(result, map);
	}

	static List<WB_Polygon> createBufferedPolygons(final WB_Polygon poly, final double d, final int n) {
		final WB_Map2D map = new WB_PlanarMap(poly.getPlane(0.0));
		final Polygon JTSpoly = toJTSPolygon(poly, map);
		final BufferParameters parameters = new BufferParameters(n, 1, (n == 0) ? 2 : 1, 5.0);
		final Geometry result = BufferOp.bufferOp(JTSpoly, d, parameters);
		return createPolygonsFromJTSGeometry(result, map);
	}

	static List<WB_Polygon> createBufferedPolygonsStraight(final WB_Polygon poly, final double d) {
		final WB_Map2D map = new WB_PlanarMap(poly.getPlane(0.0));
		final Polygon JTSpoly = toJTSPolygon(poly, map);
		final BufferParameters parameters = new BufferParameters(0, 1, 2, 5.0);
		final Geometry result = BufferOp.bufferOp(JTSpoly, d, parameters);
		return createPolygonsFromJTSGeometry(result, map);
	}

	static List<WB_Polygon> createBufferedPolygons(final Collection<? extends WB_Polygon> poly, final double d) {
		final WB_Map2D map = new WB_PlanarMap(((WB_Polygon) poly.iterator().next()).getPlane(0.0));
		final Polygon[] allPoly = new Polygon[poly.size()];
		int i = 0;
		for (final WB_Polygon pol : poly) {
			allPoly[i++] = toJTSPolygon(pol, map);
		}
		final MultiPolygon collPoly = WB_JTS.JTSgf.createMultiPolygon(allPoly);
		final Geometry result = BufferOp.bufferOp(collPoly, d);
		return createPolygonsFromJTSGeometry(result, map);
	}

	static List<WB_Polygon> createBufferedPolygons(final Collection<? extends WB_Polygon> poly, final double d, final int n) {
		final WB_Map2D map = new WB_PlanarMap(((WB_Polygon) poly.iterator().next()).getPlane(0.0));
		final Polygon[] allPoly = new Polygon[poly.size()];
		int i = 0;
		for (final WB_Polygon pol : poly) {
			allPoly[i++] = toJTSPolygon(pol, map);
		}
		final MultiPolygon collPoly = WB_JTS.JTSgf.createMultiPolygon(allPoly);
		final BufferParameters parameters = new BufferParameters(n, 1, (n == 0) ? 2 : 1, 5.0);
		final Geometry result = BufferOp.bufferOp(collPoly, d, parameters);
		return createPolygonsFromJTSGeometry(result, map);
	}

	static List<WB_Polygon> createBufferedPolygonsStraight(final Collection<? extends WB_Polygon> poly, final double d) {
		final WB_Map2D map = new WB_PlanarMap(((WB_Polygon) poly.iterator().next()).getPlane(0.0));
		final Polygon[] allPoly = new Polygon[poly.size()];
		int i = 0;
		for (final WB_Polygon pol : poly) {
			allPoly[i++] = toJTSPolygon(pol, map);
		}
		final MultiPolygon collPoly = WB_JTS.JTSgf.createMultiPolygon(allPoly);
		final BufferParameters parameters = new BufferParameters(0, 1, 2, 5.0);
		final Geometry result = BufferOp.bufferOp(collPoly, d, parameters);
		return createPolygonsFromJTSGeometry(result, map);
	}

	static List<WB_Polygon> createBoundaryPolygons(final WB_Polygon poly) {
		final WB_Map2D map = new WB_PlanarMap(poly.getPlane(0.0));
		final Polygon JTSpoly = toJTSPolygon(poly, map);
		final LineString result = JTSpoly.getExteriorRing();
		return createPolygonsFromJTSGeometry(WB_JTS.JTSgf.createPolygon(result.getCoordinates()), map);
	}

	static List<WB_Polygon> createRibbonPolygons(final WB_Polygon poly, final double d) {
		final WB_Map2D map = new WB_PlanarMap(poly.getPlane(0.0));
		final Polygon JTSpoly = toJTSPolygon(poly, map);
		final Geometry outer = BufferOp.bufferOp(JTSpoly, d * 0.5);
		final Geometry inner = BufferOp.bufferOp(JTSpoly, -d * 0.5);
		final Geometry result = outer.difference(inner);
		return createPolygonsFromJTSGeometry(result, map);
	}

	static List<WB_Polygon> createRibbonPolygons(final Collection<? extends WB_Polygon> poly, final double d) {
		final WB_Map2D map = new WB_PlanarMap(((WB_Polygon) poly.iterator().next()).getPlane(0.0));
		final Polygon[] allPoly = new Polygon[poly.size()];
		int i = 0;
		for (final WB_Polygon pol : poly) {
			allPoly[i++] = toJTSPolygon(pol, map);
		}
		final MultiPolygon collPoly = WB_JTS.JTSgf.createMultiPolygon(allPoly);
		final Geometry outer = BufferOp.bufferOp(collPoly, d * 0.5);
		final Geometry inner = BufferOp.bufferOp(collPoly, -d * 0.5);
		final Geometry result = outer.difference(inner);
		return createPolygonsFromJTSGeometry(result, map);
	}

	static List<WB_Polygon> createRibbonPolygons(final WB_Polygon poly, final double o, final double i) {
		final WB_Map2D map = new WB_PlanarMap(poly.getPlane(0.0));
		final Polygon JTSpoly = toJTSPolygon(poly, map);
		final Geometry outer = BufferOp.bufferOp(JTSpoly, o);
		final Geometry inner = BufferOp.bufferOp(JTSpoly, -i);
		final Geometry result = outer.difference(inner);
		return createPolygonsFromJTSGeometry(result, map);
	}

	static List<WB_Polygon> createRibbonPolygons(final Collection<? extends WB_Polygon> poly, final double o, final double i) {
		final WB_Map2D map = new WB_PlanarMap(((WB_Polygon) poly.iterator().next()).getPlane(0.0));
		final Polygon[] allPoly = new Polygon[poly.size()];
		int j = 0;
		for (final WB_Polygon pol : poly) {
			allPoly[j++] = toJTSPolygon(pol, map);
		}
		final MultiPolygon collPoly = WB_JTS.JTSgf.createMultiPolygon(allPoly);
		final Geometry outer = BufferOp.bufferOp(collPoly, o);
		final Geometry inner = BufferOp.bufferOp(collPoly, -i);
		final Geometry result = outer.difference(inner);
		return createPolygonsFromJTSGeometry(result, map);
	}

	static List<WB_Polygon> createSimplifiedPolygon(final WB_Polygon poly, final double tol) {
		final WB_Map2D map = new WB_PlanarMap(poly.getPlane(0.0));
		final Polygon JTSpoly = toJTSPolygon(poly, map);
		final Geometry result = TopologyPreservingSimplifier.simplify(JTSpoly, tol);
		return createPolygonsFromJTSGeometry(result, map);
	}

	static List<WB_Polygon> createDensifiedPolygon(final WB_Polygon poly, final double max) {
		final WB_Map2D map = new WB_PlanarMap(poly.getPlane(0.0));
		final Polygon JTSpoly = toJTSPolygon(poly, map);
		final Geometry result = Densifier.densify(JTSpoly, max);
		return createPolygonsFromJTSGeometry(result, map);
	}

	static List<WB_Polygon> constrainPolygons(final WB_Polygon poly, final WB_Polygon container) {
		final WB_Map2D map = new WB_PlanarMap(poly.getPlane(0.0));
		final Polygon JTSpoly1 = toJTSPolygon(poly, map);
		final Polygon JTSpoly2 = toJTSPolygon(container, map);
		final Geometry result = JTSpoly1.intersection(JTSpoly2);
		return createPolygonsFromJTSGeometry(result, map);
	}

	static List<WB_Polygon> constrainPolygons(final WB_Polygon[] polygons, final WB_Polygon container) {
		final WB_Map2D map = new WB_PlanarMap(polygons[0].getPlane(0.0));
		final List<WB_Polygon> polys = new FastList<WB_Polygon>();
		for (final WB_Polygon poly : polygons) {
			final Polygon JTSpoly1 = toJTSPolygon(poly, map);
			final Polygon JTSpoly2 = toJTSPolygon(container, map);
			final Geometry result = JTSpoly1.intersection(JTSpoly2);
			polys.addAll(createPolygonsFromJTSGeometry(result, map));
		}
		return polys;
	}

	static List<WB_Polygon> constrainPolygons(final List<WB_Polygon> polygons, final WB_Polygon container) {
		final WB_Map2D map = new WB_PlanarMap(polygons.get(0).getPlane(0.0));
		final List<WB_Polygon> polys = new FastList<WB_Polygon>();
		for (final WB_Polygon poly : polygons) {
			final Polygon JTSpoly1 = toJTSPolygon(poly, map);
			final Polygon JTSpoly2 = toJTSPolygon(container, map);
			final Geometry result = JTSpoly1.intersection(JTSpoly2);
			if (!result.isEmpty()) {
				polys.addAll(createPolygonsFromJTSGeometry(result, map));
			}
		}
		return polys;
	}

	static WB_Triangulation2D triangulate2D(final WB_CoordCollection points) {
		final int n = points.size();
		final List<Coordinate> coords = new FastList<>();
		for (int i = 0; i < n; ++i) {
			final WB_Coord c = points.get(i);
			coords.add(new Coordinate(c.xd(), c.yd(), i));
		}
		final WB_Triangulation2D result = getTriangles2D(coords);
		return result;
	}

	static WB_Triangulation2D getTriangles2D(final List<Coordinate> coords) {
		final DelaunayTriangulationBuilder dtb = new DelaunayTriangulationBuilder();
		dtb.setSites(coords);
		final QuadEdgeSubdivision qesd = dtb.getSubdivision();
		final GeometryCollection tris = (GeometryCollection) qesd.getTriangles(WB_JTS.JTSgf);
		final int ntris = tris.getNumGeometries();
		List<int[]> result = new FastList<int[]>();
		for (int i = 0; i < ntris; ++i) {
			final Polygon tri = (Polygon) tris.getGeometryN(i);
			final Coordinate[] tricoord = tri.getCoordinates();
			final int[] triind = new int[3];
			for (int j = 0; j < 3; ++j) {
				triind[j] = (int) tricoord[j].z;
			}
			result.add(triind);
		}
		final int[] T = new int[3 * result.size()];
		for (int k = 0; k < result.size(); ++k) {
			T[3 * k] = result.get(k)[0];
			T[3 * k + 1] = result.get(k)[1];
			T[3 * k + 2] = result.get(k)[2];
		}
		final MultiLineString edges = (MultiLineString) qesd.getEdges(WB_JTS.JTSgf);
		final int nedges = edges.getNumGeometries();
		result = new FastList<int[]>();
		for (int l = 0; l < nedges; ++l) {
			final LineString edge = (LineString) edges.getGeometryN(l);
			final Coordinate[] edgecoord = edge.getCoordinates();
			final int[] edgeind = new int[2];
			for (int m = 0; m < 2; ++m) {
				edgeind[m] = (int) edgecoord[m].z;
			}
			result.add(edgeind);
		}
		final int[] E = new int[2 * result.size()];
		for (int i2 = 0; i2 < result.size(); ++i2) {
			E[2 * i2] = result.get(i2)[0];
			E[2 * i2 + 1] = result.get(i2)[1];
		}
		return new WB_Triangulation2D(T, E);
	}

	static WB_Triangulation2DWithPoints getConformingTriangles2D(final Coordinate[] coords, final int[] constraints, final double tol) {
		final int m = constraints.length;
		final LineString[] constraintlines = new LineString[m / 2];
		for (int i = 0; i < m; i += 2) {
			final Coordinate[] pair = { coords[constraints[i]], coords[constraints[i + 1]] };
			constraintlines[i / 2] = WB_JTS.JTSgf.createLineString(pair);
		}
		final ConformingDelaunayTriangulationBuilder dtb = new ConformingDelaunayTriangulationBuilder();
		dtb.setTolerance(tol);
		dtb.setSites(WB_JTS.JTSgf.createMultiPoint(coords));
		dtb.setConstraints(WB_JTS.JTSgf.createMultiLineString(constraintlines));
		final QuadEdgeSubdivision qesd = dtb.getSubdivision();
		final GeometryCollection tris = (GeometryCollection) qesd.getTriangles(WB_JTS.JTSgf);
		final Coordinate[] newcoords = tris.getCoordinates();
		final List<WB_Coord> uniquePoints = new FastList<WB_Coord>();
		final WB_KDTreeInteger3D<WB_Point> tree = new WB_KDTreeInteger3D<WB_Point>();
		int currentSize = 0;
		Coordinate[] array;
		for (int length = (array = newcoords).length, n = 0; n < length; ++n) {
			final Coordinate newcoord = array[n];
			final WB_Point p = new WB_Point(newcoord.x, newcoord.y, 0.0);
			final Integer index = tree.add(p, currentSize);
			if (index == -1) {
				++currentSize;
				uniquePoints.add(p);
			}
		}
		final int ntris = tris.getNumGeometries();
		List<int[]> result = new FastList<int[]>();
		for (int j = 0; j < ntris; ++j) {
			final Polygon tri = (Polygon) tris.getGeometryN(j);
			final Coordinate[] tricoord = tri.getCoordinates();
			final int[] triind = new int[3];
			for (int k = 0; k < 3; ++k) {
				triind[k] = tree.add(new WB_Point(tricoord[k].x, tricoord[k].y, 0.0), 0);
			}
			result.add(triind);
		}
		final int[] T = new int[3 * result.size()];
		for (int l = 0; l < result.size(); ++l) {
			T[3 * l] = result.get(l)[0];
			T[3 * l + 1] = result.get(l)[1];
			T[3 * l + 2] = result.get(l)[2];
		}
		final MultiLineString edges = (MultiLineString) qesd.getEdges(WB_JTS.JTSgf);
		final int nedges = edges.getNumGeometries();
		result = new FastList<int[]>();
		for (int i2 = 0; i2 < nedges; ++i2) {
			final LineString edge = (LineString) edges.getGeometryN(i2);
			final Coordinate[] edgecoord = edge.getCoordinates();
			final int[] edgeind = new int[2];
			for (int j2 = 0; j2 < 2; ++j2) {
				edgeind[j2] = tree.add(new WB_Point(edgecoord[j2].x, edgecoord[j2].y, 0.0), 0);
			}
			result.add(edgeind);
		}
		final int[] E = new int[2 * result.size()];
		for (int i3 = 0; i3 < result.size(); ++i3) {
			E[2 * i3] = result.get(i3)[0];
			E[2 * i3 + 1] = result.get(i3)[1];
		}
		final List<WB_Coord> Points = new FastList<WB_Coord>();
		for (int i4 = 0; i4 < uniquePoints.size(); ++i4) {
			Points.add(uniquePoints.get(i4));
		}
		return new WB_Triangulation2DWithPoints(T, E, Points);
	}

	static WB_Triangulation2DWithPoints triangulateConforming2D(final WB_CoordCollection points) {
		final int n = points.size();
		final int[] constraints = new int[2 * n];
		int i = 0;
		int j = n - 1;
		while (i < n) {
			constraints[2 * i] = j;
			constraints[2 * i + 1] = i;
			j = i++;
		}
		final Coordinate[] coords = new Coordinate[n];
		for (int k = 0; k < n; ++k) {
			final WB_Coord c = points.get(k);
			coords[k] = new Coordinate(c.xd(), c.yd(), k);
		}
		return getConformingTriangles2D(coords, constraints, WB_Epsilon.EPSILON);
	}

	static WB_Triangulation2DWithPoints triangulateConforming2D(final WB_CoordCollection points, final double tolerance) {
		final int n = points.size();
		final int[] constraints = new int[2 * n];
		int i = 0;
		int j = n - 1;
		while (i < n) {
			constraints[2 * i] = j;
			constraints[2 * i + 1] = i;
			j = i++;
		}
		final Coordinate[] coords = new Coordinate[n];
		for (int k = 0; k < n; ++k) {
			final WB_Coord c = points.get(k);
			coords[k] = new Coordinate(c.xd(), c.yd(), k);
		}
		return getConformingTriangles2D(coords, constraints, tolerance);
	}

	static WB_Triangulation2DWithPoints triangulateConforming2D(final WB_CoordCollection points, final int[] constraints) {
		if (constraints == null) {
			return new WB_Triangulation2DWithPoints(WB_TriangulationFactory2D.triangulate2D(points));
		}
		final int m = constraints.length;
		if (m == 0 || m % 2 == 1) {
			return new WB_Triangulation2DWithPoints(WB_TriangulationFactory2D.triangulate2D(points));
		}
		final int n = points.size();
		final Coordinate[] coords = new Coordinate[n];
		for (int i = 0; i < n; ++i) {
			final WB_Coord c = points.get(i);
			coords[i] = new Coordinate(c.xd(), c.yd(), i);
		}
		return getConformingTriangles2D(coords, constraints, WB_Epsilon.EPSILON);
	}

	static WB_Triangulation2DWithPoints triangulateConforming2D(final WB_CoordCollection points, final int[] constraints,
			final double tolerance) {
		if (constraints == null) {
			return new WB_Triangulation2DWithPoints(WB_TriangulationFactory2D.triangulate2D(points));
		}
		final int m = constraints.length;
		if (m == 0 || m % 2 == 1) {
			return new WB_Triangulation2DWithPoints(WB_TriangulationFactory2D.triangulate2D(points));
		}
		final int n = points.size();
		final Coordinate[] coords = new Coordinate[n];
		for (int i = 0; i < n; ++i) {
			final WB_Coord c = points.get(i);
			coords[i] = new Coordinate(c.xd(), c.yd(), i);
		}
		return getConformingTriangles2D(coords, constraints, tolerance);
	}

	static WB_Triangulation2DWithPoints triangulateConforming2D(final Collection<? extends WB_Coord> points, final WB_Map2D context) {
		final int n = points.size();
		final int[] constraints = new int[2 * n];
		int i = 0;
		int j = n - 1;
		while (i < n) {
			constraints[2 * i] = j;
			constraints[2 * i + 1] = i;
			j = i++;
		}
		final Coordinate[] coords = new Coordinate[n];
		final WB_Point point = new WB_Point();
		int k = 0;
		for (final WB_Coord p : points) {
			context.mapPoint3D(p, point);
			coords[k] = new Coordinate(point.xd(), point.yd(), k);
			++k;
		}
		return getConformingTriangles2D(coords, constraints, WB_Epsilon.EPSILON);
	}

	static WB_Triangulation2DWithPoints triangulateConforming2D(final Collection<? extends WB_Coord> points, final double tol,
			final WB_Map2D context) {
		final int n = points.size();
		final int[] constraints = new int[2 * n];
		int i = 0;
		int j = n - 1;
		while (i < n) {
			constraints[2 * i] = j;
			constraints[2 * i + 1] = i;
			j = i++;
		}
		final Coordinate[] coords = new Coordinate[n];
		final WB_Point point = new WB_Point();
		int k = 0;
		for (final WB_Coord p : points) {
			context.mapPoint3D(p, point);
			coords[k] = new Coordinate(point.xd(), point.yd(), k);
			++k;
		}
		return getConformingTriangles2D(coords, constraints, tol);
	}

	static WB_Triangulation2DWithPoints triangulateConforming2D(final WB_Polygon polygon, final double tol) {
		final int n = polygon.getNumberOfPoints();
		final int[] constraints = new int[2 * n];
		int index = 0;
		int i = 0;
		int j = polygon.getNumberOfShellPoints() - 1;
		while (i < polygon.getNumberOfShellPoints()) {
			constraints[2 * index] = j;
			constraints[2 * index + 1] = i;
			++index;
			j = i++;
		}
		final int nh = polygon.getNumberOfHoles();
		final int[] npc = polygon.getNumberOfPointsPerContour();
		int offset = 0;
		for (int k = 0; k < nh; ++k) {
			offset += npc[k];
			for (int l = 0; l < npc[k + 1]; ++l) {
				constraints[2 * index] = offset + l;
				constraints[2 * index + 1] = offset + (l + 1) % npc[k + 1];
				++index;
			}
		}
		final Coordinate[] coords = new Coordinate[n];
		final WB_Point p = new WB_Point();
		final WB_Map2D context = new WB_GeometryFactory().createEmbeddedPlane(polygon.getPlane());
		for (int m = 0; m < n; ++m) {
			context.mapPoint3D(polygon.getPoint(m), p);
			coords[m] = new Coordinate(p.xd(), p.yd(), m);
		}
		final WB_Triangulation2DWithPoints tri = getConformingTriangles2D(coords, constraints, tol);
		final List<WB_Point> upoints = new FastList<WB_Point>();
		final WB_CoordCollection points = tri.getPoints();
		for (int i2 = 0; i2 < points.size(); ++i2) {
			final WB_Point q = new WB_Point();
			context.unmapPoint2D(points.get(i2), q);
			upoints.add(q);
		}
		return new WB_Triangulation2DWithPoints(tri.getTriangles(), tri.getEdges(), upoints);
	}

	private static WB_Voronoi2D voronoi2D(final ArrayList<Coordinate> coords, final WB_Map2D context) {
		final DelaunayTriangulationBuilder dtb = new DelaunayTriangulationBuilder();
		dtb.setSites(coords);
		final QuadEdgeSubdivision qes = dtb.getSubdivision();
		final GeometryCollection polys = (GeometryCollection) qes.getVoronoiDiagram(WB_JTS.JTSgf);
		final int npolys = polys.getNumGeometries();
		final List<WB_VoronoiCell2D> result = new FastList<WB_VoronoiCell2D>();
		for (int i = 0; i < npolys; ++i) {
			final Polygon poly = (Polygon) polys.getGeometryN(i);
			final Coordinate[] polycoord = poly.getCoordinates();
			final List<WB_Coord> polypoints = new FastList<WB_Coord>();
			Coordinate[] array;
			for (int length = (array = polycoord).length, j = 0; j < length; ++j) {
				final Coordinate element = array[j];
				polypoints.add(toPoint(element.x, element.y, context));
			}
			final Point centroid = poly.getCentroid();
			final WB_Point pc = (centroid == null) ? null : toPoint(centroid.getX(), centroid.getY(), context);
			final int index = (int) ((Coordinate) poly.getUserData()).z;
			final double area = poly.getArea();
			result.add(new WB_VoronoiCell2D(polypoints, index, createPoint2D(coords.get(index)), area, pc));
		}
		return new WB_Voronoi2D(result);
	}

	private static WB_Voronoi2D voronoi2D(final ArrayList<Coordinate> coords, final double d, final int c, final WB_Map2D context) {
		final DelaunayTriangulationBuilder dtb = new DelaunayTriangulationBuilder();
		dtb.setSites(coords);
		final QuadEdgeSubdivision qes = dtb.getSubdivision();
		final GeometryCollection polys = (GeometryCollection) qes.getVoronoiDiagram(WB_JTS.JTSgf);
		final int npolys = polys.getNumGeometries();
		final List<WB_VoronoiCell2D> result = new FastList<WB_VoronoiCell2D>();
		Coordinate[] coordsArray = new Coordinate[coords.size()];
		coordsArray = coords.toArray(coordsArray);
		for (int i = 0; i < npolys; ++i) {
			Geometry intersect;
			Polygon poly = (Polygon) (intersect = polys.getGeometryN(i));
			intersect = intersect.buffer(-d, c);
			final double cellindex = ((Coordinate) poly.getUserData()).z;
			for (int j = 0; j < intersect.getNumGeometries(); ++j) {
				if (intersect.getGeometryN(j).getGeometryType().equals("Polygon") && !intersect.getGeometryN(j).isEmpty()) {
					poly = (Polygon) intersect.getGeometryN(j);
					final Coordinate[] polycoord = poly.getCoordinates();
					final List<WB_Point> polypoints = new FastList<WB_Point>();
					Coordinate[] array;
					for (int length = (array = polycoord).length, k = 0; k < length; ++k) {
						final Coordinate element = array[k];
						polypoints.add(toPoint(element.x, element.y, context));
					}
					final Point centroid = poly.getCentroid();
					final WB_Point pc = (centroid == null) ? null : toPoint(centroid.getX(), centroid.getY(), context);
					final int index = (int) cellindex;
					final double area = poly.getArea();
					result.add(new WB_VoronoiCell2D(polypoints, index, createPoint2D(coords.get(index)), area, pc));
				}
			}
		}
		return new WB_Voronoi2D(result);
	}

	private static WB_Voronoi2D clippedVoronoi2D(final ArrayList<Coordinate> coords, final WB_Map2D context) {
		final DelaunayTriangulationBuilder dtb = new DelaunayTriangulationBuilder();
		dtb.setSites(coords);
		final QuadEdgeSubdivision qes = dtb.getSubdivision();
		final GeometryCollection polys = (GeometryCollection) qes.getVoronoiDiagram(WB_JTS.JTSgf);
		final int npolys = polys.getNumGeometries();
		final List<WB_VoronoiCell2D> result = new FastList<WB_VoronoiCell2D>();
		Coordinate[] coordsArray = new Coordinate[coords.size()];
		coordsArray = coords.toArray(coordsArray);
		final ConvexHull ch = new ConvexHull(coordsArray, WB_JTS.JTSgf);
		final Geometry hull = ch.getConvexHull();
		for (int i = 0; i < npolys; ++i) {
			Polygon poly = (Polygon) polys.getGeometryN(i);
			final Geometry intersect = poly.intersection(hull.getGeometryN(0));
			final double cellindex = ((Coordinate) poly.getUserData()).z;
			for (int j = 0; j < intersect.getNumGeometries(); ++j) {
				if (intersect.getGeometryN(j).getGeometryType().equals("Polygon") && !intersect.getGeometryN(j).isEmpty()) {
					poly = (Polygon) intersect.getGeometryN(j);
					final Coordinate[] polycoord = poly.getCoordinates();
					final List<WB_Point> polypoints = new FastList<WB_Point>();
					Coordinate[] array;
					for (int length = (array = polycoord).length, k = 0; k < length; ++k) {
						final Coordinate element = array[k];
						polypoints.add(toPoint(element.x, element.y, context));
					}
					final Point centroid = poly.getCentroid();
					final WB_Point pc = (centroid == null) ? null : toPoint(centroid.getX(), centroid.getY(), context);
					final int index = (int) cellindex;
					final double area = poly.getArea();
					result.add(new WB_VoronoiCell2D(polypoints, index, createPoint2D(coords.get(index)), area, pc));
				}
			}
		}
		return new WB_Voronoi2D(result);
	}

	private static WB_Voronoi2D clippedVoronoi2D(final ArrayList<Coordinate> coords, final double d, final int c, final WB_Map2D context) {
		final DelaunayTriangulationBuilder dtb = new DelaunayTriangulationBuilder();
		dtb.setSites(coords);
		final QuadEdgeSubdivision qes = dtb.getSubdivision();
		final GeometryCollection polys = (GeometryCollection) qes.getVoronoiDiagram(WB_JTS.JTSgf);
		final int npolys = polys.getNumGeometries();
		final List<WB_VoronoiCell2D> result = new FastList<WB_VoronoiCell2D>();
		Coordinate[] coordsArray = new Coordinate[coords.size()];
		coordsArray = coords.toArray(coordsArray);
		final ConvexHull ch = new ConvexHull(coordsArray, WB_JTS.JTSgf);
		final Geometry hull = ch.getConvexHull();
		for (int i = 0; i < npolys; ++i) {
			Polygon poly = (Polygon) polys.getGeometryN(i);
			Geometry intersect = poly.intersection(hull.getGeometryN(0));
			intersect = intersect.buffer(-d, c);
			final double cellindex = ((Coordinate) poly.getUserData()).z;
			for (int j = 0; j < intersect.getNumGeometries(); ++j) {
				if (intersect.getGeometryN(j).getGeometryType().equals("Polygon") && !intersect.getGeometryN(j).isEmpty()) {
					poly = (Polygon) intersect.getGeometryN(j);
					final Coordinate[] polycoord = poly.getCoordinates();
					final List<WB_Point> polypoints = new FastList<WB_Point>();
					Coordinate[] array;
					for (int length = (array = polycoord).length, k = 0; k < length; ++k) {
						final Coordinate element = array[k];
						polypoints.add(toPoint(element.x, element.y, context));
					}
					final Point centroid = poly.getCentroid();
					final WB_Point pc = (centroid == null) ? null : toPoint(centroid.getX(), centroid.getY(), context);
					final int index = (int) cellindex;
					final double area = poly.getArea();
					result.add(new WB_VoronoiCell2D(polypoints, index, createPoint2D(coords.get(index)), area, pc));
				}
			}
		}
		return new WB_Voronoi2D(result);
	}

	private static WB_Voronoi2D clippedVoronoi2D(final ArrayList<Coordinate> coords, final ArrayList<Coordinate> bdcoords,
			final WB_Map2D context) {
		final DelaunayTriangulationBuilder dtb = new DelaunayTriangulationBuilder();
		dtb.setSites(coords);
		final QuadEdgeSubdivision qes = dtb.getSubdivision();
		final GeometryCollection polys = (GeometryCollection) qes.getVoronoiDiagram(WB_JTS.JTSgf);
		final int npolys = polys.getNumGeometries();
		final List<WB_VoronoiCell2D> result = new FastList<WB_VoronoiCell2D>();
		Coordinate[] bdcoordsArray = new Coordinate[bdcoords.size()];
		bdcoordsArray = bdcoords.toArray(bdcoordsArray);
		final Polygon hull = WB_JTS.JTSgf.createPolygon(bdcoordsArray);
		for (int i = 0; i < npolys; ++i) {
			Polygon poly = (Polygon) polys.getGeometryN(i);
			final Geometry intersect = poly.intersection(hull);
			final double cellindex = ((Coordinate) poly.getUserData()).z;
			for (int j = 0; j < intersect.getNumGeometries(); ++j) {
				if (intersect.getGeometryN(j).getGeometryType().equals("Polygon") && !intersect.getGeometryN(j).isEmpty()) {
					poly = (Polygon) intersect.getGeometryN(j);
					final Coordinate[] polycoord = poly.getCoordinates();
					final List<WB_Point> polypoints = new FastList<WB_Point>();
					Coordinate[] array;
					for (int length = (array = polycoord).length, k = 0; k < length; ++k) {
						final Coordinate element = array[k];
						polypoints.add(toPoint(element.x, element.y, context));
					}
					final Point centroid = poly.getCentroid();
					final WB_Point pc = (centroid == null) ? null : toPoint(centroid.getX(), centroid.getY(), context);
					final int index = (int) cellindex;
					final double area = poly.getArea();
					result.add(new WB_VoronoiCell2D(polypoints, index, createPoint2D(coords.get(index)), area, pc));
				}
			}
		}
		return new WB_Voronoi2D(result);
	}

	private static WB_Voronoi2D clippedVoronoi2D(final ArrayList<Coordinate> coords, final WB_Polygon constraint, final WB_Map2D context) {
		final DelaunayTriangulationBuilder dtb = new DelaunayTriangulationBuilder();
		dtb.setSites(coords);
		final QuadEdgeSubdivision qes = dtb.getSubdivision();
		final GeometryCollection polys = (GeometryCollection) qes.getVoronoiDiagram(WB_JTS.JTSgf);
		final int npolys = polys.getNumGeometries();
		final List<WB_VoronoiCell2D> result = new FastList<WB_VoronoiCell2D>();
		final Polygon hull = toJTSPolygon2D(constraint);
		for (int i = 0; i < npolys; ++i) {
			Polygon poly = (Polygon) polys.getGeometryN(i);
			final Geometry intersect = poly.intersection(hull);
			final double cellindex = ((Coordinate) poly.getUserData()).z;
			for (int j = 0; j < intersect.getNumGeometries(); ++j) {
				if (intersect.getGeometryN(j).getGeometryType().equals("Polygon") && !intersect.getGeometryN(j).isEmpty()) {
					poly = (Polygon) intersect.getGeometryN(j);
					final Coordinate[] polycoord = poly.getCoordinates();
					final List<WB_Point> polypoints = new FastList<WB_Point>();
					Coordinate[] array;
					for (int length = (array = polycoord).length, k = 0; k < length; ++k) {
						final Coordinate element = array[k];
						polypoints.add(toPoint(element.x, element.y, context));
					}
					final Point centroid = poly.getCentroid();
					final WB_Point pc = (centroid == null) ? null : toPoint(centroid.getX(), centroid.getY(), context);
					final int index = (int) cellindex;
					final double area = poly.getArea();
					result.add(new WB_VoronoiCell2D(polypoints, index, createPoint2D(coords.get(index)), area, pc));
				}
			}
		}
		return new WB_Voronoi2D(result);
	}

	private static WB_Voronoi2D clippedVoronoi2D(final ArrayList<Coordinate> coords, final List<WB_Polygon> constraint,
			final WB_Map2D context) {
		final DelaunayTriangulationBuilder dtb = new DelaunayTriangulationBuilder();
		dtb.setSites(coords);
		final QuadEdgeSubdivision qes = dtb.getSubdivision();
		final GeometryCollection polys = (GeometryCollection) qes.getVoronoiDiagram(WB_JTS.JTSgf);
		final int npolys = polys.getNumGeometries();
		final List<WB_VoronoiCell2D> result = new FastList<WB_VoronoiCell2D>();
		final Geometry hull = toJTSMultiPolygon2D(constraint);
		for (int i = 0; i < npolys; ++i) {
			Polygon poly = (Polygon) polys.getGeometryN(i);
			final Geometry intersect = poly.intersection(hull);
			final double cellindex = ((Coordinate) poly.getUserData()).z;
			for (int j = 0; j < intersect.getNumGeometries(); ++j) {
				if (intersect.getGeometryN(j).getGeometryType().equals("Polygon") && !intersect.getGeometryN(j).isEmpty()) {
					poly = (Polygon) intersect.getGeometryN(j);
					final Coordinate[] polycoord = poly.getCoordinates();
					final List<WB_Point> polypoints = new FastList<WB_Point>();
					Coordinate[] array;
					for (int length = (array = polycoord).length, k = 0; k < length; ++k) {
						final Coordinate element = array[k];
						polypoints.add(toPoint(element.x, element.y, context));
					}
					final Point centroid = poly.getCentroid();
					final WB_Point pc = (centroid == null) ? null : toPoint(centroid.getX(), centroid.getY(), context);
					final int index = (int) cellindex;
					final double area = poly.getArea();
					result.add(new WB_VoronoiCell2D(polypoints, index, createPoint2D(coords.get(index)), area, pc));
				}
			}
		}
		return new WB_Voronoi2D(result);
	}

	private static WB_Voronoi2D clippedVoronoi2D(final ArrayList<Coordinate> coords, final ArrayList<Coordinate> bdcoords, final double d,
			final int c, final WB_Map2D context) {
		final DelaunayTriangulationBuilder dtb = new DelaunayTriangulationBuilder();
		dtb.setSites(coords);
		final QuadEdgeSubdivision qes = dtb.getSubdivision();
		final GeometryCollection polys = (GeometryCollection) qes.getVoronoiDiagram(WB_JTS.JTSgf);
		final int npolys = polys.getNumGeometries();
		final List<WB_VoronoiCell2D> result = new FastList<WB_VoronoiCell2D>();
		Coordinate[] bdcoordsArray = new Coordinate[bdcoords.size()];
		bdcoordsArray = bdcoords.toArray(bdcoordsArray);
		final Polygon hull = WB_JTS.JTSgf.createPolygon(bdcoordsArray);
		for (int i = 0; i < npolys; ++i) {
			Polygon poly = (Polygon) polys.getGeometryN(i);
			Geometry intersect = poly.intersection(hull);
			intersect = intersect.buffer(-d, c);
			final double cellindex = ((Coordinate) poly.getUserData()).z;
			for (int j = 0; j < intersect.getNumGeometries(); ++j) {
				if (intersect.getGeometryN(j).getGeometryType().equals("Polygon") && !intersect.getGeometryN(j).isEmpty()) {
					poly = (Polygon) intersect.getGeometryN(j);
					final Coordinate[] polycoord = poly.getCoordinates();
					final List<WB_Point> polypoints = new FastList<WB_Point>();
					Coordinate[] array;
					for (int length = (array = polycoord).length, k = 0; k < length; ++k) {
						final Coordinate element = array[k];
						polypoints.add(toPoint(element.x, element.y, context));
					}
					final Point centroid = poly.getCentroid();
					final WB_Point pc = (centroid == null) ? null : toPoint(centroid.getX(), centroid.getY(), context);
					final int index = (int) cellindex;
					final double area = poly.getArea();
					result.add(new WB_VoronoiCell2D(polypoints, index, createPoint2D(coords.get(index)), area, pc));
				}
			}
		}
		return new WB_Voronoi2D(result);
	}

	private static WB_Voronoi2D clippedVoronoi2D(final ArrayList<Coordinate> coords, final WB_Polygon constraint, final double d,
			final int c, final WB_Map2D context) {
		final DelaunayTriangulationBuilder dtb = new DelaunayTriangulationBuilder();
		dtb.setSites(coords);
		final QuadEdgeSubdivision qes = dtb.getSubdivision();
		final GeometryCollection polys = (GeometryCollection) qes.getVoronoiDiagram(WB_JTS.JTSgf);
		final int npolys = polys.getNumGeometries();
		final List<WB_VoronoiCell2D> result = new FastList<WB_VoronoiCell2D>();
		final Polygon hull = toJTSPolygon2D(constraint);
		for (int i = 0; i < npolys; ++i) {
			Polygon poly = (Polygon) polys.getGeometryN(i);
			Geometry intersect = poly.intersection(hull);
			intersect = intersect.buffer(-d, c);
			final double cellindex = ((Coordinate) poly.getUserData()).z;
			for (int j = 0; j < intersect.getNumGeometries(); ++j) {
				if (intersect.getGeometryN(j).getGeometryType().equals("Polygon") && !intersect.getGeometryN(j).isEmpty()) {
					poly = (Polygon) intersect.getGeometryN(j);
					final Coordinate[] polycoord = poly.getCoordinates();
					final List<WB_Point> polypoints = new FastList<WB_Point>();
					Coordinate[] array;
					for (int length = (array = polycoord).length, k = 0; k < length; ++k) {
						final Coordinate element = array[k];
						polypoints.add(toPoint(element.x, element.y, context));
					}
					final Point centroid = poly.getCentroid();
					final WB_Point pc = (centroid == null) ? null : toPoint(centroid.getX(), centroid.getY(), context);
					final int index = (int) cellindex;
					final double area = poly.getArea();
					result.add(new WB_VoronoiCell2D(polypoints, index, createPoint2D(coords.get(index)), area, pc));
				}
			}
		}
		return new WB_Voronoi2D(result);
	}

	private static WB_Voronoi2D clippedVoronoi2D(final ArrayList<Coordinate> coords, final List<WB_Polygon> constraint, final double d,
			final int c, final WB_Map2D context) {
		final DelaunayTriangulationBuilder dtb = new DelaunayTriangulationBuilder();
		dtb.setSites(coords);
		final QuadEdgeSubdivision qes = dtb.getSubdivision();
		final GeometryCollection polys = (GeometryCollection) qes.getVoronoiDiagram(WB_JTS.JTSgf);
		final int npolys = polys.getNumGeometries();
		final List<WB_VoronoiCell2D> result = new FastList<WB_VoronoiCell2D>();
		final Geometry hull = toJTSMultiPolygon2D(constraint);
		for (int i = 0; i < npolys; ++i) {
			Polygon poly = (Polygon) polys.getGeometryN(i);
			Geometry intersect = poly.intersection(hull);
			intersect = intersect.buffer(-d, c);
			final double cellindex = ((Coordinate) poly.getUserData()).z;
			for (int j = 0; j < intersect.getNumGeometries(); ++j) {
				if (intersect.getGeometryN(j).getGeometryType().equals("Polygon") && !intersect.getGeometryN(j).isEmpty()) {
					poly = (Polygon) intersect.getGeometryN(j);
					final Coordinate[] polycoord = poly.getCoordinates();
					final List<WB_Point> polypoints = new FastList<WB_Point>();
					Coordinate[] array;
					for (int length = (array = polycoord).length, k = 0; k < length; ++k) {
						final Coordinate element = array[k];
						polypoints.add(toPoint(element.x, element.y, context));
					}
					final Point centroid = poly.getCentroid();
					final WB_Point pc = (centroid == null) ? null : toPoint(centroid.getX(), centroid.getY(), context);
					final int index = (int) cellindex;
					final double area = poly.getArea();
					result.add(new WB_VoronoiCell2D(polypoints, index, createPoint2D(coords.get(index)), area, pc));
				}
			}
		}
		return new WB_Voronoi2D(result);
	}

	private static Coordinate toCoordinate(final WB_Coord p, final int i, final WB_Map2D context) {
		final WB_Point tmp = new WB_Point();
		context.mapPoint3D(p, tmp);
		final Coordinate c = new Coordinate(tmp.xd(), tmp.yd(), i);
		return c;
	}

	private static WB_Point toPoint(final double x, final double y, final WB_Map2D context) {
		final WB_Point tmp = new WB_Point();
		context.unmapPoint3D(x, y, 0.0, tmp);
		return tmp;
	}

	static WB_Voronoi2D getVoronoi2D(final WB_Coord[] points, final WB_Map2D context) {
		final int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		for (int i = 0; i < n; ++i) {
			coords.add(toCoordinate(points[i], i, context));
		}
		return voronoi2D(coords, context);
	}

	static WB_Voronoi2D getVoronoi2D(final Collection<? extends WB_Coord> points, final WB_Map2D context) {
		final int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, context));
			++id;
		}
		return voronoi2D(coords, context);
	}

	static WB_Voronoi2D getVoronoi2D(final WB_Coord[] points, final double d, final int c, final WB_Map2D context) {
		final int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		for (int i = 0; i < n; ++i) {
			coords.add(toCoordinate(points[i], i, context));
		}
		return voronoi2D(coords, d, c, context);
	}

	static WB_Voronoi2D getVoronoi2D(final Collection<? extends WB_Coord> points, final double d, final int c, final WB_Map2D context) {
		final int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, context));
			++id;
		}
		return voronoi2D(coords, d, c, context);
	}

	static WB_Voronoi2D getVoronoi2D(final Collection<? extends WB_Coord> points, final double d, final int c) {
		final int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, WB_JTS.XY));
			++id;
		}
		return voronoi2D(coords, d, c, WB_JTS.XY);
	}

	static WB_Voronoi2D getVoronoi2D(final Collection<? extends WB_Coord> points) {
		final int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, WB_JTS.XY));
			++id;
		}
		return voronoi2D(coords, WB_JTS.XY);
	}

	static WB_Voronoi2D getVoronoi2D(final WB_Coord[] points, final double d, final int c) {
		final int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		for (int i = 0; i < n; ++i) {
			coords.add(toCoordinate(points[i], i, WB_JTS.XY));
		}
		return voronoi2D(coords, d, c, WB_JTS.XY);
	}

	static WB_Voronoi2D getVoronoi2D(final WB_Coord[] points) {
		final int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		for (int i = 0; i < n; ++i) {
			coords.add(toCoordinate(points[i], i, WB_JTS.XY));
		}
		return voronoi2D(coords, WB_JTS.XY);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final WB_Coord[] points, final WB_Map2D context) {
		final int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		for (int i = 0; i < n; ++i) {
			coords.add(toCoordinate(points[i], i, context));
		}
		return clippedVoronoi2D(coords, context);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final WB_Coord[] points, final WB_Coord[] boundary, final WB_Map2D context) {
		int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		for (int i = 0; i < n; ++i) {
			coords.add(toCoordinate(points[i], i, context));
		}
		n = boundary.length;
		final ArrayList<Coordinate> bdcoords = new ArrayList<>(n);
		for (int j = 0; j < n; ++j) {
			bdcoords.add(toCoordinate(boundary[j], j, context));
		}
		if (!bdcoords.get(0).equals(bdcoords.get(n - 1))) {
			bdcoords.add(bdcoords.get(0));
		}
		return clippedVoronoi2D(coords, bdcoords, context);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final WB_Coord[] points, final WB_Polygon boundary, final WB_Map2D context) {
		final int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		for (int i = 0; i < n; ++i) {
			coords.add(toCoordinate(points[i], i, context));
		}
		return clippedVoronoi2D(coords, boundary, context);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final Collection<? extends WB_Coord> points, final WB_Map2D context) {
		final int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, context));
			++id;
		}
		return clippedVoronoi2D(coords, context);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final Collection<? extends WB_Coord> points, final Collection<? extends WB_Coord> boundary,
			final WB_Map2D context) {
		int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, context));
			++id;
		}
		n = boundary.size();
		final ArrayList<Coordinate> bdcoords = new ArrayList<>(n);
		id = 0;
		for (final WB_Coord p2 : boundary) {
			bdcoords.add(toCoordinate(p2, id, context));
			++id;
		}
		if (!bdcoords.get(0).equals(bdcoords.get(n - 1))) {
			bdcoords.add(bdcoords.get(0));
		}
		return clippedVoronoi2D(coords, bdcoords, context);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final Collection<? extends WB_Coord> points, final WB_Polygon boundary,
			final WB_Map2D context) {
		final int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, context));
			++id;
		}
		return clippedVoronoi2D(coords, boundary, context);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final Collection<? extends WB_Coord> points, final List<WB_Polygon> boundary,
			final WB_Map2D context) {
		final int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, context));
			++id;
		}
		return clippedVoronoi2D(coords, boundary, context);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final WB_Coord[] points, final double d, final WB_Map2D context) {
		final int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		for (int i = 0; i < n; ++i) {
			coords.add(toCoordinate(points[i], i, context));
		}
		return clippedVoronoi2D(coords, d, 2, context);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final WB_Coord[] points, final WB_Coord[] boundary, final double d, final WB_Map2D context) {
		int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		for (int i = 0; i < n; ++i) {
			coords.add(toCoordinate(points[i], i, context));
		}
		n = boundary.length;
		final ArrayList<Coordinate> bdcoords = new ArrayList<>(n);
		for (int j = 0; j < n; ++j) {
			bdcoords.add(toCoordinate(boundary[j], j, context));
		}
		if (!bdcoords.get(0).equals(bdcoords.get(n - 1))) {
			bdcoords.add(bdcoords.get(0));
		}
		return clippedVoronoi2D(coords, bdcoords, d, 2, context);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final WB_Coord[] points, final WB_Polygon boundary, final double d, final WB_Map2D context) {
		final int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		for (int i = 0; i < n; ++i) {
			coords.add(toCoordinate(points[i], i, context));
		}
		return clippedVoronoi2D(coords, boundary, d, 2, context);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final Collection<? extends WB_Coord> points, final Collection<? extends WB_Coord> boundary,
			final double d, final WB_Map2D context) {
		int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, context));
			++id;
		}
		n = boundary.size();
		final ArrayList<Coordinate> bdcoords = new ArrayList<>(n);
		id = 0;
		for (final WB_Coord p2 : boundary) {
			bdcoords.add(toCoordinate(p2, id, context));
			++id;
		}
		if (!bdcoords.get(0).equals(bdcoords.get(n - 1))) {
			bdcoords.add(bdcoords.get(0));
		}
		return clippedVoronoi2D(coords, bdcoords, d, 2, context);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final Collection<? extends WB_Coord> points, final WB_Polygon boundary, final double d,
			final WB_Map2D context) {
		final int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, context));
			++id;
		}
		return clippedVoronoi2D(coords, boundary, d, 2, context);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final Collection<? extends WB_Coord> points, final List<WB_Polygon> boundary, final double d,
			final WB_Map2D context) {
		final int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, context));
			++id;
		}
		return clippedVoronoi2D(coords, boundary, d, 2, context);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final Collection<? extends WB_Coord> points, final double d, final WB_Map2D context) {
		final int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, context));
			++id;
		}
		return clippedVoronoi2D(coords, d, 2, context);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final WB_Coord[] points, final double d, final int c, final WB_Map2D context) {
		final int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		for (int i = 0; i < n; ++i) {
			coords.add(toCoordinate(points[i], i, context));
		}
		return clippedVoronoi2D(coords, d, c, context);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final WB_Coord[] points, final WB_Coord[] boundary, final double d, final int c,
			final WB_Map2D context) {
		int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		for (int i = 0; i < n; ++i) {
			coords.add(toCoordinate(points[i], i, context));
		}
		n = boundary.length;
		final ArrayList<Coordinate> bdcoords = new ArrayList<>(n);
		for (int j = 0; j < n; ++j) {
			bdcoords.add(toCoordinate(boundary[j], j, context));
		}
		if (!bdcoords.get(0).equals(bdcoords.get(n - 1))) {
			bdcoords.add(bdcoords.get(0));
		}
		return clippedVoronoi2D(coords, bdcoords, d, c, context);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final WB_Coord[] points, final WB_Polygon boundary, final double d, final int c,
			final WB_Map2D context) {
		final int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		for (int i = 0; i < n; ++i) {
			coords.add(toCoordinate(points[i], i, context));
		}
		return clippedVoronoi2D(coords, boundary, d, c, context);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final Collection<? extends WB_Coord> points, final double d, final int c,
			final WB_Map2D context) {
		final int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, context));
			++id;
		}
		return clippedVoronoi2D(coords, d, c, context);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final Collection<? extends WB_Coord> points, final Collection<? extends WB_Coord> boundary,
			final double d, final int c, final WB_Map2D context) {
		int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, context));
			++id;
		}
		n = boundary.size();
		final ArrayList<Coordinate> bdcoords = new ArrayList<>(n);
		id = 0;
		for (final WB_Coord p2 : boundary) {
			bdcoords.add(toCoordinate(p2, id, context));
			++id;
		}
		if (!bdcoords.get(0).equals(bdcoords.get(n - 1))) {
			bdcoords.add(bdcoords.get(0));
		}
		return clippedVoronoi2D(coords, bdcoords, d, c, context);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final Collection<? extends WB_Coord> points, final WB_Polygon boundary, final double d,
			final int c, final WB_Map2D context) {
		final int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, context));
			++id;
		}
		return clippedVoronoi2D(coords, boundary, d, c, context);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final Collection<? extends WB_Coord> points, final List<WB_Polygon> boundary, final double d,
			final int c, final WB_Map2D context) {
		final int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, context));
			++id;
		}
		return clippedVoronoi2D(coords, boundary, d, c, context);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final Collection<? extends WB_Coord> points, final Collection<? extends WB_Coord> boundary,
			final double d, final int c) {
		int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, WB_JTS.XY));
			++id;
		}
		n = boundary.size();
		final ArrayList<Coordinate> bdcoords = new ArrayList<>(n);
		id = 0;
		for (final WB_Coord p2 : boundary) {
			bdcoords.add(toCoordinate(p2, id, WB_JTS.XY));
			++id;
		}
		if (!bdcoords.get(0).equals(bdcoords.get(n - 1))) {
			bdcoords.add(bdcoords.get(0));
		}
		return clippedVoronoi2D(coords, bdcoords, d, c, WB_JTS.XY);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final Collection<? extends WB_Coord> points, final Collection<? extends WB_Coord> boundary,
			final double d) {
		int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, WB_JTS.XY));
			++id;
		}
		n = boundary.size();
		final ArrayList<Coordinate> bdcoords = new ArrayList<>(n);
		id = 0;
		for (final WB_Coord p2 : boundary) {
			bdcoords.add(toCoordinate(p2, id, WB_JTS.XY));
			++id;
		}
		if (!bdcoords.get(0).equals(bdcoords.get(n - 1))) {
			bdcoords.add(bdcoords.get(0));
		}
		return clippedVoronoi2D(coords, bdcoords, d, 2, WB_JTS.XY);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final Collection<? extends WB_Coord> points, final Collection<? extends WB_Coord> boundary) {
		int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, WB_JTS.XY));
			++id;
		}
		n = boundary.size();
		final ArrayList<Coordinate> bdcoords = new ArrayList<>(n);
		id = 0;
		for (final WB_Coord p2 : boundary) {
			bdcoords.add(toCoordinate(p2, id, WB_JTS.XY));
			++id;
		}
		if (!bdcoords.get(0).equals(bdcoords.get(n - 1))) {
			bdcoords.add(bdcoords.get(0));
		}
		return clippedVoronoi2D(coords, bdcoords, WB_JTS.XY);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final Collection<? extends WB_Coord> points, final double d, final int c) {
		final int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, WB_JTS.XY));
			++id;
		}
		return clippedVoronoi2D(coords, d, c, WB_JTS.XY);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final Collection<? extends WB_Coord> points, final double d) {
		final int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, WB_JTS.XY));
			++id;
		}
		return clippedVoronoi2D(coords, d, 2, WB_JTS.XY);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final Collection<? extends WB_Coord> points) {
		final int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, WB_JTS.XY));
			++id;
		}
		return clippedVoronoi2D(coords, WB_JTS.XY);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final WB_Coord[] points, final double d, final int c) {
		final int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		for (int i = 0; i < n; ++i) {
			coords.add(toCoordinate(points[i], i, WB_JTS.XY));
		}
		return clippedVoronoi2D(coords, d, c, WB_JTS.XY);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final WB_Coord[] points, final double d) {
		final int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		for (int i = 0; i < n; ++i) {
			coords.add(toCoordinate(points[i], i, WB_JTS.XY));
		}
		return clippedVoronoi2D(coords, d, 2, WB_JTS.XY);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final WB_Coord[] points, final WB_Coord[] boundary, final double d, final int c) {
		int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		for (int i = 0; i < n; ++i) {
			coords.add(toCoordinate(points[i], i, WB_JTS.XY));
		}
		n = boundary.length;
		final ArrayList<Coordinate> bdcoords = new ArrayList<>(n);
		for (int j = 0; j < n; ++j) {
			bdcoords.add(toCoordinate(boundary[j], j, WB_JTS.XY));
		}
		if (!bdcoords.get(0).equals(bdcoords.get(n - 1))) {
			bdcoords.add(bdcoords.get(0));
		}
		return clippedVoronoi2D(coords, bdcoords, d, c, WB_JTS.XY);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final WB_Coord[] points, final WB_Coord[] boundary, final double d) {
		int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		for (int i = 0; i < n; ++i) {
			coords.add(toCoordinate(points[i], i, WB_JTS.XY));
		}
		n = boundary.length;
		final ArrayList<Coordinate> bdcoords = new ArrayList<>(n);
		for (int j = 0; j < n; ++j) {
			bdcoords.add(toCoordinate(boundary[j], j, WB_JTS.XY));
		}
		if (!bdcoords.get(0).equals(bdcoords.get(n - 1))) {
			bdcoords.add(bdcoords.get(0));
		}
		return clippedVoronoi2D(coords, bdcoords, d, 2, WB_JTS.XY);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final WB_Coord[] points, final WB_Coord[] boundary) {
		int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		for (int i = 0; i < n; ++i) {
			coords.add(toCoordinate(points[i], i, WB_JTS.XY));
		}
		n = boundary.length;
		final ArrayList<Coordinate> bdcoords = new ArrayList<>(n);
		for (int j = 0; j < n; ++j) {
			bdcoords.add(toCoordinate(boundary[j], j, WB_JTS.XY));
		}
		if (!bdcoords.get(0).equals(bdcoords.get(n - 1))) {
			bdcoords.add(bdcoords.get(0));
		}
		return clippedVoronoi2D(coords, bdcoords, WB_JTS.XY);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final WB_Coord[] points) {
		final int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		for (int i = 0; i < n; ++i) {
			coords.add(toCoordinate(points[i], i, WB_JTS.XY));
		}
		return clippedVoronoi2D(coords, WB_JTS.XY);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final WB_Coord[] points, final WB_Polygon boundary) {
		final int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		for (int i = 0; i < n; ++i) {
			coords.add(toCoordinate(points[i], i, WB_JTS.XY));
		}
		return clippedVoronoi2D(coords, boundary, WB_JTS.XY);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final WB_Coord[] points, final List<WB_Polygon> boundary) {
		final int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		for (int i = 0; i < n; ++i) {
			coords.add(toCoordinate(points[i], i, WB_JTS.XY));
		}
		return clippedVoronoi2D(coords, boundary, WB_JTS.XY);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final Collection<? extends WB_Coord> points, final WB_Polygon boundary) {
		final int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, WB_JTS.XY));
			++id;
		}
		return clippedVoronoi2D(coords, boundary, WB_JTS.XY);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final Collection<? extends WB_Coord> points, final List<WB_Polygon> boundary) {
		final int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, WB_JTS.XY));
			++id;
		}
		return clippedVoronoi2D(coords, boundary, WB_JTS.XY);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final WB_Coord[] points, final WB_Polygon boundary, final double d) {
		final int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		for (int i = 0; i < n; ++i) {
			coords.add(toCoordinate(points[i], i, WB_JTS.XY));
		}
		return clippedVoronoi2D(coords, boundary, d, 2, WB_JTS.XY);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final WB_Coord[] points, final List<WB_Polygon> boundary, final double d) {
		final int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		for (int i = 0; i < n; ++i) {
			coords.add(toCoordinate(points[i], i, WB_JTS.XY));
		}
		return clippedVoronoi2D(coords, boundary, d, 2, WB_JTS.XY);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final Collection<? extends WB_Coord> points, final WB_Polygon boundary, final double d) {
		final int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, WB_JTS.XY));
			++id;
		}
		return clippedVoronoi2D(coords, boundary, d, 2, WB_JTS.XY);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final Collection<? extends WB_Coord> points, final List<WB_Polygon> boundary, final double d) {
		final int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, WB_JTS.XY));
			++id;
		}
		return clippedVoronoi2D(coords, boundary, d, 2, WB_JTS.XY);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final WB_Coord[] points, final WB_Polygon boundary, final double d, final int c) {
		final int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		for (int i = 0; i < n; ++i) {
			coords.add(toCoordinate(points[i], i, WB_JTS.XY));
		}
		return clippedVoronoi2D(coords, boundary, d, c, WB_JTS.XY);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final WB_Coord[] points, final List<WB_Polygon> boundary, final double d, final int c) {
		final int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		for (int i = 0; i < n; ++i) {
			coords.add(toCoordinate(points[i], i, WB_JTS.XY));
		}
		return clippedVoronoi2D(coords, boundary, d, c, WB_JTS.XY);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final Collection<? extends WB_Coord> points, final WB_Polygon boundary, final double d,
			final int c) {
		final int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, WB_JTS.XY));
			++id;
		}
		return clippedVoronoi2D(coords, boundary, d, c, WB_JTS.XY);
	}

	static WB_Voronoi2D getClippedVoronoi2D(final Collection<? extends WB_Coord> points, final List<WB_Polygon> boundary, final double d,
			final int c) {
		final int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, WB_JTS.XY));
			++id;
		}
		return clippedVoronoi2D(coords, boundary, d, c, WB_JTS.XY);
	}

	public static class PlanarPathTriangulator {
		static final WB_ProgressReporter.WB_ProgressTracker tracker;
		private static WB_GeometryFactory2D geometryfactory;

		static {
			tracker = WB_ProgressReporter.WB_ProgressTracker.instance();
			PlanarPathTriangulator.geometryfactory = new WB_GeometryFactory2D();
		}

		public static long[] getTriangleKeys(final List<? extends HE_Path> paths, final WB_Plane P) {
			PlanarPathTriangulator.tracker.setStartStatus("HET_PlanarPathTriangulator", "Starting planar path triangulation.");
			final WB_Map2D emb = new WB_OrthoProject(P);
			final RingTree ringtree = new RingTree();
			final WB_KDTree3D<WB_Point, Long> vertextree = new WB_KDTree3D<WB_Point, Long>();
			PlanarPathTriangulator.tracker.setDuringStatus("HET_PlanarPathTriangulator", "Building contours tree.");
			for (HE_Path path2 : paths) {
				final HE_Path path = path2;
				if (path.isLoop() && path.getPathOrder() > 2) {
					final List<HE_Vertex> vertices = path.getPathVertices();
					final Coordinate[] pts = new Coordinate[vertices.size() + 1];
					for (int j = 0; j < vertices.size(); ++j) {
						final WB_Point proj = new WB_Point();
						emb.mapPoint3D(vertices.get(j), proj);
						vertextree.add(proj, vertices.get(j).getKey());
						pts[vertices.size() - j] = new Coordinate(proj.xd(), proj.yd(), 0.0);
					}
					final WB_Point proj2 = new WB_Point();
					emb.mapPoint3D(vertices.get(0), proj2);
					pts[0] = new Coordinate(proj2.xd(), proj2.yd(), 0.0);
					ringtree.add(WB_JTS.JTSgf.createLinearRing(pts));
				}
			}
			PlanarPathTriangulator.tracker.setDuringStatus("HET_PlanarPathTriangulator", "Extracting polygons from contours tree.");
			final List<WB_Polygon> polygons = ringtree.extractPolygons();
			final List<WB_Coord[]> triangles = new FastList<WB_Coord[]>();
			final WB_ProgressReporter.WB_ProgressCounter counter = new WB_ProgressReporter.WB_ProgressCounter(polygons.size(), 10);
			PlanarPathTriangulator.tracker.setCounterStatus("HET_PlanarPathTriangulator", "Triangulating polygons.", counter);
			for (final WB_Polygon poly : polygons) {
				final int[] tris = poly.getTriangles();
				for (int k = 0; k < tris.length; k += 3) {
					triangles.add(new WB_Coord[] { poly.getPoint(tris[k]), poly.getPoint(tris[k + 1]), poly.getPoint(tris[k + 2]) });
				}
				counter.increment();
			}
			final long[] trianglekeys = new long[3 * triangles.size()];
			for (int l = 0; l < triangles.size(); ++l) {
				final WB_Coord[] tri = triangles.get(l);
				final long key0 = vertextree.getNearestNeighbor(tri[0]).value;
				final long key2 = vertextree.getNearestNeighbor(tri[1]).value;
				final long key3 = vertextree.getNearestNeighbor(tri[2]).value;
				trianglekeys[3 * l] = key0;
				trianglekeys[3 * l + 1] = key3;
				trianglekeys[3 * l + 2] = key2;
			}
			PlanarPathTriangulator.tracker.setStopStatus("HET_PlanarPathTriangulator", "All paths triangulated.");
			return trianglekeys;
		}

		public static HE_Mesh toMesh(final WB_Triangulation2DWithPoints tri) {
			final HEC_FromFacelist ffl = new HEC_FromFacelist().setFaces(tri.getTriangles()).setVertices(tri.getPoints());
			return new HE_Mesh(ffl);
		}

		private static class RingNode {
			RingNode parent;
			List<RingNode> children;
			LinearRing ring;
			Polygon poly;
			boolean hole;

			RingNode() {
				this.parent = null;
				this.ring = null;
				this.children = new ArrayList<>();
				this.hole = true;
			}

			RingNode(final RingNode parent, final LinearRing ring) {
				this.parent = parent;
				this.ring = ring;
				final Coordinate[] coords = ring.getCoordinates();
				this.poly = WB_JTS.JTSgf.createPolygon(coords);
				this.hole = !CGAlgorithms.isCCW(coords);
				this.children = new ArrayList<>();
			}
		}

		private static class RingTree {
			RingNode root;

			RingTree() {
				this.root = new RingNode();
			}

			void add(final LinearRing ring) {
				final Polygon poly = WB_JTS.JTSgf.createPolygon(ring);
				RingNode currentParent = this.root;
				RingNode foundParent;
				do {
					foundParent = null;
					for (final RingNode node : currentParent.children) {
						final Polygon other = node.poly;
						if (poly.within(other)) {
							foundParent = node;
							currentParent = node;
							break;
						}
					}
				} while (foundParent != null);
				final RingNode newNode = new RingNode(currentParent, ring);
				final List<RingNode> nodesToRemove = new ArrayList<>();
				for (int j = 0; j < currentParent.children.size(); ++j) {
					final RingNode node2 = currentParent.children.get(j);
					final Polygon other2 = node2.poly;
					if (other2.within(poly)) {
						newNode.children.add(node2);
						nodesToRemove.add(node2);
					}
				}
				currentParent.children.removeAll(nodesToRemove);
				currentParent.children.add(newNode);
			}

			List<WB_Polygon> extractPolygons() {
				final List<WB_Polygon> polygons = new ArrayList<>();
				final List<RingNode> shellNodes = new ArrayList<>();
				this.addExteriorNodes(this.root, shellNodes);
				for (final RingNode node : shellNodes) {
					int count = 0;
					for (RingNode element : node.children) {
						if (element.hole) {
							++count;
						}
					}
					final LinearRing[] holes = new LinearRing[count];
					int index = 0;
					for (RingNode element : node.children) {
						if (element.hole) {
							holes[index++] = element.ring;
						}
					}
					final Geometry result = WB_JTS.JTSgf.createPolygon(node.ring, holes);
					if (result.getGeometryType().equals("Polygon")) {
						polygons.add(WB_JTS.createPolygonFromJTSPolygon2D((Polygon) result));
					} else {
						if (!result.getGeometryType().equals("MultiPolygon")) {
							continue;
						}
						for (int k = 0; k < result.getNumGeometries(); ++k) {
							final Geometry ggeo = result.getGeometryN(k);
							polygons.add(WB_JTS.createPolygonFromJTSPolygon2D((Polygon) ggeo));
						}
					}
				}
				return polygons;
			}

			void addExteriorNodes(final RingNode parent, final List<RingNode> shellNodes) {
				for (final RingNode node : parent.children) {
					if (!node.hole) {
						shellNodes.add(node);
					}
					this.addExteriorNodes(node, shellNodes);
				}
			}
		}
	}

	public static class PolygonTriangulatorJTS {
		private List<Coordinate> shellCoords;
		private boolean[] shellCoordAvailable;

		public static int[] triangulateQuad(final WB_Coord p0, final WB_Coord p1, final WB_Coord p2, final WB_Coord p3) {
			final boolean p0inside = WB_GeometryOp.pointInTriangleBary3D(p0, p1, p2, p3);
			if (p0inside) {
				return new int[] { 0, 1, 2, 0, 2, 3 };
			}
			final boolean p2inside = WB_GeometryOp.pointInTriangleBary3D(p2, p1, p0, p3);
			if (p2inside) {
				return new int[] { 0, 1, 2, 0, 2, 3 };
			}
			return new int[] { 0, 1, 3, 1, 2, 3 };
		}

		public WB_Triangulation2D triangulatePolygon2D(final WB_Polygon polygon, final boolean optimize) {
			final List<WB_Point> pts = new FastList<WB_Point>();
			for (int i = 0; i < polygon.numberOfShellPoints; ++i) {
				pts.add(polygon.getPoint(i));
			}
			int index = polygon.numberOfShellPoints;
			final List[] hpts = new FastList[polygon.numberOfContours - 1];
			for (int j = 0; j < polygon.numberOfContours - 1; ++j) {
				hpts[j] = new FastList<Object>();
				for (int k = 0; k < polygon.numberOfPointsPerContour[j + 1]; ++k) {
					hpts[j].add(polygon.points.get(index++));
				}
			}
			final WB_Plane P = polygon.getPlane(0.0);
			final WB_Triangulation2DWithPoints triangulation = this.triangulatePolygon2D(pts, hpts, optimize,
					new WB_GeometryFactory().createEmbeddedPlane(P));
			final WB_KDTreeInteger3D<WB_Point> pointmap = new WB_KDTreeInteger3D<WB_Point>();
			for (int l = 0; l < polygon.numberOfPoints; ++l) {
				pointmap.add(polygon.getPoint(l), l);
			}
			final int[] triangles = triangulation.getTriangles();
			final int[] edges = triangulation.getEdges();
			final WB_CoordCollection tripoints = triangulation.getPoints();
			final int[] intmap = new int[tripoints.size()];
			index = 0;
			for (int m = 0; m < tripoints.size(); ++m) {
				final int found = pointmap.getNearestNeighbor(tripoints.get(m)).value;
				intmap[index++] = found;
			}
			for (int m = 0; m < triangles.length; ++m) {
				triangles[m] = intmap[triangles[m]];
			}
			for (int m = 0; m < edges.length; ++m) {
				edges[m] = intmap[edges[m]];
			}
			return new WB_Triangulation2D(triangles, edges);
		}

		public WB_Triangulation2DWithPoints triangulatePolygon2D(final List<? extends WB_Coord> outerPolygon,
				final List<? extends WB_Coord>[] innerPolygons, final boolean optimize, final WB_Map2D context) {
			final Coordinate[] coords = new Coordinate[outerPolygon.size() + 1];
			WB_Point point = new WB_Point();
			for (int i = 0; i < outerPolygon.size(); ++i) {
				context.mapPoint3D(outerPolygon.get(i), point);
				coords[i] = new Coordinate(point.xd(), point.yd(), 0.0);
			}
			context.mapPoint3D(outerPolygon.get(0), point);
			coords[outerPolygon.size()] = new Coordinate(point.xd(), point.yd(), 0.0);
			LinearRing[] holes = null;
			if (innerPolygons != null) {
				holes = new LinearRing[innerPolygons.length];
				for (int j = 0; j < innerPolygons.length; ++j) {
					final Coordinate[] icoords = new Coordinate[innerPolygons[j].size() + 1];
					for (int k = 0; k < innerPolygons[j].size(); ++k) {
						context.mapPoint3D(innerPolygons[j].get(k), point);
						icoords[k] = new Coordinate(point.xd(), point.yd(), 0.0);
					}
					context.mapPoint3D(innerPolygons[j].get(0), point);
					icoords[innerPolygons[j].size()] = new Coordinate(point.xd(), point.yd(), 0.0);
					final LinearRing hole = WB_JTS.JTSgf.createLinearRing(icoords);
					holes[j] = hole;
				}
			}
			final LinearRing shell = WB_JTS.JTSgf.createLinearRing(coords);
			final Polygon inputPolygon = WB_JTS.JTSgf.createPolygon(shell, holes);
			final int[] ears = this.triangulate(inputPolygon, optimize);
			final int[] E = extractEdgesTri(ears);
			final List<WB_Point> Points = new FastList<WB_Point>();
			for (int l = 0; l < this.shellCoords.size() - 1; ++l) {
				point = new WB_Point();
				context.unmapPoint2D(this.shellCoords.get(l).x, this.shellCoords.get(l).y, point);
				Points.add(point);
			}
			return new WB_Triangulation2DWithPoints(ears, E, Points);
		}

		public WB_Triangulation2DWithPoints triangulatePolygon2D(final int[] polygon, final WB_Coord[] points, final boolean optimize,
				final WB_Map2D context) {
			final Coordinate[] coords = new Coordinate[polygon.length + 1];
			WB_Point point = new WB_Point();
			for (int i = 0; i < polygon.length; ++i) {
				context.mapPoint3D(points[polygon[i]], point);
				coords[i] = new Coordinate(point.xd(), point.yd(), i);
			}
			context.mapPoint3D(points[polygon[0]], point);
			coords[polygon.length] = new Coordinate(point.xd(), point.yd(), 0.0);
			final Polygon inputPolygon = WB_JTS.JTSgf.createPolygon(coords);
			final int[] ears = this.triangulate(inputPolygon, optimize);
			for (int j = 0; j < ears.length; ++j) {
				ears[j] = polygon[ears[j]];
			}
			final int[] E = extractEdgesTri(ears);
			final List<WB_Point> Points = new FastList<WB_Point>();
			for (int k = 0; k < this.shellCoords.size() - 1; ++k) {
				point = new WB_Point();
				context.unmapPoint2D(this.shellCoords.get(k).x, this.shellCoords.get(k).y, point);
				Points.add(point);
			}
			return new WB_Triangulation2DWithPoints(ears, E, Points);
		}

		public WB_Triangulation2D triangulatePolygon2D(final int[] polygon, final List<? extends WB_Coord> points, final boolean optimize,
				final WB_Map2D context) {
			final Coordinate[] coords = new Coordinate[polygon.length + 1];
			final WB_Point point = new WB_Point();
			for (int i = 0; i < polygon.length; ++i) {
				context.mapPoint3D(points.get(polygon[i]), point);
				coords[i] = new Coordinate(point.xd(), point.yd(), polygon[i]);
			}
			context.mapPoint3D(points.get(polygon[0]), point);
			coords[polygon.length] = new Coordinate(point.xd(), point.yd(), polygon[0]);
			final Polygon inputPolygon = WB_JTS.JTSgf.createPolygon(coords);
			final int[] ears = this.triangulate(inputPolygon, optimize);
			for (int j = 0; j < ears.length; ++j) {
				ears[j] = (int) this.shellCoords.get(ears[j]).z;
			}
			final int[] E = extractEdgesTri(ears);
			return new WB_Triangulation2D(ears, E);
		}

		public int[] triangulateFace(final HE_Face face, final boolean optimize) {
			final int fo = face.getFaceDegree();
			final Coordinate[] coords = new Coordinate[fo + 1];
			final WB_Coord normal = HE_MeshOp.getFaceNormal(face);
			WB_Swizzle coordViewer;
			if (Math.abs(normal.xd()) > Math.abs(normal.yd())) {
				coordViewer = ((Math.abs(normal.xd()) > Math.abs(normal.zd())) ? WB_Swizzle.YZ : WB_Swizzle.XY);
			} else {
				coordViewer = ((Math.abs(normal.yd()) > Math.abs(normal.zd())) ? WB_Swizzle.ZX : WB_Swizzle.XY);
			}
			final WB_KDTreeInteger3D<WB_Point> pointmap = new WB_KDTreeInteger3D<WB_Point>();
			int i = 0;
			HE_Halfedge he = face.getHalfedge();
			do {
				coords[i] = new Coordinate(coordViewer.xd(he.getVertex()), coordViewer.yd(he.getVertex()), 0.0);
				pointmap.add(new WB_Point(coordViewer.xd(he.getVertex()), coordViewer.yd(he.getVertex())), i);
				++i;
				he = he.getNextInFace();
			} while (he != face.getHalfedge());
			coords[i] = new Coordinate(coordViewer.xd(he.getVertex()), coordViewer.yd(he.getVertex()), 0.0);
			final LinearRing shell = WB_JTS.JTSgf.createLinearRing(coords);
			final Polygon inputPolygon = WB_JTS.JTSgf.createPolygon(shell);
			final int[] ears = this.triangulate(inputPolygon, optimize);
			final List<WB_Point> tripoints = new FastList<WB_Point>();
			for (int j = 0; j < this.shellCoords.size() - 1; ++j) {
				tripoints.add(new WB_Point(this.shellCoords.get(j).x, this.shellCoords.get(j).y));
			}
			final int[] intmap = new int[tripoints.size()];
			i = 0;
			for (final WB_Coord point : tripoints) {
				final int found = pointmap.getNearestNeighbor(point).value;
				intmap[i++] = found;
			}
			for (int k = 0; k < ears.length; ++k) {
				ears[k] = intmap[ears[k]];
			}
			return ears;
		}

		public int[] triangulateSimplePolygon(final WB_Polygon polygon, final boolean optimize) {
			final int fo = polygon.getNumberOfShellPoints();
			final Coordinate[] coords = new Coordinate[fo + 1];
			final WB_Coord normal = polygon.getNormal();
			WB_Swizzle coordViewer;
			if (Math.abs(normal.xd()) > Math.abs(normal.yd())) {
				coordViewer = ((Math.abs(normal.xd()) > Math.abs(normal.zd())) ? WB_Swizzle.YZ : WB_Swizzle.XY);
			} else {
				coordViewer = ((Math.abs(normal.yd()) > Math.abs(normal.zd())) ? WB_Swizzle.ZX : WB_Swizzle.XY);
			}
			final WB_KDTreeInteger3D<WB_Point> pointmap = new WB_KDTreeInteger3D<WB_Point>();
			for (int i = 0; i <= polygon.getNumberOfShellPoints(); ++i, ++i) {
				final WB_Coord c = polygon.getPoint(i);
				coords[i] = new Coordinate(coordViewer.xd(c), coordViewer.yd(c), 0.0);
				pointmap.add(new WB_Point(coordViewer.xd(c), coordViewer.yd(c)), i);
			}
			final LinearRing shell = WB_JTS.JTSgf.createLinearRing(coords);
			final Polygon inputPolygon = WB_JTS.JTSgf.createPolygon(shell);
			final int[] ears = this.triangulate(inputPolygon, optimize);
			final List<WB_Point> tripoints = new FastList<WB_Point>();
			for (int j = 0; j < this.shellCoords.size() - 1; ++j) {
				tripoints.add(new WB_Point(this.shellCoords.get(j).x, this.shellCoords.get(j).y));
			}
			final int[] intmap = new int[tripoints.size()];
			int k = 0;
			for (final WB_Coord point : tripoints) {
				final int found = pointmap.getNearestNeighbor(point).value;
				intmap[k++] = found;
			}
			for (int l = 0; l < ears.length; ++l) {
				ears[l] = intmap[ears[l]];
			}
			return ears;
		}

		private static int[] extractEdgesTri(final int[] ears) {
			final int f = ears.length;
			final UnifiedMap<Long, int[]> map = new UnifiedMap<Long, int[]>();
			for (int i = 0; i < ears.length; i += 3) {
				final int v0 = ears[i];
				final int v2 = ears[i + 1];
				final int v3 = ears[i + 2];
				long index = getIndex(v0, v2, f);
				map.put(index, new int[] { v0, v2 });
				index = getIndex(v2, v3, f);
				map.put(index, new int[] { v2, v3 });
				index = getIndex(v3, v0, f);
				map.put(index, new int[] { v3, v0 });
			}
			final int[] edges = new int[2 * map.size()];
			final Collection<int[]> values = map.values();
			int j = 0;
			for (final int[] value : values) {
				edges[2 * j] = value[0];
				edges[2 * j + 1] = value[1];
				++j;
			}
			return edges;
		}

		private static long getIndex(final int i, final int j, final int f) {
			return (i > j) ? (j + i * f) : (i + j * f);
		}

		private int[] triangulate(final Polygon inputPolygon, final boolean improve) {
			final List<IndexedTriangle> earList = new ArrayList<>();
			this.createShell(inputPolygon);
			final Geometry test = inputPolygon.buffer(0.0);
			int N = this.shellCoords.size() - 1;
			Arrays.fill(this.shellCoordAvailable = new boolean[N], true);
			boolean finished = false;
			boolean found = false;
			int k0 = 0;
			int k2 = 1;
			int k3 = 2;
			int firstK = 0;
			do {
				found = false;
				while (CGAlgorithms.computeOrientation(this.shellCoords.get(k0), this.shellCoords.get(k2), this.shellCoords.get(k3)) == 0) {
					k0 = k2;
					if (k0 == firstK) {
						finished = true;
						break;
					}
					k2 = k3;
					k3 = this.nextShellCoord(k3 + 1);
				}
				if (!finished && this.isValidEdge(k0, k3)) {
					final LineString ls = WB_JTS.JTSgf
							.createLineString(new Coordinate[] { this.shellCoords.get(k0), this.shellCoords.get(k3) });
					if (test.covers(ls)) {
						final Polygon earPoly = WB_JTS.JTSgf.createPolygon(WB_JTS.JTSgf.createLinearRing(new Coordinate[] {
								this.shellCoords.get(k0), this.shellCoords.get(k2), this.shellCoords.get(k3), this.shellCoords.get(k0) }),
								null);
						if (test.covers(earPoly)) {
							found = true;
							final IndexedTriangle ear = new IndexedTriangle(k0, k2, k3);
							earList.add(ear);
							this.shellCoordAvailable[k2] = false;
							--N;
							k0 = this.nextShellCoord(0);
							k2 = this.nextShellCoord(k0 + 1);
							k3 = this.nextShellCoord(k2 + 1);
							firstK = k0;
							if (N < 3) {
								finished = true;
							}
						}
					}
				}
				if (!finished && !found) {
					k0 = k2;
					if (k0 == firstK) {
						finished = true;
					} else {
						k2 = k3;
						k3 = this.nextShellCoord(k3 + 1);
					}
				}
			} while (!finished);
			if (improve) {
				this.doImprove(earList);
			}
			final int[] tris = new int[3 * earList.size()];
			for (int i = 0; i < earList.size(); ++i) {
				final int[] tri = earList.get(i).getVertices();
				tris[3 * i] = tri[0];
				tris[3 * i + 1] = tri[1];
				tris[3 * i + 2] = tri[2];
			}
			return tris;
		}

		protected WB_Polygon makeSimplePolygon(final WB_Polygon polygon) {
			final Polygon poly = WB_JTS.toJTSPolygon2D(polygon);
			this.createShell(poly);
			final Coordinate[] coords = new Coordinate[this.shellCoords.size()];
			return WB_JTS.createPolygonFromJTSPolygon2D(WB_JTS.JTSgf.createPolygon(this.shellCoords.toArray(coords)));
		}

		protected void createShell(final Polygon inputPolygon) {
			final Polygon poly = (Polygon) inputPolygon.clone();
			this.shellCoords = new ArrayList<>();
			final List<Geometry> orderedHoles = getOrderedHoles(poly);
			final Coordinate[] coords = poly.getExteriorRing().getCoordinates();
			this.shellCoords.addAll(Arrays.asList(coords));
			for (Geometry orderedHole : orderedHoles) {
				this.joinHoleToShell(inputPolygon, orderedHole);
			}
		}

		private boolean isValidEdge(final int index0, final int index1) {
			final Coordinate[] line = { this.shellCoords.get(index0), this.shellCoords.get(index1) };
			for (int index2 = this.nextShellCoord(index0 + 1); index2 != index0; index2 = this.nextShellCoord(index2 + 1)) {
				if (index2 != index1) {
					final Coordinate c = this.shellCoords.get(index2);
					if (!c.equals2D(line[0]) && !c.equals2D(line[1]) && CGAlgorithms.isOnLine(c, line)) {
						return false;
					}
				}
			}
			return true;
		}

		private int nextShellCoord(final int pos) {
			int pnew;
			for (pnew = pos % this.shellCoordAvailable.length; !this.shellCoordAvailable[pnew]; pnew = (pnew + 1)
					% this.shellCoordAvailable.length) {
			}
			return pnew;
		}

		private void doImprove(final List<IndexedTriangle> earList) {
			final EdgeFlipper ef = new EdgeFlipper(this.shellCoords);
			boolean changed;
			do {
				changed = false;
				for (int i = 0; i < earList.size() - 1 && !changed; ++i) {
					final IndexedTriangle ear0 = earList.get(i);
					for (int j = i + 1; j < earList.size() && !changed; ++j) {
						final IndexedTriangle ear2 = earList.get(j);
						final int[] sharedVertices = ear0.getSharedVertices(ear2);
						if (sharedVertices != null && sharedVertices.length == 2 && ef.flip(ear0, ear2, sharedVertices)) {
							changed = true;
						}
					}
				}
			} while (changed);
		}

		private static List<Geometry> getOrderedHoles(final Polygon poly) {
			final List<Geometry> holes = new ArrayList<>();
			final List<IndexedEnvelope> bounds = new ArrayList<>();
			if (poly.getNumInteriorRing() > 0) {
				for (int i = 0; i < poly.getNumInteriorRing(); ++i) {
					bounds.add(new IndexedEnvelope(i, poly.getInteriorRingN(i).getEnvelopeInternal()));
				}
				Collections.sort(bounds, new IndexedEnvelopeComparator());
				for (int i = 0; i < bounds.size(); ++i) {
					holes.add(poly.getInteriorRingN(bounds.get(i).index));
				}
			}
			return holes;
		}

		private void joinHoleToShell(final Polygon inputPolygon, final Geometry hole) {
			double minD2 = Double.MAX_VALUE;
			int shellVertexIndex = -1;
			final int Ns = this.shellCoords.size() - 1;
			final int holeVertexIndex = getLowestVertex(hole);
			final Coordinate[] holeCoords = hole.getCoordinates();
			final Coordinate ch = holeCoords[holeVertexIndex];
			final List<IndexedDouble> distanceList = new ArrayList<>();
			for (int i = Ns - 1; i >= 0; --i) {
				final Coordinate cs = this.shellCoords.get(i);
				final double d2 = (ch.x - cs.x) * (ch.x - cs.x) + (ch.y - cs.y) * (ch.y - cs.y);
				if (d2 < minD2) {
					minD2 = d2;
					shellVertexIndex = i;
				}
				distanceList.add(new IndexedDouble(i, d2));
			}
			LineString join = WB_JTS.JTSgf.createLineString(new Coordinate[] { ch, this.shellCoords.get(shellVertexIndex) });
			if (inputPolygon.covers(join)) {
				this.doJoinHole(shellVertexIndex, holeCoords, holeVertexIndex);
				return;
			}
			Collections.sort(distanceList, new IndexedDoubleComparator());
			for (int j = 1; j < distanceList.size(); ++j) {
				join = WB_JTS.JTSgf.createLineString(new Coordinate[] { ch, this.shellCoords.get(distanceList.get(j).index) });
				if (inputPolygon.covers(join)) {
					shellVertexIndex = distanceList.get(j).index;
					this.doJoinHole(shellVertexIndex, holeCoords, holeVertexIndex);
					return;
				}
			}
		}

		private void doJoinHole(final int shellVertexIndex, final Coordinate[] holeCoords, final int holeVertexIndex) {
			final List<Coordinate> newCoords = new ArrayList<>();
			newCoords.add(new Coordinate(this.shellCoords.get(shellVertexIndex)));
			final int N = holeCoords.length - 1;
			int i = holeVertexIndex;
			do {
				newCoords.add(new Coordinate(holeCoords[i]));
				i = (i + 1) % N;
			} while (i != holeVertexIndex);
			newCoords.add(new Coordinate(holeCoords[holeVertexIndex]));
			this.shellCoords.addAll(shellVertexIndex, newCoords);
		}

		private static int getLowestVertex(final Geometry geom) {
			final Coordinate[] coords = geom.getCoordinates();
			final double minY = geom.getEnvelopeInternal().getMinY();
			for (int i = 0; i < coords.length; ++i) {
				if (Math.abs(coords[i].y - minY) < WB_Epsilon.EPSILON) {
					return i;
				}
			}
			throw new IllegalStateException("Failed to find lowest vertex");
		}

		private static class IndexedDoubleComparator implements Comparator<IndexedDouble> {
			@Override
			public int compare(final IndexedDouble o1, final IndexedDouble o2) {
				final double delta = o1.value - o2.value;
				if (Math.abs(delta) < WB_Epsilon.EPSILON) {
					return 0;
				}
				return (delta > 0.0) ? 1 : -1;
			}
		}

		private static class IndexedEnvelopeComparator implements Comparator<IndexedEnvelope> {
			@Override
			public int compare(final IndexedEnvelope o1, final IndexedEnvelope o2) {
				double delta = o1.envelope.getMinY() - o2.envelope.getMinY();
				if (Math.abs(delta) < WB_Epsilon.EPSILON) {
					delta = o1.envelope.getMinX() - o2.envelope.getMinX();
					if (Math.abs(delta) < WB_Epsilon.EPSILON) {
						return 0;
					}
				}
				return (delta > 0.0) ? 1 : -1;
			}
		}

		private static class EdgeFlipper {
			private final List<Coordinate> shellCoords;

			EdgeFlipper(final List<Coordinate> shellCoords) {
				this.shellCoords = Collections.unmodifiableList((List<? extends Coordinate>) shellCoords);
			}

			boolean flip(final IndexedTriangle ear0, final IndexedTriangle ear1, final int[] sharedVertices) {
				if (sharedVertices == null || sharedVertices.length != 2) {
					return false;
				}
				final Coordinate shared0 = this.shellCoords.get(sharedVertices[0]);
				final Coordinate shared2 = this.shellCoords.get(sharedVertices[1]);
				int[] vertices;
				int i;
				for (vertices = ear0.getVertices(), i = 0; vertices[i] == sharedVertices[0] || vertices[i] == sharedVertices[1]; ++i) {
				}
				final int v0 = vertices[i];
				boolean reverse = false;
				if (vertices[(i + 1) % 3] == sharedVertices[0]) {
					reverse = true;
				}
				final Coordinate c0 = this.shellCoords.get(v0);
				for (i = 0, vertices = ear1.getVertices(); vertices[i] == sharedVertices[0] || vertices[i] == sharedVertices[1]; ++i) {
				}
				final int v2 = vertices[i];
				final Coordinate c2 = this.shellCoords.get(v2);
				final int dir0 = CGAlgorithms.orientationIndex(c0, c2, shared0);
				final int dir2 = CGAlgorithms.orientationIndex(c0, c2, shared2);
				if (dir0 == -dir2 && c0.distance(c2) < shared0.distance(shared2)) {
					if (reverse) {
						ear0.setPoints(sharedVertices[0], v2, v0);
						ear1.setPoints(v0, v2, sharedVertices[1]);
					} else {
						ear0.setPoints(sharedVertices[0], v0, v2);
						ear1.setPoints(v2, v0, sharedVertices[1]);
					}
					return true;
				}
				return false;
			}
		}

		private static class IndexedTriangle {
			private final int[] points;

			IndexedTriangle(final int v0, final int v1, final int v2) {
				this.points = new int[3];
				this.setPoints(v0, v1, v2);
			}

			void setPoints(final int v0, final int v1, final int v2) {
				this.points[0] = v0;
				this.points[1] = v1;
				this.points[2] = v2;
			}

			int[] getVertices() {
				final int[] copy = new int[3];
				for (int i = 0; i < 3; ++i) {
					copy[i] = this.points[i];
				}
				return copy;
			}

			int[] getSharedVertices(final IndexedTriangle other) {
				int count = 0;
				final boolean[] shared = new boolean[3];
				for (int i = 0; i < 3; ++i) {
					for (int j = 0; j < 3; ++j) {
						if (this.points[i] == other.points[j]) {
							++count;
							shared[i] = true;
						}
					}
				}
				int[] common = null;
				if (count > 0) {
					common = new int[count];
					int k = 0;
					int l = 0;
					while (k < 3) {
						if (shared[k]) {
							common[l++] = this.points[k];
						}
						++k;
					}
				}
				return common;
			}
		}

		private static class IndexedEnvelope {
			int index;
			Envelope envelope;

			IndexedEnvelope(final int i, final Envelope env) {
				this.index = i;
				this.envelope = env;
			}
		}

		private static class IndexedDouble {
			int index;
			double value;

			IndexedDouble(final int i, final double v) {
				this.index = i;
				this.value = v;
			}
		}
	}

	static class ShapeReader {
		private static AffineTransform INVERT_Y;
		private WB_GeometryFactory2D geometryfactory;

		static {
			ShapeReader.INVERT_Y = AffineTransform.getScaleInstance(1.0, -1.0);
		}

		ShapeReader() {
			this.geometryfactory = new WB_GeometryFactory2D();
		}

		List<WB_Polygon> read(final Shape shp, final double flatness) {
			final PathIterator pathIt = shp.getPathIterator(ShapeReader.INVERT_Y, flatness);
			return this.read(pathIt);
		}

		List<WB_Polygon> read(final PathIterator pathIt) {
			final List<Coordinate[]> pathPtSeq = toCoordinates(pathIt);
			final RingTree tree = new RingTree();
			for (final Coordinate[] pts : pathPtSeq) {
				final LinearRing ring = WB_JTS.JTSgf.createLinearRing(pts);
				tree.add(ring);
			}
			return tree.extractPolygons();
		}

		private static List<Coordinate[]> toCoordinates(final PathIterator pathIt) {
			final List<Coordinate[]> coordArrays = new ArrayList<>();
			while (!pathIt.isDone()) {
				final Coordinate[] pts = nextCoordinateArray(pathIt);
				if (pts == null) {
					break;
				}
				coordArrays.add(pts);
			}
			return coordArrays;
		}

		private static Coordinate[] nextCoordinateArray(final PathIterator pathIt) {
			final double[] pathPt = new double[6];
			CoordinateList coordList = null;
			boolean isDone = false;
			while (!pathIt.isDone()) {
				final int segType = pathIt.currentSegment(pathPt);
				switch (segType) {
					case 0 : {
						if (coordList != null) {
							isDone = true;
							break;
						}
						coordList = new CoordinateList();
						coordList.add(new Coordinate(pathPt[0], pathPt[1]));
						pathIt.next();
						break;
					}
					case 1 : {
						coordList.add(new Coordinate(pathPt[0], pathPt[1]));
						pathIt.next();
						break;
					}
					case 4 : {
						coordList.closeRing();
						pathIt.next();
						isDone = true;
						break;
					}
					default : {
						throw new IllegalArgumentException("unhandled (non-linear) segment type encountered");
					}
				}
				if (isDone) {
					break;
				}
			}
			return coordList.toCoordinateArray();
		}

		private static class RingNode {
			RingNode parent;
			List<RingNode> children;
			LinearRing ring;
			Polygon poly;
			boolean hole;

			RingNode() {
				this.parent = null;
				this.ring = null;
				this.children = new ArrayList<>();
				this.hole = true;
			}

			RingNode(final RingNode parent, final LinearRing ring) {
				this.parent = parent;
				this.ring = ring;
				final Coordinate[] coords = ring.getCoordinates();
				this.poly = WB_JTS.JTSgf.createPolygon(coords);
				this.hole = CGAlgorithms.isCCW(coords);
				this.children = new ArrayList<>();
			}
		}

		private class RingTree {
			RingNode root;

			RingTree() {
				this.root = new RingNode();
			}

			void add(final LinearRing ring) {
				final Polygon poly = WB_JTS.JTSgf.createPolygon(ring);
				RingNode currentParent = this.root;
				RingNode foundParent;
				do {
					foundParent = null;
					for (final RingNode node : currentParent.children) {
						final Polygon other = node.poly;
						if (poly.within(other)) {
							foundParent = node;
							currentParent = node;
							break;
						}
					}
				} while (foundParent != null);
				final RingNode newNode = new RingNode(currentParent, ring);
				final List<RingNode> nodesToRemove = new ArrayList<>();
				for (int j = 0; j < currentParent.children.size(); ++j) {
					final RingNode node2 = currentParent.children.get(j);
					final Polygon other2 = node2.poly;
					if (other2.within(poly)) {
						newNode.children.add(node2);
						nodesToRemove.add(node2);
					}
				}
				currentParent.children.removeAll(nodesToRemove);
				currentParent.children.add(newNode);
			}

			List<WB_Polygon> extractPolygons() {
				final List<WB_Polygon> polygons = new ArrayList<>();
				final List<RingNode> shellNodes = new ArrayList<>();
				this.addExteriorNodes(this.root, shellNodes);
				for (final RingNode node : shellNodes) {
					int count = 0;
					for (RingNode element : node.children) {
						if (element.hole) {
							++count;
						}
					}
					final LinearRing[] holes = new LinearRing[count];
					int index = 0;
					for (RingNode element : node.children) {
						if (element.hole) {
							holes[index++] = element.ring;
						}
					}
					final Geometry result = WB_JTS.JTSgf.createPolygon(node.ring, holes);
					if (result.getGeometryType().equals("Polygon")) {
						polygons.add(WB_JTS.createPolygonFromJTSPolygon2D((Polygon) result));
					} else {
						if (!result.getGeometryType().equals("MultiPolygon")) {
							continue;
						}
						for (int k = 0; k < result.getNumGeometries(); ++k) {
							final Geometry ggeo = result.getGeometryN(k);
							polygons.add(WB_JTS.createPolygonFromJTSPolygon2D((Polygon) ggeo));
						}
					}
				}
				return polygons;
			}

			void addExteriorNodes(final RingNode parent, final List<RingNode> shellNodes) {
				for (final RingNode node : parent.children) {
					if (!node.hole) {
						shellNodes.add(node);
					}
					this.addExteriorNodes(node, shellNodes);
				}
			}
		}
	}
}
