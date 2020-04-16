package wblut.processing;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import processing.core.PApplet;
import processing.core.PConstants;
import processing.core.PImage;
import processing.core.PMatrix3D;
import processing.opengl.PGraphicsOpenGL;
import wblut.core.WB_ProgressReporter.WB_ProgressTracker;
import wblut.geom.WB_AABB;
import wblut.geom.WB_AABBTree3D;
import wblut.geom.WB_AABBTree3D.WB_AABBNode3D;
import wblut.geom.WB_BinaryGrid3D;
import wblut.geom.WB_Circle;
import wblut.geom.WB_Classification;
import wblut.geom.WB_Coord;
import wblut.geom.WB_CoordCollection;
import wblut.geom.WB_Curve;
import wblut.geom.WB_Danzer3D;
import wblut.geom.WB_Danzer3D.WB_DanzerTile3D;
import wblut.geom.WB_GeometryFactory3D;
import wblut.geom.WB_GeometryOp3D;
import wblut.geom.WB_Hexagon;
import wblut.geom.WB_IndexedAABBTree3D;
import wblut.geom.WB_IndexedAABBTree3D.WB_IndexedAABBNode3D;
import wblut.geom.WB_Line;
import wblut.geom.WB_Map;
import wblut.geom.WB_Map2D;
import wblut.geom.WB_Network;
import wblut.geom.WB_Network.Connection;
import wblut.geom.WB_Network.Node;
import wblut.geom.WB_Octagon;
import wblut.geom.WB_OctreeInteger;
import wblut.geom.WB_Pentagon;
import wblut.geom.WB_Plane;
import wblut.geom.WB_Point;
import wblut.geom.WB_PolyLine;
import wblut.geom.WB_Polygon;
import wblut.geom.WB_Quad;
import wblut.geom.WB_QuadtreeInteger;
import wblut.geom.WB_Ray;
import wblut.geom.WB_Ring;
import wblut.geom.WB_Segment;
import wblut.geom.WB_SimpleMesh;
import wblut.geom.WB_Tetrahedron;
import wblut.geom.WB_Transform3D;
import wblut.geom.WB_Triangle;
import wblut.geom.WB_TriangleFactory;
import wblut.geom.WB_Triangulation2D;
import wblut.geom.WB_Triangulation3D;
import wblut.geom.WB_Vector;
import wblut.hemesh.HE_Face;
import wblut.hemesh.HE_FaceEdgeCirculator;
import wblut.hemesh.HE_FaceHalfedgeInnerCirculator;
import wblut.hemesh.HE_FaceIterator;
import wblut.hemesh.HE_FaceVertexCirculator;
import wblut.hemesh.HE_Halfedge;
import wblut.hemesh.HE_HalfedgeStructure;
import wblut.hemesh.HE_Mesh;
import wblut.hemesh.HE_MeshCollection;
import wblut.hemesh.HE_MeshIterator;
import wblut.hemesh.HE_MeshOp;
import wblut.hemesh.HE_MeshOp.HE_FaceLineIntersection;
import wblut.hemesh.HE_Path;
import wblut.hemesh.HE_Selection;
import wblut.hemesh.HE_TextureCoordinate;
import wblut.hemesh.HE_Vertex;
import wblut.hemesh.HE_VertexIterator;

public class WB_Render3D extends WB_Render2D {
	private final WB_GeometryFactory3D geometryfactory = new WB_GeometryFactory3D();
	public static final WB_ProgressTracker tracker = WB_ProgressTracker.instance();
	protected final PGraphicsOpenGL home;

	public WB_Render3D(final PApplet home) {
		super(home);
		if (home.g == null) {
			throw new IllegalArgumentException("WB_Render3D can only be used after size()");
		}
		if (!(home.g instanceof PGraphicsOpenGL)) {
			throw new IllegalArgumentException(
					"WB_Render3D can only be used with P3D, OPENGL or derived ProcessingPGraphics object");
		}
		this.home = (PGraphicsOpenGL) home.g;
	}

	public WB_Render3D(final PGraphicsOpenGL home) {
		super(home);
		this.home = home;
	}

	public void drawAABB(final WB_AABB AABB) {
		home.pushMatrix();
		translate(AABB.getCenter());
		home.box((float) AABB.getWidth(), (float) AABB.getHeight(), (float) AABB.getDepth());
		home.popMatrix();
	}

	public void drawAABBNode(final WB_AABBNode3D node) {
		drawAABB(node.getAABB());
		if (node.getChildA() != null) {
			drawAABBNode(node.getChildA());
		}
		if (node.getChildB() != null) {
			drawAABBNode(node.getChildB());
		}
	}

	private void drawAABBNode(final WB_AABBNode3D node, final int level) {
		if (node.getLevel() == level) {
			drawAABB(node.getAABB());
		}
		if (node.getLevel() < level) {
			if (node.getChildA() != null) {
				drawAABBNode(node.getChildA(), level);
			}
			if (node.getChildB() != null) {
				drawAABBNode(node.getChildB(), level);
			}
		}
	}

	private void drawAABBLeafNode(final WB_AABBNode3D node) {
		if (node.isLeaf()) {
			drawAABB(node.getAABB());
		} else {
			if (node.getChildA() != null) {
				drawAABBLeafNode(node.getChildA());
			}
			if (node.getChildB() != null) {
				drawAABBLeafNode(node.getChildB());
			}
		}
	}

	public void drawAABBLeafNodes(final WB_AABBTree3D tree) {
		drawAABBLeafNode(tree.getRoot());
	}

	public void drawAABBNode(final WB_IndexedAABBNode3D node) {
		drawAABB(node.getAABB());
		if (node.getChildA() != null) {
			drawAABBNode(node.getChildA());
		}
		if (node.getChildB() != null) {
			drawAABBNode(node.getChildB());
		}
	}

	private void drawAABBNode(final WB_IndexedAABBNode3D node, final int level) {
		if (node.getLevel() == level) {
			drawAABB(node.getAABB());
		}
		if (node.getLevel() < level) {
			if (node.getChildA() != null) {
				drawAABBNode(node.getChildA(), level);
			}
			if (node.getChildB() != null) {
				drawAABBNode(node.getChildB(), level);
			}
		}
	}

	private void drawAABBLeafNode(final WB_IndexedAABBNode3D node) {
		if (node.isLeaf()) {
			drawAABB(node.getAABB());
		} else {
			if (node.getChildA() != null) {
				drawAABBLeafNode(node.getChildA());
			}
			if (node.getChildB() != null) {
				drawAABBLeafNode(node.getChildB());
			}
		}
	}

	public void drawAABBLeafNodes(final WB_IndexedAABBTree3D tree) {
		drawAABBLeafNode(tree.getRoot());
	}

	public void drawAABBTree(final WB_AABBTree3D tree) {
		drawAABBNode(tree.getRoot());
	}

	public void drawAABBTree(final WB_AABBTree3D tree, final int level) {
		drawAABBNode(tree.getRoot(), level);
	}

	public void drawAABBTree(final WB_IndexedAABBTree3D tree) {
		drawAABBNode(tree.getRoot());
	}

	public void drawAABBTree(final WB_IndexedAABBTree3D tree, final int level) {
		drawAABBNode(tree.getRoot(), level);
	}

	public void drawBezierEdges(final HE_HalfedgeStructure mesh) {
		HE_Halfedge he;
		WB_Coord p0;
		WB_Coord p1;
		WB_Coord p2;
		WB_Coord p3;
		HE_Face f;
		final HE_FaceIterator fItr = mesh.fItr();
		while (fItr.hasNext()) {
			f = fItr.next();
			if (f.isVisible()) {
				home.beginShape();
				he = f.getHalfedge();
				p0 = HE_MeshOp.getHalfedgeCenter(he.getPrevInFace());
				vertex(p0);
				do {
					p1 = he.getVertex();
					p2 = he.getVertex();
					p3 = HE_MeshOp.getHalfedgeCenter(he);
					home.bezierVertex(p1.xf(), p1.yf(), p1.zf(), p2.xf(), p2.yf(), p2.zf(), p3.xf(), p3.yf(), p3.zf());
					he = he.getNextInFace();
				} while (he != f.getHalfedge());
				home.endShape();
			}
		}
	}

	public void drawBoundaryEdges(final HE_HalfedgeStructure mesh) {
		HE_Halfedge he;
		final Iterator<HE_Halfedge> heItr = mesh.heItr();
		while (heItr.hasNext()) {
			he = heItr.next();
			if (he.isVisible()) {
				if (he.getFace() == null) {
					if (he.getPair() != null) {
						line(he.getVertex(), he.getPair().getVertex());
					} else if (he.getNextInFace() != null) {
						line(he.getVertex(), he.getNextInFace().getVertex());
					}
				}
			}
		}
	}

	public void drawBoundaryVertices(final HE_HalfedgeStructure mesh) {
		HE_Vertex v;
		final Iterator<HE_Vertex> vItr = mesh.vItr();
		while (vItr.hasNext()) {
			v = vItr.next();
			if (v.isVisible()) {
				if (v.isBoundary()) {
					drawPoint(v);
				}
			}
		}
	}

	public void drawBoundaryVertices(final HE_HalfedgeStructure mesh, final double r) {
		HE_Vertex v;
		final Iterator<HE_Vertex> vItr = mesh.vItr();
		while (vItr.hasNext()) {
			v = vItr.next();
			if (v.isVisible()) {
				if (v.isBoundary()) {
					drawVertex(v, r);
				}
			}
		}
	}

	public void drawBoundaryFaces(final HE_HalfedgeStructure mesh) {
		HE_Face f;
		final HE_FaceIterator fItr = mesh.fItr();
		while (fItr.hasNext()) {
			f = fItr.next();
			if (f.isVisible()) {
				if (f.isBoundary()) {
					drawFace(f);
				}
			}
		}
	}

	public void drawBoundaryHalfedges(final HE_HalfedgeStructure mesh) {
		drawBoundaryEdges(mesh);
	}

	public void drawCircle(final WB_Circle C) {
		home.pushMatrix();
		translate(C.getCenter());
		final WB_Transform3D T = new WB_Transform3D(geometryfactory.Z(), C.getNormal());
		final WB_Vector angles = T.getEulerAnglesXYZ();
		home.rotateZ(angles.zf());
		home.rotateY(angles.yf());
		home.rotateX(angles.xf());
		home.ellipse(0, 0, 2 * (float) C.getRadius(), 2 * (float) C.getRadius());
		home.popMatrix();
	}

	public void drawCircle(final Collection<WB_Circle> circles) {
		final Iterator<WB_Circle> citr = circles.iterator();
		while (citr.hasNext()) {
			drawCircle(citr.next());
		}
	}

	public void drawCircle(final WB_Circle[] circles) {
		for (final WB_Circle circle : circles) {
			drawCircle(circle);
		}
	}

	public void drawCurve(final Collection<WB_Curve> curves, final int steps) {
		final Iterator<WB_Curve> citr = curves.iterator();
		while (citr.hasNext()) {
			drawCurve(citr.next(), steps);
		}
	}

	public void drawCurve(final WB_Curve C, final double minU, final double maxU, final int steps) {
		final int n = Math.max(1, steps);
		WB_Point p0 = C.getPointOnCurve(minU);
		WB_Point p1;
		final double du = (maxU - minU) / n;
		for (int i = 0; i < n; i++) {
			p1 = C.getPointOnCurve(minU + (i + 1) * du);
			line(p0, p1);
			p0 = p1;
		}
	}

	public void drawCurve(final WB_Curve C, final int steps) {
		final int n = Math.max(1, steps);
		WB_Point p0 = C.getPointOnCurve(0);
		WB_Point p1;
		final double du = 1.0 / n;
		for (int i = 0; i < n; i++) {
			p1 = C.getPointOnCurve((i + 1) * du);
			line(p0, p1);
			p0 = p1;
		}
	}

	public void drawEdge(final HE_Halfedge e) {
		if (e.isVisible()) {
			line(e.getStartVertex(), e.getEndVertex());
		}
	}

	public void drawEdge(final long key, final HE_Mesh mesh) {
		final HE_Halfedge e = mesh.getHalfedgeWithKey(key);
		if (e != null && e.isVisible()) {
			drawEdge(e);
		}
	}

	public void drawEdges(final Collection<? extends HE_HalfedgeStructure> meshes) {
		final Iterator<? extends HE_HalfedgeStructure> mItr = meshes.iterator();
		while (mItr.hasNext()) {
			drawEdges(mItr.next());
		}
	}

	public void drawEdges(final HE_MeshCollection meshes) {
		final Iterator<? extends HE_HalfedgeStructure> mItr = meshes.mItr();
		while (mItr.hasNext()) {
			drawEdges(mItr.next());
		}
	}

	public void drawEdges(final HE_Face f) {
		HE_Halfedge e = f.getHalfedge();
		do {
			if (e.isVisible()) {
				line(e.getVertex(), e.getEndVertex());
			}
			e = e.getNextInFace();
		} while (e != f.getHalfedge());
	}

	public void drawEdges(final HE_HalfedgeStructure mesh) {
		if (mesh == null) {
			return;
		}
		final Iterator<HE_Halfedge> eItr = mesh.eItr();
		HE_Halfedge e;
		while (eItr.hasNext()) {
			e = eItr.next();
			if (e.isVisible()) {
				line(e.getVertex(), e.getEndVertex());
			}
		}
	}

	public void drawEdges(final HE_Selection selection) {
		if (selection == null) {
			return;
		}
		final Iterator<HE_Halfedge> heItr = selection.heItr();
		HE_Halfedge he;
		while (heItr.hasNext()) {
			he = heItr.next();
			if (he.isEdge() || !selection.contains(he.getPair())) {
				if (he.isVisible()) {
					line(he.getVertex(), he.getEndVertex());
				}
			}
		}
	}

	public void drawEdgesWithInternalLabel(final int label, final HE_HalfedgeStructure mesh) {
		if (mesh == null) {
			return;
		}
		final Iterator<HE_Halfedge> eItr = mesh.eItr();
		HE_Halfedge e;
		while (eItr.hasNext()) {
			e = eItr.next();
			if (e.isVisible()) {
				if (e.getInternalLabel() == label) {
					line(e.getVertex(), e.getEndVertex());
				}
			}
		}
	}

	public void drawEdgesWithLabel(final int label, final HE_HalfedgeStructure mesh) {
		if (mesh == null) {
			return;
		}
		final Iterator<HE_Halfedge> eItr = mesh.eItr();
		HE_Halfedge e;
		while (eItr.hasNext()) {
			e = eItr.next();
			if (e.isVisible()) {
				if (e.getLabel() == label) {
					line(e.getVertex(), e.getEndVertex());
				}
			}
		}
	}

	public void drawFace(final HE_Face f) {
		drawFace(f, false);
	}

	private void drawFaceInt(final HE_Face f) {
		drawFaceInt(f, false);
	}

	public void drawFace(final HE_Face f, final boolean smooth) {
		if (!f.isVisible()) {
			return;
		}
		final int fo = f.getFaceDegree();
		final List<HE_Vertex> vertices = f.getFaceVertices();
		if (fo < 3 || vertices.size() < 3) {
		} else if (fo == 3) {
			final int[] tri = new int[] { 0, 1, 2 };
			HE_Vertex v0, v1, v2;
			WB_Coord n0, n1, n2;
			if (smooth) {
				home.beginShape(PConstants.TRIANGLES);
				v0 = vertices.get(tri[0]);
				n0 = HE_MeshOp.getVertexNormal(v0);
				v1 = vertices.get(tri[1]);
				n1 = HE_MeshOp.getVertexNormal(v1);
				v2 = vertices.get(tri[2]);
				n2 = HE_MeshOp.getVertexNormal(v2);
				normal(n0);
				vertex(v0);
				normal(n1);
				vertex(v1);
				normal(n2);
				vertex(v2);
				home.endShape();
			} else {
				home.beginShape(PConstants.TRIANGLES);
				v0 = vertices.get(tri[0]);
				v1 = vertices.get(tri[1]);
				v2 = vertices.get(tri[2]);
				vertex(v0);
				vertex(v1);
				vertex(v2);
				home.endShape();
			}
		} else {
			final int[] tris = f.getTriangles();
			HE_Vertex v0, v1, v2;
			WB_Coord n0, n1, n2;
			if (smooth) {
				for (int i = 0; i < tris.length; i += 3) {
					home.beginShape(PConstants.TRIANGLES);
					v0 = vertices.get(tris[i]);
					n0 = HE_MeshOp.getVertexNormal(v0);
					v1 = vertices.get(tris[i + 1]);
					n1 = HE_MeshOp.getVertexNormal(v1);
					v2 = vertices.get(tris[i + 2]);
					n2 = HE_MeshOp.getVertexNormal(v2);
					normal(n0);
					vertex(v0);
					normal(n1);
					vertex(v1);
					normal(n2);
					vertex(v2);
					home.endShape();
				}
			} else {
				for (int i = 0; i < tris.length; i += 3) {
					home.beginShape(PConstants.TRIANGLES);
					v0 = vertices.get(tris[i]);
					v1 = vertices.get(tris[i + 1]);
					v2 = vertices.get(tris[i + 2]);
					vertex(v0);
					vertex(v1);
					vertex(v2);
					home.endShape();
				}
			}
		}
	}

	private void drawFaceInt(final HE_Face f, final boolean smooth) {
		if (!f.isVisible()) {
			return;
		}
		final int fo = f.getFaceDegree();
		final List<HE_Vertex> vertices = f.getFaceVertices();
		if (fo < 3 || vertices.size() < 3) {
		} else if (fo == 3) {
			final int[] tri = new int[] { 0, 1, 2 };
			HE_Vertex v0, v1, v2;
			WB_Coord n0, n1, n2;
			if (smooth) {
				v0 = vertices.get(tri[0]);
				n0 = HE_MeshOp.getVertexNormal(v0);
				v1 = vertices.get(tri[1]);
				n1 = HE_MeshOp.getVertexNormal(v1);
				v2 = vertices.get(tri[2]);
				n2 = HE_MeshOp.getVertexNormal(v2);
				normal(n0);
				vertex(v0);
				normal(n1);
				vertex(v1);
				normal(n2);
				vertex(v2);
			} else {
				v0 = vertices.get(tri[0]);
				v1 = vertices.get(tri[1]);
				v2 = vertices.get(tri[2]);
				vertex(v0);
				vertex(v1);
				vertex(v2);
			}
		} else {
			final int[] tris = f.getTriangles();
			HE_Vertex v0, v1, v2;
			WB_Coord n0, n1, n2;
			if (smooth) {
				for (int i = 0; i < tris.length; i += 3) {
					v0 = vertices.get(tris[i]);
					n0 = HE_MeshOp.getVertexNormal(v0);
					v1 = vertices.get(tris[i + 1]);
					n1 = HE_MeshOp.getVertexNormal(v1);
					v2 = vertices.get(tris[i + 2]);
					n2 = HE_MeshOp.getVertexNormal(v2);
					normal(n0);
					vertex(v0);
					normal(n1);
					vertex(v1);
					normal(n2);
					vertex(v2);
				}
			} else {
				for (int i = 0; i < tris.length; i += 3) {
					v0 = vertices.get(tris[i]);
					v1 = vertices.get(tris[i + 1]);
					v2 = vertices.get(tris[i + 2]);
					vertex(v0);
					vertex(v1);
					vertex(v2);
				}
			}
		}
	}

	public void drawFace(final HE_Face f, final PImage texture) {
		drawFace(f, texture, false);
	}

	private void drawFaceInt(final HE_Face f, final PImage texture) {
		drawFaceInt(f, texture, false);
	}

	public void drawFace(final HE_Face f, final PImage texture, final boolean smooth) {
		if (!f.isVisible()) {
			return;
		}
		final int fo = f.getFaceDegree();
		final List<HE_Halfedge> halfedges = f.getFaceHalfedges();
		if (fo < 3 || halfedges.size() < 3) {
		} else if (fo == 3) {
			final int[] tri = new int[] { 0, 1, 2 };
			HE_Halfedge he0, he1, he2;
			HE_Vertex v0, v1, v2;
			WB_Coord n0, n1, n2;
			he0 = halfedges.get(tri[0]);
			v0 = he0.getVertex();
			n0 = HE_MeshOp.getVertexNormal(v0);
			he1 = halfedges.get(tri[1]);
			v1 = he1.getVertex();
			n1 = HE_MeshOp.getVertexNormal(v1);
			he2 = halfedges.get(tri[2]);
			v2 = he2.getVertex();
			n2 = HE_MeshOp.getVertexNormal(v2);
			if (smooth) {
				home.beginShape(PConstants.TRIANGLES);
				home.texture(texture);
				normal(n0);
				home.vertex(v0.xf(), v0.yf(), v0.zf(), getU(he0.getUVW()), getV(he0.getUVW()));
				normal(n1);
				home.vertex(v1.xf(), v1.yf(), v1.zf(), getU(he1.getUVW()), getV(he1.getUVW()));
				normal(n2);
				home.vertex(v2.xf(), v2.yf(), v2.zf(), getU(he2.getUVW()), getV(he2.getUVW()));
				home.endShape();
			} else {
				home.beginShape(PConstants.TRIANGLES);
				home.texture(texture);
				home.vertex(v0.xf(), v0.yf(), v0.zf(), getU(he0.getUVW()), getV(he0.getUVW()));
				home.vertex(v1.xf(), v1.yf(), v1.zf(), getU(he1.getUVW()), getV(he1.getUVW()));
				home.vertex(v2.xf(), v2.yf(), v2.zf(), getU(he2.getUVW()), getV(he2.getUVW()));
				home.endShape();
			}
		} else {
			final int[] tris = f.getTriangles();
			HE_Halfedge he0, he1, he2;
			HE_Vertex v0, v1, v2;
			WB_Coord n0, n1, n2;
			if (smooth) {
				for (int i = 0; i < tris.length; i += 3) {
					home.beginShape(PConstants.TRIANGLES);
					home.texture(texture);
					he0 = halfedges.get(tris[i]);
					v0 = he0.getVertex();
					n0 = HE_MeshOp.getVertexNormal(v0);
					he1 = halfedges.get(tris[i + 1]);
					v1 = he1.getVertex();
					n1 = HE_MeshOp.getVertexNormal(v1);
					he2 = halfedges.get(tris[i + 2]);
					v2 = he2.getVertex();
					n2 = HE_MeshOp.getVertexNormal(v2);
					normal(n0);
					home.vertex(v0.xf(), v0.yf(), v0.zf(), getU(he0.getUVW()), getV(he0.getUVW()));
					normal(n1);
					home.vertex(v1.xf(), v1.yf(), v1.zf(), getU(he1.getUVW()), getV(he1.getUVW()));
					normal(n2);
					home.vertex(v2.xf(), v2.yf(), v2.zf(), getU(he2.getUVW()), getV(he2.getUVW()));
					home.endShape();
				}
			} else {
				for (int i = 0; i < tris.length; i += 3) {
					home.beginShape(PConstants.TRIANGLES);
					home.texture(texture);
					he0 = halfedges.get(tris[i]);
					v0 = he0.getVertex();
					// n0 = HE_MeshOp.getVertexNormal(v0);
					he1 = halfedges.get(tris[i + 1]);
					v1 = he1.getVertex();
					// n1 = HE_MeshOp.getVertexNormal(v1);
					he2 = halfedges.get(tris[i + 2]);
					v2 = he2.getVertex();
					// n2 = HE_MeshOp.getVertexNormal(v2);
					home.vertex(v0.xf(), v0.yf(), v0.zf(), getU(he0.getUVW()), getV(he0.getUVW()));
					home.vertex(v1.xf(), v1.yf(), v1.zf(), getU(he1.getUVW()), getV(he1.getUVW()));
					home.vertex(v2.xf(), v2.yf(), v2.zf(), getU(he2.getUVW()), getV(he2.getUVW()));
					home.endShape();
				}
			}
		}
	}

	private void drawFaceInt(final HE_Face f, final PImage texture, final boolean smooth) {
		if (!f.isVisible()) {
			return;
		}
		final int fo = f.getFaceDegree();
		final List<HE_Halfedge> halfedges = f.getFaceHalfedges();
		if (fo < 3 || halfedges.size() < 3) {
		} else if (fo == 3) {
			final int[] tri = new int[] { 0, 1, 2 };
			HE_Halfedge he0, he1, he2;
			HE_Vertex v0, v1, v2;
			WB_Coord n0, n1, n2;
			he0 = halfedges.get(tri[0]);
			v0 = he0.getVertex();
			n0 = HE_MeshOp.getVertexNormal(v0);
			he1 = halfedges.get(tri[1]);
			v1 = he1.getVertex();
			n1 = HE_MeshOp.getVertexNormal(v1);
			he2 = halfedges.get(tri[2]);
			v2 = he2.getVertex();
			n2 = HE_MeshOp.getVertexNormal(v2);
			if (smooth) {
				home.texture(texture);
				normal(n0);
				home.vertex(v0.xf(), v0.yf(), v0.zf(), getU(he0.getUVW()), getV(he0.getUVW()));
				normal(n1);
				home.vertex(v1.xf(), v1.yf(), v1.zf(), getU(he1.getUVW()), getV(he1.getUVW()));
				normal(n2);
				home.vertex(v2.xf(), v2.yf(), v2.zf(), getU(he2.getUVW()), getV(he2.getUVW()));
			} else {
				home.texture(texture);
				home.vertex(v0.xf(), v0.yf(), v0.zf(), getU(he0.getUVW()), getV(he0.getUVW()));
				home.vertex(v1.xf(), v1.yf(), v1.zf(), getU(he1.getUVW()), getV(he1.getUVW()));
				home.vertex(v2.xf(), v2.yf(), v2.zf(), getU(he2.getUVW()), getV(he2.getUVW()));
			}
		} else {
			final int[] tris = f.getTriangles();
			HE_Halfedge he0, he1, he2;
			HE_Vertex v0, v1, v2;
			WB_Coord n0, n1, n2;
			if (smooth) {
				for (int i = 0; i < tris.length; i += 3) {
					home.texture(texture);
					he0 = halfedges.get(tris[i]);
					v0 = he0.getVertex();
					n0 = HE_MeshOp.getVertexNormal(v0);
					he1 = halfedges.get(tris[i + 1]);
					v1 = he1.getVertex();
					n1 = HE_MeshOp.getVertexNormal(v1);
					he2 = halfedges.get(tris[i + 2]);
					v2 = he2.getVertex();
					n2 = HE_MeshOp.getVertexNormal(v2);
					normal(n0);
					home.vertex(v0.xf(), v0.yf(), v0.zf(), getU(he0.getUVW()), getV(he0.getUVW()));
					normal(n1);
					home.vertex(v1.xf(), v1.yf(), v1.zf(), getU(he1.getUVW()), getV(he1.getUVW()));
					normal(n2);
					home.vertex(v2.xf(), v2.yf(), v2.zf(), getU(he2.getUVW()), getV(he2.getUVW()));
				}
			} else {
				for (int i = 0; i < tris.length; i += 3) {
					home.texture(texture);
					he0 = halfedges.get(tris[i]);
					v0 = he0.getVertex();
					// n0 = HE_MeshOp.getVertexNormal(v0);
					he1 = halfedges.get(tris[i + 1]);
					v1 = he1.getVertex();
					// n1 = HE_MeshOp.getVertexNormal(v1);
					he2 = halfedges.get(tris[i + 2]);
					v2 = he2.getVertex();
					// n2 = HE_MeshOp.getVertexNormal(v2);
					home.vertex(v0.xf(), v0.yf(), v0.zf(), getU(he0.getUVW()), getV(he0.getUVW()));
					home.vertex(v1.xf(), v1.yf(), v1.zf(), getU(he1.getUVW()), getV(he1.getUVW()));
					home.vertex(v2.xf(), v2.yf(), v2.zf(), getU(he2.getUVW()), getV(he2.getUVW()));
				}
			}
		}
	}

	public void drawFace(final HE_Face f, final PImage[] textures) {
		drawFace(f, textures, false);
	}

	public void drawFace(final HE_Face f, final PImage[] textures, final boolean smooth) {
		if (!f.isVisible()) {
			return;
		}
		final int fo = f.getFaceDegree();
		final int fti = f.getTextureId();
		final List<HE_Halfedge> halfedges = f.getFaceHalfedges();
		if (fo < 3 || halfedges.size() < 3) {
		} else if (fo == 3) {
			final int[] tri = new int[] { 0, 1, 2 };
			HE_Halfedge he0, he1, he2;
			HE_Vertex v0, v1, v2;
			WB_Coord n0, n1, n2;
			he0 = halfedges.get(tri[0]);
			v0 = he0.getVertex();
			he1 = halfedges.get(tri[1]);
			v1 = he1.getVertex();
			he2 = halfedges.get(tri[2]);
			v2 = he2.getVertex();
			if (smooth) {
				n0 = HE_MeshOp.getVertexNormal(v0);
				n1 = HE_MeshOp.getVertexNormal(v1);
				n2 = HE_MeshOp.getVertexNormal(v2);
				home.beginShape(PConstants.TRIANGLES);
				if (fti >= 0 && fti < textures.length) {
					home.texture(textures[fti]);
					normal(n0);
					home.vertex(v0.xf(), v0.yf(), v0.zf(), getU(he0.getUVW()), getV(he0.getUVW()));
					normal(n1);
					home.vertex(v1.xf(), v1.yf(), v1.zf(), getU(he1.getUVW()), getV(he1.getUVW()));
					normal(n2);
					home.vertex(v2.xf(), v2.yf(), v2.zf(), getU(he2.getUVW()), getV(he2.getUVW()));
					home.endShape();
				} else {
					normal(n0);
					home.vertex(v0.xf(), v0.yf(), v0.zf());
					normal(n1);
					home.vertex(v1.xf(), v1.yf(), v1.zf());
					normal(n2);
					home.vertex(v2.xf(), v2.yf(), v2.zf());
					home.endShape();
				}
			} else {
				home.beginShape(PConstants.TRIANGLES);
				if (fti >= 0 && fti < textures.length) {
					home.texture(textures[fti]);
					home.vertex(v0.xf(), v0.yf(), v0.zf(), getU(he0.getUVW()), getV(he0.getUVW()));
					home.vertex(v1.xf(), v1.yf(), v1.zf(), getU(he1.getUVW()), getV(he1.getUVW()));
					home.vertex(v2.xf(), v2.yf(), v2.zf(), getU(he2.getUVW()), getV(he2.getUVW()));
					home.endShape();
				} else {
					home.vertex(v0.xf(), v0.yf(), v0.zf());
					home.vertex(v1.xf(), v1.yf(), v1.zf());
					home.vertex(v2.xf(), v2.yf(), v2.zf());
					home.endShape();
				}
			}
		} else {
			final int[] tris = f.getTriangles();
			HE_Halfedge he0, he1, he2;
			HE_Vertex v0, v1, v2;
			WB_Coord n0, n1, n2;
			if (smooth) {
				for (int i = 0; i < tris.length; i += 3) {
					he0 = halfedges.get(tris[i]);
					v0 = he0.getVertex();
					n0 = HE_MeshOp.getVertexNormal(v0);
					he1 = halfedges.get(tris[i + 1]);
					v1 = he1.getVertex();
					n1 = HE_MeshOp.getVertexNormal(v1);
					he2 = halfedges.get(tris[i + 2]);
					v2 = he2.getVertex();
					n2 = HE_MeshOp.getVertexNormal(v2);
					home.beginShape(PConstants.TRIANGLES);
					if (fti >= 0 && fti < textures.length) {
						home.texture(textures[fti]);
						normal(n0);
						home.vertex(v0.xf(), v0.yf(), v0.zf(), getU(he0.getUVW()), getV(he0.getUVW()));
						normal(n1);
						home.vertex(v1.xf(), v1.yf(), v1.zf(), getU(he1.getUVW()), getV(he1.getUVW()));
						normal(n2);
						home.vertex(v2.xf(), v2.yf(), v2.zf(), getU(he2.getUVW()), getV(he2.getUVW()));
						home.endShape();
					} else {
						normal(n0);
						home.vertex(v0.xf(), v0.yf(), v0.zf());
						normal(n1);
						home.vertex(v1.xf(), v1.yf(), v1.zf());
						normal(n2);
						home.vertex(v2.xf(), v2.yf(), v2.zf());
						home.endShape();
					}
				}
			} else {
				for (int i = 0; i < tris.length; i += 3) {
					home.beginShape(PConstants.TRIANGLES);
					he0 = halfedges.get(tris[i]);
					v0 = he0.getVertex();
					he1 = halfedges.get(tris[i + 1]);
					v1 = he1.getVertex();
					he2 = halfedges.get(tris[i + 2]);
					v2 = he2.getVertex();
					if (fti >= 0 && fti < textures.length) {
						home.texture(textures[fti]);
						home.vertex(v0.xf(), v0.yf(), v0.zf(), getU(he0.getUVW()), getV(he0.getUVW()));
						home.vertex(v1.xf(), v1.yf(), v1.zf(), getU(he1.getUVW()), getV(he1.getUVW()));
						home.vertex(v2.xf(), v2.yf(), v2.zf(), getU(he2.getUVW()), getV(he2.getUVW()));
						home.endShape();
					} else {
						home.vertex(v0.xf(), v0.yf(), v0.zf());
						home.vertex(v1.xf(), v1.yf(), v1.zf());
						home.vertex(v2.xf(), v2.yf(), v2.zf());
						home.endShape();
					}
				}
			}
		}
	}

	public void drawFace(final long key, final boolean smooth, final HE_Mesh mesh) {
		final HE_Face f = mesh.getFaceWithKey(key);
		if (f != null) {
			drawFace(f, smooth);
		}
	}

	public void drawFace(final long key, final HE_Mesh mesh) {
		final HE_Face f = mesh.getFaceWithKey(key);
		if (f != null) {
			drawFace(f, false);
		}
	}

	public void drawFaceFC(final HE_Face f) {
		drawFaceFC(f, false);
	}

	private void drawFaceFCInt(final HE_Face f) {
		drawFaceFCInt(f, false);
	}

	public void drawFaceFC(final HE_Face f, final boolean smooth) {
		if (!f.isVisible()) {
			return;
		}
		if (f.getFaceDegree() > 2) {
			home.pushStyle();
			home.fill(f.getColor());
			final int[] tris = f.getTriangles();
			final List<HE_Vertex> vertices = f.getFaceVertices();
			HE_Vertex v0, v1, v2;
			WB_Coord n0, n1, n2;
			if (smooth) {
				for (int i = 0; i < tris.length; i += 3) {
					home.beginShape(PConstants.TRIANGLES);
					v0 = vertices.get(tris[i]);
					n0 = HE_MeshOp.getVertexNormal(v0);
					v1 = vertices.get(tris[i + 1]);
					n1 = HE_MeshOp.getVertexNormal(v1);
					v2 = vertices.get(tris[i + 2]);
					n2 = HE_MeshOp.getVertexNormal(v2);
					normal(n0);
					vertex(v0);
					normal(n1);
					vertex(v1);
					normal(n2);
					vertex(v2);
					home.endShape();
				}
			} else {
				for (int i = 0; i < tris.length; i += 3) {
					home.beginShape(PConstants.TRIANGLES);
					v0 = vertices.get(tris[i]);
					v1 = vertices.get(tris[i + 1]);
					v2 = vertices.get(tris[i + 2]);
					vertex(v0);
					vertex(v1);
					vertex(v2);
					home.endShape();
				}
			}
			home.popStyle();
		}
	}

	private void drawFaceFCInt(final HE_Face f, final boolean smooth) {
		if (!f.isVisible()) {
			return;
		}
		if (f.getFaceDegree() > 2) {
			home.pushStyle();
			home.fill(f.getColor());
			final int[] tris = f.getTriangles();
			final List<HE_Vertex> vertices = f.getFaceVertices();
			HE_Vertex v0, v1, v2;
			WB_Coord n0, n1, n2;
			if (smooth) {
				for (int i = 0; i < tris.length; i += 3) {
					v0 = vertices.get(tris[i]);
					n0 = HE_MeshOp.getVertexNormal(v0);
					v1 = vertices.get(tris[i + 1]);
					n1 = HE_MeshOp.getVertexNormal(v1);
					v2 = vertices.get(tris[i + 2]);
					n2 = HE_MeshOp.getVertexNormal(v2);
					normal(n0);
					vertex(v0);
					normal(n1);
					vertex(v1);
					normal(n2);
					vertex(v2);
				}
			} else {
				for (int i = 0; i < tris.length; i += 3) {
					v0 = vertices.get(tris[i]);
					v1 = vertices.get(tris[i + 1]);
					v2 = vertices.get(tris[i + 2]);
					vertex(v0);
					vertex(v1);
					vertex(v2);
				}
			}
			home.popStyle();
		}
	}

	public void drawFaceHC(final HE_Face f) {
		drawFaceHC(f, false);
	}

	private void drawFaceHCInt(final HE_Face f) {
		drawFaceHCInt(f, false);
	}

	public void drawFaceHC(final HE_Face f, final boolean smooth) {
		if (!f.isVisible()) {
			return;
		}
		if (f.getFaceDegree() > 2) {
			final int[] tris = f.getTriangles();
			final List<HE_Vertex> vertices = f.getFaceVertices();
			final List<HE_Halfedge> halfedges = f.getFaceHalfedges();
			HE_Vertex v0, v1, v2;
			WB_Coord n0, n1, n2;
			if (smooth) {
				for (int i = 0; i < tris.length; i += 3) {
					home.beginShape(PConstants.TRIANGLES);
					v0 = vertices.get(tris[i]);
					n0 = HE_MeshOp.getVertexNormal(v0);
					v1 = vertices.get(tris[i + 1]);
					n1 = HE_MeshOp.getVertexNormal(v1);
					v2 = vertices.get(tris[i + 2]);
					n2 = HE_MeshOp.getVertexNormal(v2);
					home.fill(halfedges.get(tris[i]).getColor());
					normal(n0);
					vertex(v0);
					home.fill(halfedges.get(tris[i + 1]).getColor());
					normal(n1);
					vertex(v1);
					home.fill(halfedges.get(tris[i + 2]).getColor());
					normal(n2);
					vertex(v2);
					home.endShape();
				}
			} else {
				for (int i = 0; i < tris.length; i += 3) {
					home.beginShape(PConstants.TRIANGLES);
					v0 = vertices.get(tris[i]);
					v1 = vertices.get(tris[i + 1]);
					v2 = vertices.get(tris[i + 2]);
					home.fill(halfedges.get(tris[i]).getColor());
					vertex(v0);
					home.fill(halfedges.get(tris[i + 1]).getColor());
					vertex(v1);
					home.fill(halfedges.get(tris[i + 2]).getColor());
					vertex(v2);
					home.endShape();
				}
			}
		}
	}

	private void drawFaceHCInt(final HE_Face f, final boolean smooth) {
		if (!f.isVisible()) {
			return;
		}
		if (f.getFaceDegree() > 2) {
			final int[] tris = f.getTriangles();
			final List<HE_Vertex> vertices = f.getFaceVertices();
			final List<HE_Halfedge> halfedges = f.getFaceHalfedges();
			HE_Vertex v0, v1, v2;
			WB_Coord n0, n1, n2;
			if (smooth) {
				for (int i = 0; i < tris.length; i += 3) {
					v0 = vertices.get(tris[i]);
					n0 = HE_MeshOp.getVertexNormal(v0);
					v1 = vertices.get(tris[i + 1]);
					n1 = HE_MeshOp.getVertexNormal(v1);
					v2 = vertices.get(tris[i + 2]);
					n2 = HE_MeshOp.getVertexNormal(v2);
					home.fill(halfedges.get(tris[i]).getColor());
					normal(n0);
					vertex(v0);
					home.fill(halfedges.get(tris[i + 1]).getColor());
					normal(n1);
					vertex(v1);
					home.fill(halfedges.get(tris[i + 2]).getColor());
					normal(n2);
					vertex(v2);
				}
			} else {
				for (int i = 0; i < tris.length; i += 3) {
					v0 = vertices.get(tris[i]);
					v1 = vertices.get(tris[i + 1]);
					v2 = vertices.get(tris[i + 2]);
					home.fill(halfedges.get(tris[i]).getColor());
					vertex(v0);
					home.fill(halfedges.get(tris[i + 1]).getColor());
					vertex(v1);
					home.fill(halfedges.get(tris[i + 2]).getColor());
					vertex(v2);
				}
			}
		}
	}

	public void drawFaceNormal(final HE_Face f, final double d) {
		final WB_Coord p1 = HE_MeshOp.getFaceCenter(f);
		final WB_Point p2 = WB_Point.mul(HE_MeshOp.getFaceNormal(f), d).addSelf(p1);
		line(p1, p2);
	}

	public void drawFaceNormals(final HE_HalfedgeStructure mesh, final double d) {
		final HE_FaceIterator fItr = mesh.fItr();
		WB_Coord fc;
		WB_Coord fn;
		HE_Face f;
		while (fItr.hasNext()) {
			f = fItr.next();
			fc = HE_MeshOp.getFaceCenter(f);
			fn = HE_MeshOp.getFaceNormal(f);
			line(fc, WB_Point.addMul(fc, d, fn));
		}
	}

	public void drawFaceOffset(final HE_Face f, final double d) {
		if (!f.isVisible()) {
			return;
		}
		final int fo = f.getFaceDegree();
		final List<HE_Vertex> vertices = f.getFaceVertices();
		if (fo < 3 || vertices.size() < 3) {
		} else if (fo == 3) {
			final WB_Coord fn = HE_MeshOp.getFaceNormal(f);
			final int[] tri = new int[] { 0, 1, 2 };
			HE_Vertex v0, v1, v2;
			final float df = (float) d;
			home.beginShape(PConstants.TRIANGLES);
			v0 = vertices.get(tri[0]);
			v1 = vertices.get(tri[1]);
			v2 = vertices.get(tri[2]);
			home.vertex(v0.xf() + df * fn.xf(), v0.yf() + df * fn.yf(), v0.zf() + df * fn.zf());
			home.vertex(v1.xf() + df * fn.xf(), v1.yf() + df * fn.yf(), v1.zf() + df * fn.zf());
			home.vertex(v2.xf() + df * fn.xf(), v2.yf() + df * fn.yf(), v2.zf() + df * fn.zf());
			home.endShape();
		} else {
			final int[] tris = f.getTriangles();
			final WB_Coord fn = HE_MeshOp.getFaceNormal(f);
			HE_Vertex v0, v1, v2;
			final float df = (float) d;
			for (int i = 0; i < tris.length; i += 3) {
				home.beginShape(PConstants.TRIANGLES);
				v0 = vertices.get(tris[i]);
				v1 = vertices.get(tris[i + 1]);
				v2 = vertices.get(tris[i + 2]);
				home.vertex(v0.xf() + df * fn.xf(), v0.yf() + df * fn.yf(), v0.zf() + df * fn.zf());
				home.vertex(v1.xf() + df * fn.xf(), v1.yf() + df * fn.yf(), v1.zf() + df * fn.zf());
				home.vertex(v2.xf() + df * fn.xf(), v2.yf() + df * fn.yf(), v2.zf() + df * fn.zf());
				home.endShape();
			}
		}
	}

	public void drawFaceSmooth(final HE_Face f) {
		drawFace(f, true);
	}

	public void drawFaceSmooth(final long key, final HE_Mesh mesh) {
		final HE_Face f = mesh.getFaceWithKey(key);
		if (f != null) {
			drawFace(f, true);
		}
	}

	public void drawFaceSmoothFC(final HE_Face f) {
		drawFaceFC(f, true);
	}

	public void drawFaceSmoothFC(final long key, final HE_Mesh mesh) {
		final HE_Face f = mesh.getFaceWithKey(key);
		if (f != null) {
			drawFaceFC(f, true);
		}
	}

	public void drawFaceSmoothHC(final HE_Face f) {
		drawFaceHC(f, true);
	}

	public void drawFaceSmoothHC(final long key, final HE_Mesh mesh) {
		final HE_Face f = mesh.getFaceWithKey(key);
		if (f != null) {
			drawFaceHC(f, true);
		}
	}

	public void drawFaceSmoothVC(final HE_Face f) {
		drawFaceVC(f, true);
	}

	public void drawFaceSmoothVC(final long key, final HE_Mesh mesh) {
		final HE_Face f = mesh.getFaceWithKey(key);
		if (f != null) {
			drawFaceVC(f, true);
		}
	}

	public void drawFaceVC(final HE_Face f) {
		drawFaceVC(f, false);
	}

	private void drawFaceVCInt(final HE_Face f) {
		drawFaceVCInt(f, false);
	}

	public void drawFaceVC(final HE_Face f, final boolean smooth) {
		if (!f.isVisible()) {
			return;
		}
		if (f.getFaceDegree() > 2) {
			final int[] tris = f.getTriangles();
			final List<HE_Vertex> vertices = f.getFaceVertices();
			HE_Vertex v0, v1, v2;
			WB_Coord n0, n1, n2;
			if (smooth) {
				for (int i = 0; i < tris.length; i += 3) {
					home.beginShape(PConstants.TRIANGLES);
					v0 = vertices.get(tris[i]);
					n0 = HE_MeshOp.getVertexNormal(v0);
					v1 = vertices.get(tris[i + 1]);
					n1 = HE_MeshOp.getVertexNormal(v1);
					v2 = vertices.get(tris[i + 2]);
					n2 = HE_MeshOp.getVertexNormal(v2);
					home.fill(v0.getColor());
					normal(n0);
					vertex(v0);
					home.fill(v1.getColor());
					normal(n1);
					vertex(v1);
					home.fill(v2.getColor());
					normal(n2);
					vertex(v2);
					home.endShape();
				}
			} else {
				for (int i = 0; i < tris.length; i += 3) {
					home.beginShape(PConstants.TRIANGLES);
					v0 = vertices.get(tris[i]);
					v1 = vertices.get(tris[i + 1]);
					v2 = vertices.get(tris[i + 2]);
					home.fill(v0.getColor());
					vertex(v0);
					home.fill(v1.getColor());
					vertex(v1);
					home.fill(v2.getColor());
					vertex(v2);
					home.endShape();
				}
			}
		}
	}

	private void drawFaceVCInt(final HE_Face f, final boolean smooth) {
		if (!f.isVisible()) {
			return;
		}
		if (f.getFaceDegree() > 2) {
			final int[] tris = f.getTriangles();
			final List<HE_Vertex> vertices = f.getFaceVertices();
			HE_Vertex v0, v1, v2;
			WB_Coord n0, n1, n2;
			if (smooth) {
				for (int i = 0; i < tris.length; i += 3) {
					v0 = vertices.get(tris[i]);
					n0 = HE_MeshOp.getVertexNormal(v0);
					v1 = vertices.get(tris[i + 1]);
					n1 = HE_MeshOp.getVertexNormal(v1);
					v2 = vertices.get(tris[i + 2]);
					n2 = HE_MeshOp.getVertexNormal(v2);
					home.fill(v0.getColor());
					normal(n0);
					vertex(v0);
					home.fill(v1.getColor());
					normal(n1);
					vertex(v1);
					home.fill(v2.getColor());
					normal(n2);
					vertex(v2);
				}
			} else {
				for (int i = 0; i < tris.length; i += 3) {
					v0 = vertices.get(tris[i]);
					v1 = vertices.get(tris[i + 1]);
					v2 = vertices.get(tris[i + 2]);
					home.fill(v0.getColor());
					vertex(v0);
					home.fill(v1.getColor());
					vertex(v1);
					home.fill(v2.getColor());
					vertex(v2);
				}
			}
		}
	}

	public void drawFaces(final Collection<? extends HE_HalfedgeStructure> meshes) {
		final Iterator<? extends HE_HalfedgeStructure> mItr = meshes.iterator();
		while (mItr.hasNext()) {
			drawFaces(mItr.next());
		}
	}

	public void drawFaces(final HE_MeshCollection meshes) {
		final HE_MeshIterator mItr = meshes.mItr();
		while (mItr.hasNext()) {
			drawFaces(mItr.next());
		}
	}

	public void drawFaces(final HE_HalfedgeStructure mesh) {
		if (mesh == null) {
			return;
		}
		final HE_FaceIterator fItr = mesh.fItr();
		home.beginShape(PConstants.TRIANGLES);
		while (fItr.hasNext()) {
			drawFaceInt(fItr.next());
		}
		home.endShape();
	}

	public void drawFaces(final HE_HalfedgeStructure mesh, final PImage texture) {
		if (mesh == null) {
			return;
		}
		final HE_FaceIterator fItr = mesh.fItr();
		home.beginShape(PConstants.TRIANGLES);
		while (fItr.hasNext()) {
			drawFaceInt(fItr.next(), texture);
		}
		home.endShape();
	}

	public void drawFaces(final HE_HalfedgeStructure mesh, final PImage[] textures) {
		if (mesh == null) {
			return;
		}
		final HE_FaceIterator fItr = mesh.fItr();
		while (fItr.hasNext()) {
			drawFace(fItr.next(), textures);
		}
	}

	public void drawFacesFC(final HE_HalfedgeStructure mesh) {
		if (mesh == null) {
			return;
		}
		final HE_FaceIterator fItr = mesh.fItr();
		home.beginShape(PConstants.TRIANGLES);
		while (fItr.hasNext()) {
			try {
				drawFaceFCInt(fItr.next());
			} catch (final Exception e) {
			}
		}
		home.endShape();
	}

	public void drawFacesHC(final HE_HalfedgeStructure mesh) {
		if (mesh == null) {
			return;
		}
		final HE_FaceIterator fItr = mesh.fItr();
		home.beginShape(PConstants.TRIANGLES);
		while (fItr.hasNext()) {
			drawFaceHCInt(fItr.next());
		}
		home.endShape();
	}

	public void drawFacesSmooth(final HE_HalfedgeStructure mesh) {
		if (mesh == null) {
			return;
		}
		final HE_FaceIterator fItr = mesh.fItr();
		home.beginShape(PConstants.TRIANGLES);
		while (fItr.hasNext()) {
			drawFaceInt(fItr.next(), true);
		}
		home.endShape();
	}

	public void drawFacesSmooth(final HE_HalfedgeStructure mesh, final PImage texture) {
		if (mesh == null) {
			return;
		}
		final HE_FaceIterator fItr = mesh.fItr();
		home.beginShape(PConstants.TRIANGLES);
		while (fItr.hasNext()) {
			drawFaceInt(fItr.next(), texture, true);
		}
		home.endShape();
	}

	public void drawFacesSmooth(final HE_HalfedgeStructure mesh, final PImage[] textures) {
		if (mesh == null) {
			return;
		}
		final HE_FaceIterator fItr = mesh.fItr();
		while (fItr.hasNext()) {
			drawFace(fItr.next(), textures, true);
		}
	}

	public void drawFacesSmoothFC(final HE_HalfedgeStructure mesh) {
		if (mesh == null) {
			return;
		}
		final HE_FaceIterator fItr = mesh.fItr();
		home.beginShape(PConstants.TRIANGLES);
		while (fItr.hasNext()) {
			drawFaceFCInt(fItr.next(), true);
		}
		home.endShape();
	}

	public void drawFacesSmoothHC(final HE_HalfedgeStructure mesh) {
		if (mesh == null) {
			return;
		}
		final HE_FaceIterator fItr = mesh.fItr();
		home.beginShape(PConstants.TRIANGLES);
		while (fItr.hasNext()) {
			drawFaceHCInt(fItr.next(), true);
		}
		home.endShape();
	}

	public void drawFacesSmoothVC(final HE_HalfedgeStructure mesh) {
		if (mesh == null) {
			return;
		}
		final HE_FaceIterator fItr = mesh.fItr();
		home.beginShape(PConstants.TRIANGLES);
		while (fItr.hasNext()) {
			drawFaceVCInt(fItr.next(), true);
		}
		home.endShape();
	}

	public void drawFacesVC(final HE_HalfedgeStructure mesh) {
		if (mesh == null) {
			return;
		}
		final HE_FaceIterator fItr = mesh.fItr();
		home.beginShape(PConstants.TRIANGLES);
		while (fItr.hasNext()) {
			drawFaceVCInt(fItr.next());
		}
		home.endShape();
	}

	public void drawFacesWithInternalLabel(final int label, final HE_HalfedgeStructure mesh) {
		if (mesh == null) {
			return;
		}
		final HE_FaceIterator fItr = mesh.fItr();
		home.beginShape(PConstants.TRIANGLES);
		HE_Face f;
		while (fItr.hasNext()) {
			f = fItr.next();
			if (f.getInternalLabel() == label) {
				drawFaceInt(f);
			}
		}
		home.endShape();
	}

	public void drawFacesWithLabel(final int label, final HE_HalfedgeStructure mesh) {
		if (mesh == null) {
			return;
		}
		final HE_FaceIterator fItr = mesh.fItr();
		home.beginShape(PConstants.TRIANGLES);
		HE_Face f;
		while (fItr.hasNext()) {
			f = fItr.next();
			if (f.getLabel() == label) {
				drawFaceInt(f);
			}
		}
		home.endShape();
	}

	public void drawNetwork(final WB_Network network) {
		final List<Connection> connections = network.getConnections();
		for (int i = 0; i < network.getNumberOfConnections(); i++) {
			drawNetworkConnection(connections.get(i));
		}
	}

	public void drawNetworkNode(final Node node, final double s) {
		home.pushMatrix();
		translate(node);
		home.box((float) s);
		home.popMatrix();
	}

	public void drawNetworkConnection(final Connection connection) {
		line(connection.start(), connection.end());
	}

	public void drawHalfedgesWithInternalLabel(final int label, final HE_HalfedgeStructure mesh) {
		final Iterator<HE_Halfedge> heItr = mesh.heItr();
		HE_Halfedge he;
		while (heItr.hasNext()) {
			he = heItr.next();
			if (he.isVisible()) {
				if (he.getInternalLabel() == label) {
					line(he.getVertex(), he.getEndVertex());
				}
			}
		}
	}

	public void drawHalfedgesWithLabel(final int label, final HE_HalfedgeStructure mesh) {
		final Iterator<HE_Halfedge> heItr = mesh.heItr();
		HE_Halfedge he;
		while (heItr.hasNext()) {
			he = heItr.next();
			if (he.isVisible()) {
				if (he.getLabel() == label) {
					line(he.getVertex(), he.getEndVertex().getPosition().mulAddMul(0.8, 0.2, he.getVertex()));
				}
			}
		}
	}

	public void drawHalfedge(final HE_Halfedge he, final double d, final double s) {
		if (!he.isVisible()) {
			return;
		}
		final WB_Point c = new WB_Point(HE_MeshOp.getHalfedgeCenter(he));
		c.addMulSelf(d, HE_MeshOp.getHalfedgeNormal(he));
		home.stroke(255, 0, 0);
		line(he.getVertex(), c);
		if (HE_MeshOp.getHalfedgeType(he) == WB_Classification.CONVEX) {
			home.stroke(0, 255, 0);
		} else if (HE_MeshOp.getHalfedgeType(he) == WB_Classification.CONCAVE) {
			home.stroke(255, 0, 0);
		} else {
			home.stroke(0, 0, 255);
		}
		home.pushMatrix();
		translate(c);
		home.box((float) s);
		home.popMatrix();
	}

	public void drawHalfedge(final HE_Halfedge he, final double d, final double s, final double f) {
		if (!he.isVisible()) {
			return;
		}
		final WB_Point c = geometryfactory.createInterpolatedPoint(he.getVertex(), he.getEndVertex(), f);
		c.addMulSelf(d, HE_MeshOp.getHalfedgeNormal(he));
		line(he.getVertex(), c);
		if (HE_MeshOp.getHalfedgeType(he) == WB_Classification.CONVEX) {
			home.stroke(0, 255, 0);
		} else if (HE_MeshOp.getHalfedgeType(he) == WB_Classification.CONCAVE) {
			home.stroke(255, 0, 0);
		} else {
			home.stroke(0, 0, 255);
		}
		home.pushMatrix();
		translate(c);
		home.box((float) s);
		home.popMatrix();
	}

	public void drawHalfedge(final long key, final double d, final double s, final HE_Mesh mesh) {
		final HE_Halfedge he = mesh.getHalfedgeWithKey(key);
		drawHalfedge(he, d, s);
	}

	public void drawHalfedges(final HE_HalfedgeStructure mesh, final double d) {
		WB_Point c;
		HE_Halfedge he;
		final Iterator<HE_Halfedge> heItr = mesh.heItr();
		home.pushStyle();
		while (heItr.hasNext()) {
			he = heItr.next();
			if (he.isVisible()) {
				if (he.getFace() != null) {
					c = new WB_Point(HE_MeshOp.getHalfedgeCenter(he));
					c.addMulSelf(d, HE_MeshOp.getHalfedgeNormal(he));
					home.stroke(255, 0, 0);
					line(he.getVertex(), c);
					if (HE_MeshOp.getHalfedgeType(he) == WB_Classification.CONVEX) {
						home.stroke(0, 255, 0);
						home.fill(0, 255, 0);
					} else if (HE_MeshOp.getHalfedgeType(he) == WB_Classification.CONCAVE) {
						home.stroke(255, 0, 0);
						home.fill(255, 0, 0);
					} else {
						home.stroke(0, 0, 255);
						home.fill(0, 0, 255);
					}
					home.pushMatrix();
					translate(c);
					home.box((float) d);
					home.popMatrix();
				} else {
					c = new WB_Point(HE_MeshOp.getHalfedgeCenter(he));
					c.addMulSelf(-d, HE_MeshOp.getHalfedgeNormal(he.getPair()));
					home.stroke(255, 0, 0);
					line(he.getVertex(), c);
					home.stroke(0, 255, 255);
					home.pushMatrix();
					translate(c);
					home.box((float) d);
					home.popMatrix();
				}
			}
		}
		home.popStyle();
	}

	public void drawHalfedges(final HE_HalfedgeStructure mesh, final double d, final double f) {
		WB_Point c;
		HE_Halfedge he;
		final Iterator<HE_Halfedge> heItr = mesh.heItr();
		home.pushStyle();
		while (heItr.hasNext()) {
			he = heItr.next();
			if (he.isVisible()) {
				if (he.getFace() != null) {
					c = geometryfactory.createInterpolatedPoint(he.getVertex(), he.getEndVertex(), f);
					c.addMulSelf(d, HE_MeshOp.getHalfedgeNormal(he));
					home.stroke(255, 0, 0);
					line(he.getVertex(), c);
					if (HE_MeshOp.getHalfedgeType(he) == WB_Classification.CONVEX) {
						home.stroke(0, 255, 0);
						home.fill(0, 255, 0);
					} else if (HE_MeshOp.getHalfedgeType(he) == WB_Classification.CONCAVE) {
						home.stroke(255, 0, 0);
						home.fill(255, 0, 0);
					} else {
						home.stroke(0, 0, 255);
						home.fill(0, 0, 255);
					}
					home.pushMatrix();
					translate(c);
					home.box((float) d);
					home.popMatrix();
				} else {
					c = geometryfactory.createInterpolatedPoint(he.getVertex(), he.getEndVertex(), f);
					c.addMulSelf(-d, HE_MeshOp.getHalfedgeNormal(he.getPair()));
					home.stroke(255, 0, 0);
					line(he.getVertex(), c);
					home.stroke(0, 255, 255);
					home.pushMatrix();
					translate(c);
					home.box((float) d);
					home.popMatrix();
				}
			}
		}
		home.popStyle();
	}

	public void drawHalfedgeSimple(final HE_Halfedge he, final double d, final double s) {
		if (!he.isVisible()) {
			return;
		}
		final WB_Point c = new WB_Point(HE_MeshOp.getHalfedgeCenter(he));
		c.addMulSelf(d, HE_MeshOp.getHalfedgeNormal(he));
		line(he.getVertex(), c);
		home.pushMatrix();
		translate(c);
		home.box((float) s);
		home.popMatrix();
	}

	public void drawLine(final WB_Line L, final double d) {
		home.line((float) (L.getOrigin().xd() - d * L.getDirection().xd()),
				(float) (L.getOrigin().yd() - d * L.getDirection().yd()),
				(float) (L.getOrigin().zd() - d * L.getDirection().zd()),
				(float) (L.getOrigin().xd() + d * L.getDirection().xd()),
				(float) (L.getOrigin().yd() + d * L.getDirection().yd()),
				(float) (L.getOrigin().zd() + d * L.getDirection().zd()));
	}

	public void drawLineEmbedded2D(final WB_Line L, final double d, final WB_Map2D map) {
		drawSegmentEmbedded2D(WB_Point.addMul(L.getOrigin(), -d, L.getDirection()),
				WB_Point.addMul(L.getOrigin(), d, L.getDirection()), map);
	}

	public void drawLineMapped(final WB_Line L, final double d, final WB_Map map) {
		drawSegmentMapped(WB_Point.addMul(L.getOrigin(), -d, L.getDirection()),
				WB_Point.addMul(L.getOrigin(), d, L.getDirection()), map);
	}

	public void drawLineUnmapped(final WB_Line L, final double d, final WB_Map map) {
		drawSegmentUnmapped(WB_Point.addMul(L.getOrigin(), -d, L.getDirection()),
				WB_Point.addMul(L.getOrigin(), d, L.getDirection()), map);
	}

	public void drawMesh(final WB_SimpleMesh mesh) {
		if (mesh == null) {
			return;
		}
		if (mesh.getNumberOfVertices() == 0) {
			return;
		}
		for (final int[] face : mesh.getFacesAsInt()) {
			drawPolygon(face, mesh.getPoints());
		}
	}

	public void drawMeshEdges(final WB_SimpleMesh mesh) {
		if (mesh == null) {
			return;
		}
		if (mesh.getNumberOfVertices() == 0) {
			return;
		}
		for (final int[] face : mesh.getFacesAsInt()) {
			drawPolygonEdges(face, mesh.getPoints());
		}
	}

	public void drawMeshFaces(final WB_SimpleMesh mesh) {
		if (mesh == null) {
			return;
		}
		if (mesh.getNumberOfVertices() == 0) {
			return;
		}
		for (final int[] face : mesh.getFacesAsInt()) {
			drawPolygon(face, mesh.getPoints());
		}
	}

	public void drawNetworkNodes(final WB_Network network, final double s) {
		final List<Node> nodes = network.getNodes();
		for (int i = 0; i < network.getNumberOfNodes(); i++) {
			drawNetworkNode(nodes.get(i), s);
		}
	}

	public void drawPath(final HE_Path path) {
		if (path.getPathOrder() == 1) {
			drawSegment(path.getPathHalfedge().getStartVertex(), path.getPathHalfedge().getEndVertex());
		}
		home.beginShape();
		for (final HE_Vertex v : path.getPathVertices()) {
			home.vertex(v.xf(), v.yf(), v.zf());
		}
		if (path.isLoop()) {
			home.endShape(PConstants.CLOSE);
		} else {
			home.endShape();
		}
	}

	public void drawPlane(final WB_Plane P, final double d) {
		home.beginShape(PConstants.QUAD);
		home.vertex((float) (P.getOrigin().xd() - d * P.getU().xd() - d * P.getV().xd()),
				(float) (P.getOrigin().yd() - d * P.getU().yd() - d * P.getV().yd()),
				(float) (P.getOrigin().zd() - d * P.getU().zd() - d * P.getV().zd()));
		home.vertex((float) (P.getOrigin().xd() - d * P.getU().xd() + d * P.getV().xd()),
				(float) (P.getOrigin().yd() - d * P.getU().yd() + d * P.getV().yd()),
				(float) (P.getOrigin().zd() - d * P.getU().zd() + d * P.getV().zd()));
		home.vertex((float) (P.getOrigin().xd() + d * P.getU().xd() + d * P.getV().xd()),
				(float) (P.getOrigin().yd() + d * P.getU().yd() + d * P.getV().yd()),
				(float) (P.getOrigin().zd() + d * P.getU().zd() + d * P.getV().zd()));
		home.vertex((float) (P.getOrigin().xd() + d * P.getU().xd() - d * P.getV().xd()),
				(float) (P.getOrigin().yd() + d * P.getU().yd() - d * P.getV().yd()),
				(float) (P.getOrigin().zd() + d * P.getU().zd() - d * P.getV().zd()));
		home.endShape();
	}

	public void drawPlaneHatch(final WB_Plane P, final double d, final int h) {
		home.beginShape(PConstants.QUAD);
		home.vertex((float) (P.getOrigin().xd() - d * P.getU().xd() - d * P.getV().xd()),
				(float) (P.getOrigin().yd() - d * P.getU().yd() - d * P.getV().yd()),
				(float) (P.getOrigin().zd() - d * P.getU().zd() - d * P.getV().zd()));
		home.vertex((float) (P.getOrigin().xd() - d * P.getU().xd() + d * P.getV().xd()),
				(float) (P.getOrigin().yd() - d * P.getU().yd() + d * P.getV().yd()),
				(float) (P.getOrigin().zd() - d * P.getU().zd() + d * P.getV().zd()));
		home.vertex((float) (P.getOrigin().xd() + d * P.getU().xd() + d * P.getV().xd()),
				(float) (P.getOrigin().yd() + d * P.getU().yd() + d * P.getV().yd()),
				(float) (P.getOrigin().zd() + d * P.getU().zd() + d * P.getV().zd()));
		home.vertex((float) (P.getOrigin().xd() + d * P.getU().xd() - d * P.getV().xd()),
				(float) (P.getOrigin().yd() + d * P.getU().yd() - d * P.getV().yd()),
				(float) (P.getOrigin().zd() + d * P.getU().zd() - d * P.getV().zd()));
		home.endShape();
		if (h > 0) {
			home.beginShape(PConstants.LINES);
			home.vertex((float) (P.getOrigin().xd() - d * P.getU().xd() + d * P.getV().xd()),
					(float) (P.getOrigin().yd() - d * P.getU().yd() + d * P.getV().yd()),
					(float) (P.getOrigin().zd() - d * P.getU().zd() + d * P.getV().zd()));
			home.vertex((float) (P.getOrigin().xd() + d * P.getU().xd() - d * P.getV().xd()),
					(float) (P.getOrigin().yd() + d * P.getU().yd() - d * P.getV().yd()),
					(float) (P.getOrigin().zd() + d * P.getU().zd() - d * P.getV().zd()));
			home.endShape();
		}
		if (h > 1) {
			home.beginShape(PConstants.LINES);
			final double dh = 2 * d / h;
			for (int i = 1; i < h; i++) {
				home.vertex((float) (P.getOrigin().xd() - d * P.getU().xd() + (-d + i * dh) * P.getV().xd()),
						(float) (P.getOrigin().yd() - d * P.getU().yd() + (-d + i * dh) * P.getV().yd()),
						(float) (P.getOrigin().zd() - d * P.getU().zd() + (-d + i * dh) * P.getV().zd()));
				home.vertex((float) (P.getOrigin().xd() + (-d + i * dh) * P.getU().xd() - d * P.getV().xd()),
						(float) (P.getOrigin().yd() + (-d + i * dh) * P.getU().yd() - d * P.getV().yd()),
						(float) (P.getOrigin().zd() + (-d + i * dh) * P.getU().zd() - d * P.getV().zd()));
				home.vertex((float) (P.getOrigin().xd() + d * P.getU().xd() - (-d + i * dh) * P.getV().xd()),
						(float) (P.getOrigin().yd() + d * P.getU().yd() - (-d + i * dh) * P.getV().yd()),
						(float) (P.getOrigin().zd() + d * P.getU().zd() - (-d + i * dh) * P.getV().zd()));
				home.vertex((float) (P.getOrigin().xd() - (-d + i * dh) * P.getU().xd() + d * P.getV().xd()),
						(float) (P.getOrigin().yd() - (-d + i * dh) * P.getU().yd() + d * P.getV().yd()),
						(float) (P.getOrigin().zd() - (-d + i * dh) * P.getU().zd() + d * P.getV().zd()));
			}
			home.endShape();
		}
	}

	public void drawPoint(final Collection<? extends WB_Coord> points, final double d) {
		for (final WB_Coord v : points) {
			drawPoint(v, d);
		}
	}

	public void drawPoint(final WB_Coord p) {
		home.point(p.xf(), p.yf(), p.zf());
	}

	public void drawPoint(final WB_Coord[] points) {
		for (final WB_Coord v : points) {
			drawPoint(v);
		}
	}

	public void drawPoint(final Collection<? extends WB_Coord> points) {
		for (final WB_Coord v : points) {
			drawPoint(v);
		}
	}

	public void drawPoint(final WB_CoordCollection points) {
		for (int i = 0; i < points.size(); i++) {
			drawPoint(points.get(i));
		}
	}

	public void drawPoint(final WB_CoordCollection points, final double r) {
		for (int i = 0; i < points.size(); i++) {
			drawPoint(points.get(i), r);
		}
	}

	public void drawPoint(final WB_Coord p, final double r) {
		home.pushMatrix();
		translate(p);
		home.box((float) r);
		home.popMatrix();
	}

	public void drawPoint(final WB_Coord[] points, final double d) {
		for (final WB_Coord v : points) {
			home.pushMatrix();
			translate(v);
			home.box((float) d);
			home.popMatrix();
		}
	}

	public void drawPointEmbedded2D(final Collection<? extends WB_Coord> points, final double r, final WB_Map2D map) {
		for (final WB_Coord p : points) {
			drawPointEmbedded2D(p, r, map);
		}
	}

	public void drawPointEmbedded2D(final Collection<? extends WB_Coord> points, final WB_Map2D map) {
		for (final WB_Coord p : points) {
			drawPointEmbedded2D(p, map);
		}
	}

	public void drawPointEmbedded2D(final WB_Coord p, final double r, final WB_Map2D map) {
		final WB_Point q = new WB_Point();
		map.mapPoint3D(p, q);
		map.unmapPoint2D(q.xd(), q.yd(), q);
		drawPoint(q, r);
	}

	public void drawPointEmbedded2D(final WB_Coord p, final WB_Map2D map) {
		final WB_Point q = new WB_Point();
		map.mapPoint3D(p, q);
		map.unmapPoint2D(q.xd(), q.yd(), q);
		drawPoint(q);
	}

	public void drawPointEmbedded2D(final WB_Coord[] points, final double r, final WB_Map2D map) {
		for (final WB_Coord p : points) {
			drawPointEmbedded2D(p, r, map);
		}
	}

	public void drawPointEmbedded2D(final WB_Coord[] points, final WB_Map2D map) {
		for (final WB_Coord p : points) {
			drawPointEmbedded2D(p, map);
		}
	}

	public void drawPointMapped(final Collection<? extends WB_Coord> points, final double r, final WB_Map map) {
		for (final WB_Coord p : points) {
			drawPointMapped(p, r, map);
		}
	}

	public void drawPointMapped(final Collection<? extends WB_Coord> points, final WB_Map map) {
		for (final WB_Coord p : points) {
			drawPointMapped(p, map);
		}
	}

	public void drawPointMapped(final WB_Coord p, final double r, final WB_Map map) {
		final WB_Point q = new WB_Point();
		map.mapPoint3D(p, q);
		drawPoint(q, r);
	}

	public void drawPointMapped(final WB_Coord p, final WB_Map map) {
		final WB_Point q = new WB_Point();
		map.mapPoint3D(p, q);
		drawPoint(q);
	}

	public void drawPointMapped(final WB_Coord[] points, final double r, final WB_Map map) {
		for (final WB_Coord p : points) {
			drawPointMapped(p, r, map);
		}
	}

	public void drawPointMapped(final WB_Coord[] points, final WB_Map map) {
		for (final WB_Coord p : points) {
			drawPointMapped(p, map);
		}
	}

	public void drawPointUnmapped(final Collection<? extends WB_Coord> points, final double r, final WB_Map map) {
		for (final WB_Coord p : points) {
			drawPointUnmapped(p, r, map);
		}
	}

	public void drawPointUnmapped(final Collection<? extends WB_Coord> points, final WB_Map map) {
		for (final WB_Coord p : points) {
			drawPointUnmapped(p, map);
		}
	}

	public void drawPointUnmapped(final WB_Coord p, final double r, final WB_Map map) {
		final WB_Point q = new WB_Point();
		map.mapPoint3D(p, q);
		drawPoint(q, r);
	}

	public void drawPointUnmapped(final WB_Coord p, final WB_Map map) {
		final WB_Point q = new WB_Point();
		map.mapPoint3D(p, q);
		drawPoint(q);
	}

	public void drawPointUnmapped(final WB_Coord[] points, final double r, final WB_Map map) {
		for (final WB_Coord p : points) {
			drawPointUnmapped(p, r, map);
		}
	}

	public void drawPointUnmapped(final WB_Coord[] points, final WB_Map map) {
		for (final WB_Coord p : points) {
			drawPointUnmapped(p, map);
		}
	}

	public void drawPolygon(final Collection<? extends WB_Polygon> polygons) {
		final Iterator<? extends WB_Polygon> polyItr = polygons.iterator();
		while (polyItr.hasNext()) {
			drawPolygon(polyItr.next());
		}
	}

	public void drawPolygon(final int[] indices, final List<? extends WB_Coord> points) {
		if (points != null && indices != null) {
			home.beginShape(PConstants.POLYGON);
			for (final int indice : indices) {
				home.vertex(points.get(indice).xf(), points.get(indice).yf(), points.get(indice).zf());
			}
			home.endShape(PConstants.CLOSE);
		}
	}

	public void drawPolygon(final int[] indices, final WB_CoordCollection points) {
		if (points != null && indices != null) {
			home.beginShape(PConstants.POLYGON);
			for (final int indice : indices) {
				home.vertex(points.get(indice).xf(), points.get(indice).yf(), points.get(indice).zf());
			}
			home.endShape(PConstants.CLOSE);
		}
	}

	public void drawPolygon(final WB_Polygon P) {
		final int[] tris = P.getTriangles();
		for (int i = 0; i < tris.length; i += 3) {
			drawTriangle(P.getPoint(tris[i]), P.getPoint(tris[i + 1]), P.getPoint(tris[i + 2]));
		}
	}

	public void drawPolygonEdges(final Collection<? extends WB_Polygon> polygons) {
		final Iterator<? extends WB_Polygon> polyItr = polygons.iterator();
		while (polyItr.hasNext()) {
			drawPolygonEdges(polyItr.next());
		}
	}

	public void drawPolygonEdges(final int[] indices, final List<? extends WB_Coord> points) {
		if (points != null && indices != null) {
			home.beginShape();
			for (final int indice : indices) {
				home.vertex(points.get(indice).xf(), points.get(indice).yf(), points.get(indice).zf());
			}
			home.endShape(PConstants.CLOSE);
		}
	}

	public void drawPolygonEdges(final int[] indices, final WB_CoordCollection points) {
		if (points != null && indices != null) {
			home.beginShape();
			for (final int indice : indices) {
				home.vertex(points.get(indice).xf(), points.get(indice).yf(), points.get(indice).zf());
			}
			home.endShape(PConstants.CLOSE);
		}
	}

	public void drawPolygonEdges(final WB_Polygon P) {
		final int[] npc = P.getNumberOfPointsPerContour();
		int index = 0;
		for (int i = 0; i < P.getNumberOfContours(); i++) {
			home.beginShape();
			for (int j = 0; j < npc[i]; j++) {
				vertex(P.getPoint(index++));
			}
			home.endShape(PConstants.CLOSE);
		}
	}

	public void drawPolygonEdgesEmbedded2D(final WB_Polygon P, final WB_Map2D map) {
		final int[] npc = P.getNumberOfPointsPerContour();
		int index = 0;
		for (int i = 0; i < P.getNumberOfContours(); i++) {
			home.beginShape();
			for (int j = 0; j < npc[i]; j++) {
				vertexEmbedded2D(P.getPoint(index++), map);
			}
			home.endShape(PConstants.CLOSE);
		}
	}

	public void drawPolygonEdgesMapped(final WB_Polygon P, final WB_Map map) {
		final int[] npc = P.getNumberOfPointsPerContour();
		int index = 0;
		for (int i = 0; i < P.getNumberOfContours(); i++) {
			home.beginShape();
			for (int j = 0; j < npc[i]; j++) {
				vertexMapped(P.getPoint(index++), map);
			}
			home.endShape(PConstants.CLOSE);
		}
	}

	public void drawPolygonEdgesUnmapped(final WB_Polygon P, final WB_Map map) {
		final int[] npc = P.getNumberOfPointsPerContour();
		int index = 0;
		for (int i = 0; i < P.getNumberOfContours(); i++) {
			home.beginShape();
			for (int j = 0; j < npc[i]; j++) {
				vertexUnmapped(P.getPoint(index++), map);
			}
			home.endShape(PConstants.CLOSE);
		}
	}

	public void drawPolygonEmbedded2D(final WB_Polygon P, final WB_Map2D map) {
		final int[] tris = P.getTriangles();
		for (int i = 0; i < tris.length; i += 3) {
			drawTriangleEmbedded2D(P.getPoint(tris[i]), P.getPoint(tris[i + 1]), P.getPoint(tris[i + 2]), map);
		}
	}

	public void drawPolygonMapped(final WB_Polygon P, final WB_Map map) {
		final int[] tris = P.getTriangles();
		for (int i = 0; i < tris.length; i += 3) {
			drawTriangleMapped(P.getPoint(tris[i]), P.getPoint(tris[i + 1]), P.getPoint(tris[i + 2]), map);
		}
	}

	public void drawPolygonUnmapped(final WB_Polygon P, final WB_Map map) {
		final int[] tris = P.getTriangles();
		for (int i = 0; i < tris.length; i += 3) {
			drawTriangleUnmapped(P.getPoint(tris[i]), P.getPoint(tris[i + 1]), P.getPoint(tris[i + 2]), map);
		}
	}

	public void drawPolygonVertices(final Collection<WB_Polygon> polygons, final double d) {
		final Iterator<WB_Polygon> polyItr = polygons.iterator();
		while (polyItr.hasNext()) {
			drawPolygonVertices(polyItr.next(), d);
		}
	}

	public void drawPolygonVertices(final WB_Polygon polygon, final double d) {
		WB_Coord v1;
		final int n = polygon.getNumberOfPoints();
		for (int i = 0; i < n; i++) {
			v1 = polygon.getPoint(i);
			home.pushMatrix();
			translate(v1);
			home.box((float) d);
			home.popMatrix();
		}
	}

	public void drawPolyLine(final WB_PolyLine P) {
		for (int i = 0; i < P.getNumberOfPoints() - 1; i++) {
			line(P.getPoint(i), P.getPoint(i + 1));
		}
	}

	public void drawPolylineEdges(final Collection<WB_PolyLine> polylines) {
		final Iterator<WB_PolyLine> polyItr = polylines.iterator();
		while (polyItr.hasNext()) {
			drawPolylineEdges(polyItr.next());
		}
	}

	public void drawPolylineEdges(final WB_PolyLine P) {
		for (int i = 0; i < P.getNumberOfPoints() - 1; i++) {
			line(P.getPoint(i), P.getPoint(i + 1));
		}
	}

	public void drawPolyLineEmbedded2D(final WB_PolyLine P, final WB_Map2D map) {
		for (int i = 0; i < P.getNumberOfPoints() - 1; i++) {
			drawSegmentEmbedded2D(P.getPoint(i), P.getPoint(i + 1), map);
		}
	}

	public void drawPolyLineMapped(final WB_PolyLine P, final WB_Map map) {
		for (int i = 0; i < P.getNumberOfPoints() - 1; i++) {
			drawSegmentMapped(P.getPoint(i), P.getPoint(i + 1), map);
		}
	}

	public void drawPolyLineUnmapped(final WB_PolyLine P, final WB_Map map) {
		for (int i = 0; i < P.getNumberOfPoints() - 1; i++) {
			drawSegmentUnmapped(P.getPoint(i), P.getPoint(i + 1), map);
		}
	}

	public void drawPolylineVertices(final Collection<WB_PolyLine> polylines, final double d) {
		final Iterator<WB_PolyLine> polyItr = polylines.iterator();
		while (polyItr.hasNext()) {
			drawPolylineVertices(polyItr.next(), d);
		}
	}

	public void drawPolylineVertices(final WB_PolyLine P, final double d) {
		WB_Point v1;
		for (int i = 0; i < P.getNumberOfPoints(); i++) {
			v1 = P.getPoint(i);
			home.pushMatrix();
			translate(v1);
			home.box((float) d);
			home.popMatrix();
		}
	}

	public void drawRay(final WB_Ray R, final double d) {
		home.line((float) R.getOrigin().xd(), (float) R.getOrigin().yd(), (float) R.getOrigin().zd(),
				(float) (R.getOrigin().xd() + d * R.getDirection().xd()),
				(float) (R.getOrigin().yd() + d * R.getDirection().yd()),
				(float) (R.getOrigin().zd() + d * R.getDirection().zd()));
	}

	public void drawRayEmbedded2D(final WB_Ray R, final double d, final WB_Map2D map) {
		drawSegmentEmbedded2D(R.getOrigin(), WB_Point.addMul(R.getOrigin(), d, R.getDirection()), map);
	}

	public void drawRayMapped(final WB_Ray R, final double d, final WB_Map map) {
		drawSegmentMapped(R.getOrigin(), WB_Point.addMul(R.getOrigin(), d, R.getDirection()), map);
	}

	public void drawRayUnmapped(final WB_Ray R, final double d, final WB_Map map) {
		drawSegmentUnmapped(R.getOrigin(), WB_Point.addMul(R.getOrigin(), d, R.getDirection()), map);
	}

	public void drawRing(final WB_Ring P) {
		for (int i = 0, j = P.getNumberOfPoints() - 1; i < P.getNumberOfPoints(); j = i++) {
			line(P.getPoint(j), P.getPoint(i));
		}
	}

	public void drawRingEmbedded2D(final WB_Ring P, final WB_Map2D map) {
		for (int i = 0, j = P.getNumberOfPoints() - 1; i < P.getNumberOfPoints(); j = i++) {
			drawSegmentEmbedded2D(P.getPoint(j), P.getPoint(i), map);
		}
	}

	public void drawRingMapped(final WB_Ring P, final WB_Map map) {
		for (int i = 0, j = P.getNumberOfPoints() - 1; i < P.getNumberOfPoints(); j = i++) {
			drawSegmentMapped(P.getPoint(j), P.getPoint(i), map);
		}
	}

	public void drawRingUnmapped(final WB_Ring P, final WB_Map map) {
		for (int i = 0, j = P.getNumberOfPoints() - 1; i < P.getNumberOfPoints(); j = i++) {
			drawSegmentUnmapped(P.getPoint(j), P.getPoint(i), map);
		}
	}

	public void drawSegment(final Collection<? extends WB_Segment> segments) {
		final Iterator<? extends WB_Segment> segItr = segments.iterator();
		while (segItr.hasNext()) {
			drawSegment(segItr.next());
		}
	}

	public void drawSegment(final List<? extends WB_Coord> points) {
		for (int i = 0; i < points.size(); i += 2) {
			line(points.get(i), points.get(i + 1));
		}
	}

	public void drawSegment(final WB_Coord p, final WB_Coord q) {
		line(p, q);
	}

	public void drawSegment(final int[] segs, final List<? extends WB_Coord> coords) {
		for (int i = 0; i < segs.length; i += 2) {
			line(coords.get(segs[i]), coords.get(segs[i + 1]));
		}
	}

	public void drawSegment(final WB_Segment S) {
		line(S.getOrigin(), S.getEndpoint());
	}

	public void drawSegmentEmbedded2D(final Collection<? extends WB_Segment> segments, final WB_Map2D map) {
		final Iterator<? extends WB_Segment> segItr = segments.iterator();
		while (segItr.hasNext()) {
			drawSegmentEmbedded2D(segItr.next(), map);
		}
	}

	public void drawSegmentEmbedded2D(final WB_Coord p, final WB_Coord q, final WB_Map2D map) {
		home.beginShape();
		vertexEmbedded2D(p, map);
		vertexEmbedded2D(q, map);
		home.endShape();
	}

	public void drawSegmentEmbedded2D(final WB_Segment segment, final WB_Map2D map) {
		drawSegmentEmbedded2D(segment.getOrigin(), segment.getEndpoint(), map);
	}

	public void drawSegmentEmbedded2D(final WB_Segment[] segments, final WB_Map2D map) {
		for (final WB_Segment segment : segments) {
			drawSegmentEmbedded2D(segment, map);
		}
	}

	public void drawSegmentMapped(final Collection<? extends WB_Segment> segments, final WB_Map map) {
		final Iterator<? extends WB_Segment> segItr = segments.iterator();
		while (segItr.hasNext()) {
			drawSegmentMapped(segItr.next(), map);
		}
	}

	public void drawSegmentMapped(final WB_Coord p, final WB_Coord q, final WB_Map map) {
		home.beginShape();
		vertexMapped(p, map);
		vertexMapped(q, map);
		home.endShape();
	}

	public void drawSegmentMapped(final WB_Segment segment, final WB_Map map) {
		drawSegmentMapped(segment.getOrigin(), segment.getEndpoint(), map);
	}

	public void drawSegmentMapped(final WB_Segment[] segments, final WB_Map map) {
		for (final WB_Segment segment : segments) {
			drawSegmentMapped(segment, map);
		}
	}

	public void drawSegmentUnmapped(final Collection<? extends WB_Segment> segments, final WB_Map map) {
		final Iterator<? extends WB_Segment> segItr = segments.iterator();
		while (segItr.hasNext()) {
			drawSegmentUnmapped(segItr.next(), map);
		}
	}

	public void drawSegmentUnmapped(final WB_Coord p, final WB_Coord q, final WB_Map map) {
		home.beginShape();
		vertexUnmapped(p, map);
		vertexUnmapped(q, map);
		home.endShape();
	}

	public void drawSegmentUnmapped(final WB_Segment segment, final WB_Map map) {
		drawSegmentUnmapped(segment.getOrigin(), segment.getEndpoint(), map);
	}

	public void drawSegmentUnmapped(final WB_Segment[] segments, final WB_Map map) {
		for (final WB_Segment segment : segments) {
			drawSegmentUnmapped(segment, map);
		}
	}

	public void drawSimplePolygon(final WB_Polygon P) {
		{
			home.beginShape(PConstants.POLYGON);
			for (int i = 0; i < P.getNumberOfPoints(); i++) {
				vertex(P.getPoint(i));
			}
		}
		home.endShape();
	}

	public void drawTetrahedron(final int[] indices, final List<? extends WB_Coord> points) {
		if (points != null && indices != null) {
			for (int i = 0; i < indices.length; i += 4) {
				drawTetrahedron(points.get(indices[i]), points.get(indices[i + 1]), points.get(indices[i + 2]),
						points.get(indices[i + 3]));
			}
		}
	}

	public void drawTetrahedron(final int[] indices, final WB_Coord[] points) {
		if (points != null && indices != null) {
			for (int i = 0; i < indices.length; i += 4) {
				drawTetrahedron(points[indices[i]], points[indices[i + 1]], points[indices[i + 2]],
						points[indices[i + 3]]);
			}
		}
	}

	public void drawTetrahedron(final WB_Coord p0, final WB_Coord p1, final WB_Coord p2, final WB_Coord p3) {
		home.beginShape(PConstants.TRIANGLES);
		vertex(p0);
		vertex(p1);
		vertex(p2);
		vertex(p1);
		vertex(p0);
		vertex(p3);
		vertex(p2);
		vertex(p1);
		vertex(p3);
		vertex(p0);
		vertex(p2);
		vertex(p3);
		home.endShape();
	}

	public void drawTetrahedron(final WB_Tetrahedron tetra) {
		home.beginShape(PConstants.TRIANGLES);
		vertex(tetra.p1());
		vertex(tetra.p2());
		vertex(tetra.p3());
		vertex(tetra.p2());
		vertex(tetra.p1());
		vertex(tetra.p4());
		vertex(tetra.p3());
		vertex(tetra.p2());
		vertex(tetra.p4());
		vertex(tetra.p1());
		vertex(tetra.p3());
		vertex(tetra.p4());
		home.endShape();
	}

	public void drawTetrahedron(final Collection<? extends WB_Tetrahedron> tetras) {
		for (final WB_Tetrahedron tetra : tetras) {
			home.beginShape(PConstants.TRIANGLES);
			vertex(tetra.p1());
			vertex(tetra.p2());
			vertex(tetra.p3());
			vertex(tetra.p2());
			vertex(tetra.p1());
			vertex(tetra.p4());
			vertex(tetra.p3());
			vertex(tetra.p2());
			vertex(tetra.p4());
			vertex(tetra.p1());
			vertex(tetra.p3());
			vertex(tetra.p4());
			home.endShape();
		}
	}

	public void drawTetrahedron(final WB_Tetrahedron[] tetras) {
		for (final WB_Tetrahedron tetra : tetras) {
			home.beginShape(PConstants.TRIANGLES);
			vertex(tetra.p1());
			vertex(tetra.p2());
			vertex(tetra.p3());
			vertex(tetra.p2());
			vertex(tetra.p1());
			vertex(tetra.p4());
			vertex(tetra.p3());
			vertex(tetra.p2());
			vertex(tetra.p4());
			vertex(tetra.p1());
			vertex(tetra.p3());
			vertex(tetra.p4());
			home.endShape();
		}
	}

	public void drawTree(final WB_OctreeInteger tree) {
		if (tree.getNumNodes() == 0) {
			drawAABB(tree.getBox());
		} else {
			final List<WB_OctreeInteger> nodes = tree.getNodes();
			for (final WB_OctreeInteger node : nodes) {
				drawTree(node);
			}
		}
	}

	public void drawTree(final WB_QuadtreeInteger tree) {
		if (tree.getNumNodes() == 0) {
			drawAABB2D(tree.getBox());
		} else {
			final List<WB_QuadtreeInteger> nodes = tree.getNodes();
			for (final WB_QuadtreeInteger node : nodes) {
				drawTree(node);
			}
		}
	}

	public void drawTriangle(final Collection<? extends WB_Triangle> triangles) {
		final Iterator<? extends WB_Triangle> triItr = triangles.iterator();
		while (triItr.hasNext()) {
			drawTriangle(triItr.next());
		}
	}

	public void drawTriangle(final int[] tri, final List<? extends WB_Coord> points) {
		for (int i = 0; i < tri.length; i += 3) {
			home.beginShape(PConstants.TRIANGLES);
			vertex(points.get(tri[i]));
			vertex(points.get(tri[i + 1]));
			vertex(points.get(tri[i + 2]));
			home.endShape();
		}
	}

	public void drawTriangle(final List<? extends WB_Coord> points) {
		for (int i = 0; i < points.size(); i += 3) {
			home.beginShape(PConstants.TRIANGLES);
			vertex(points.get(i));
			vertex(points.get(i + 1));
			vertex(points.get(i + 2));
			home.endShape();
		}
	}

	public void drawTriangle(final double[] raw) {
		for (int i = 0; i < raw.length;) {
			home.beginShape(PConstants.TRIANGLES);
			vertex(raw[i++], raw[i++], raw[i++]);
			vertex(raw[i++], raw[i++], raw[i++]);
			vertex(raw[i++], raw[i++], raw[i++]);
			home.endShape();
		}
	}

	public void drawTriangle(final int[] tri, final WB_Coord[] points) {
		for (int i = 0; i < tri.length; i += 3) {
			home.beginShape(PConstants.TRIANGLES);
			vertex(points[tri[i]]);
			vertex(points[tri[i + 1]]);
			vertex(points[tri[i + 2]]);
			home.endShape();
		}
	}

	public void drawTriangle(final WB_Coord p1, final WB_Coord p2, final WB_Coord p3) {
		home.beginShape(PConstants.TRIANGLES);
		vertex(p1);
		vertex(p2);
		vertex(p3);
		home.endShape();
	}

	public void drawTriangle(final WB_Triangle triangle) {
		home.beginShape();
		vertex(triangle.p1());
		vertex(triangle.p2());
		vertex(triangle.p3());
		home.endShape(PConstants.CLOSE);
	}

	public void drawTriangleEdges(final Collection<? extends WB_Triangle> triangles) {
		final Iterator<? extends WB_Triangle> triItr = triangles.iterator();
		while (triItr.hasNext()) {
			drawTriangleEdges(triItr.next());
		}
	}

	public void drawTriangleEdges(final WB_Triangle triangle) {
		line(triangle.p1(), triangle.p2());
		line(triangle.p3(), triangle.p2());
		line(triangle.p1(), triangle.p3());
	}

	public void drawTriangleEdgesEmbedded2D(final Collection<? extends WB_Triangle> triangles, final WB_Map2D map) {
		final Iterator<? extends WB_Triangle> triItr = triangles.iterator();
		while (triItr.hasNext()) {
			drawTriangleEdgesEmbedded2D(triItr.next(), map);
		}
	}

	public void drawTriangleEdgesEmbedded2D(final WB_Triangle triangle, final WB_Map2D map) {
		drawSegmentEmbedded2D(triangle.p1(), triangle.p2(), map);
		drawSegmentEmbedded2D(triangle.p2(), triangle.p3(), map);
		drawSegmentEmbedded2D(triangle.p3(), triangle.p1(), map);
	}

	public void drawTriangleEdgesEmbedded2D(final WB_Triangle[] triangles, final WB_Map2D map) {
		for (final WB_Triangle triangle : triangles) {
			drawTriangleEdgesEmbedded2D(triangle, map);
		}
	}

	public void drawTriangleEdgesMapped(final Collection<? extends WB_Triangle> triangles, final WB_Map map) {
		final Iterator<? extends WB_Triangle> triItr = triangles.iterator();
		while (triItr.hasNext()) {
			drawTriangleEdgesMapped(triItr.next(), map);
		}
	}

	public void drawTriangleEdgesMapped(final WB_Triangle triangle, final WB_Map map) {
		drawSegmentMapped(triangle.p1(), triangle.p2(), map);
		drawSegmentMapped(triangle.p2(), triangle.p3(), map);
		drawSegmentMapped(triangle.p3(), triangle.p1(), map);
	}

	public void drawTriangleEdgesMapped(final WB_Triangle[] triangles, final WB_Map map) {
		for (final WB_Triangle triangle : triangles) {
			drawTriangleEdgesMapped(triangle, map);
		}
	}

	public void drawTriangleEdgesUnmapped(final Collection<? extends WB_Triangle> triangles, final WB_Map map) {
		final Iterator<? extends WB_Triangle> triItr = triangles.iterator();
		while (triItr.hasNext()) {
			drawTriangleEdgesUnmapped(triItr.next(), map);
		}
	}

	public void drawTriangleEdgesUnmapped(final WB_Triangle triangle, final WB_Map map) {
		drawSegmentUnmapped(triangle.p1(), triangle.p2(), map);
		drawSegmentUnmapped(triangle.p2(), triangle.p3(), map);
		drawSegmentUnmapped(triangle.p3(), triangle.p1(), map);
	}

	public void drawTriangleEdgesUnmapped(final WB_Triangle[] triangles, final WB_Map map) {
		for (final WB_Triangle triangle : triangles) {
			drawTriangleEdgesUnmapped(triangle, map);
		}
	}

	public void drawTriangleEmbedded2D(final Collection<? extends WB_Triangle> triangles, final WB_Map2D map) {
		final Iterator<? extends WB_Triangle> triItr = triangles.iterator();
		while (triItr.hasNext()) {
			drawTriangleEmbedded2D(triItr.next(), map);
		}
	}

	public void drawTriangleEmbedded2D(final int[] tri, final List<? extends WB_Coord> points, final WB_Map2D map) {
		for (int i = 0; i < tri.length; i += 3) {
			home.beginShape(PConstants.TRIANGLES);
			vertexEmbedded2D(points.get(tri[i]), map);
			vertexEmbedded2D(points.get(tri[i + 1]), map);
			vertexEmbedded2D(points.get(tri[i + 2]), map);
			home.endShape();
		}
	}

	public void drawTriangleEmbedded2D(final int[] tri, final WB_Coord[] points, final WB_Map2D map) {
		for (int i = 0; i < tri.length; i += 3) {
			home.beginShape(PConstants.TRIANGLES);
			vertexEmbedded2D(points[tri[i]], map);
			vertexEmbedded2D(points[tri[i + 1]], map);
			vertexEmbedded2D(points[tri[i + 2]], map);
			home.endShape();
		}
	}

	public void drawTriangleEmbedded2D(final WB_Coord p1, final WB_Coord p2, final WB_Coord p3, final WB_Map2D map) {
		home.beginShape(PConstants.TRIANGLES);
		vertexEmbedded2D(p1, map);
		vertexEmbedded2D(p2, map);
		vertexEmbedded2D(p3, map);
		home.endShape();
	}

	public void drawTriangleEmbedded2D(final WB_Triangle T, final WB_Map2D map) {
		home.beginShape(PConstants.TRIANGLES);
		vertexEmbedded2D(T.p1(), map);
		vertexEmbedded2D(T.p2(), map);
		vertexEmbedded2D(T.p3(), map);
		home.endShape();
	}

	public void drawTriangleEmbedded2D(final WB_Triangle[] triangles, final WB_Map2D map) {
		for (final WB_Triangle triangle : triangles) {
			drawTriangleEmbedded2D(triangle, map);
		}
	}

	public void drawTriangleMapped(final Collection<? extends WB_Triangle> triangles, final WB_Map map) {
		final Iterator<? extends WB_Triangle> triItr = triangles.iterator();
		while (triItr.hasNext()) {
			drawTriangleMapped(triItr.next(), map);
		}
	}

	public void drawTriangleMapped(final int[] tri, final List<? extends WB_Coord> points, final WB_Map map) {
		for (int i = 0; i < tri.length; i += 3) {
			home.beginShape(PConstants.TRIANGLES);
			vertexMapped(points.get(tri[i]), map);
			vertexMapped(points.get(tri[i + 1]), map);
			vertexMapped(points.get(tri[i + 2]), map);
			home.endShape();
		}
	}

	public void drawTriangleMapped(final int[] tri, final WB_Coord[] points, final WB_Map map) {
		for (int i = 0; i < tri.length; i += 3) {
			home.beginShape(PConstants.TRIANGLES);
			vertexMapped(points[tri[i]], map);
			vertexMapped(points[tri[i + 1]], map);
			vertexMapped(points[tri[i + 2]], map);
			home.endShape();
		}
	}

	public void drawTriangleMapped(final WB_Coord p1, final WB_Coord p2, final WB_Coord p3, final WB_Map map) {
		home.beginShape(PConstants.TRIANGLES);
		vertexMapped(p1, map);
		vertexMapped(p2, map);
		vertexMapped(p3, map);
		home.endShape();
	}

	public void drawTriangleMapped(final WB_Triangle T, final WB_Map map) {
		home.beginShape(PConstants.TRIANGLES);
		vertexMapped(T.p1(), map);
		vertexMapped(T.p2(), map);
		vertexMapped(T.p3(), map);
		home.endShape();
	}

	public void drawTriangleMapped(final WB_Triangle[] triangles, final WB_Map map) {
		for (final WB_Triangle triangle : triangles) {
			drawTriangleMapped(triangle, map);
		}
	}

	public void drawTriangleUnmapped(final Collection<? extends WB_Triangle> triangles, final WB_Map map) {
		final Iterator<? extends WB_Triangle> triItr = triangles.iterator();
		while (triItr.hasNext()) {
			drawTriangleUnmapped(triItr.next(), map);
		}
	}

	public void drawTriangleUnmapped(final int[] tri, final List<? extends WB_Coord> points, final WB_Map map) {
		for (int i = 0; i < tri.length; i += 3) {
			home.beginShape(PConstants.TRIANGLES);
			vertexUnmapped(points.get(tri[i]), map);
			vertexUnmapped(points.get(tri[i + 1]), map);
			vertexUnmapped(points.get(tri[i + 2]), map);
			home.endShape();
		}
	}

	public void drawTriangleUnmapped(final int[] tri, final WB_Coord[] points, final WB_Map map) {
		for (int i = 0; i < tri.length; i += 3) {
			home.beginShape(PConstants.TRIANGLES);
			vertexUnmapped(points[tri[i]], map);
			vertexUnmapped(points[tri[i + 1]], map);
			vertexUnmapped(points[tri[i + 2]], map);
			home.endShape();
		}
	}

	public void drawTriangleUnmapped(final WB_Coord p1, final WB_Coord p2, final WB_Coord p3, final WB_Map map) {
		home.beginShape(PConstants.TRIANGLES);
		vertexUnmapped(p1, map);
		vertexUnmapped(p2, map);
		vertexUnmapped(p3, map);
		home.endShape();
	}

	public void drawTriangleUnmapped(final WB_Triangle T, final WB_Map map) {
		home.beginShape(PConstants.TRIANGLES);
		vertexUnmapped(T.p1(), map);
		vertexUnmapped(T.p2(), map);
		vertexUnmapped(T.p3(), map);
		home.endShape();
	}

	public void drawTriangleUnmapped(final WB_Triangle[] triangles, final WB_Map map) {
		for (final WB_Triangle triangle : triangles) {
			drawTriangleUnmapped(triangle, map);
		}
	}

	public void drawTriangulation(final WB_Triangulation2D tri, final List<? extends WB_Coord> points) {
		final int[] triangles = tri.getTriangles();
		home.beginShape(PConstants.TRIANGLES);
		for (int i = 0; i < triangles.length; i += 3) {
			vertex(points.get(triangles[i]));
			vertex(points.get(triangles[i + 1]));
			vertex(points.get(triangles[i + 2]));
		}
		home.endShape();
	}

	public void drawTriangulation(final WB_Triangulation2D tri, final WB_Coord[] points) {
		final int[] triangles = tri.getTriangles();
		home.beginShape(PConstants.TRIANGLES);
		for (int i = 0; i < triangles.length; i += 3) {
			vertex(points[triangles[i]]);
			vertex(points[triangles[i + 1]]);
			vertex(points[triangles[i + 2]]);
		}
		home.endShape();
	}

	public void drawTriangulation(final WB_Triangulation3D tri, final List<? extends WB_Coord> points) {
		drawTetrahedron(tri.getTetrahedra(), points);
	}

	public void drawTriangulation(final WB_Triangulation3D tri, final WB_Coord[] points) {
		drawTetrahedron(tri.getTetrahedra(), points);
	}

	public void drawTriangulationEdges(final WB_Triangulation2D tri, final List<? extends WB_Coord> points) {
		final int[] edges = tri.getEdges();
		for (int i = 0; i < edges.length; i += 2) {
			drawSegment(points.get(edges[i]), points.get(edges[i + 1]));
		}
	}

	public void drawTriangulationEdges(final WB_Triangulation2D tri, final WB_Coord[] points) {
		final int[] edges = tri.getEdges();
		for (int i = 0; i < edges.length; i += 2) {
			drawSegment(points[edges[i]], points[edges[i + 1]]);
		}
	}

	public void drawTriangulationEdgesEmbedded2D(final WB_Triangulation2D tri, final List<? extends WB_Coord> points,
			final WB_Map2D map) {
		final int[] edges = tri.getEdges();
		for (int i = 0; i < edges.length; i += 2) {
			drawSegmentEmbedded2D(points.get(edges[i]), points.get(edges[i + 1]), map);
		}
	}

	public void drawTriangulationEdgesEmbedded2D(final WB_Triangulation2D tri, final WB_Coord[] points,
			final WB_Map2D map) {
		final int[] edges = tri.getEdges();
		for (int i = 0; i < edges.length; i += 2) {
			drawSegmentEmbedded2D(points[edges[i]], points[edges[i + 1]], map);
		}
	}

	public void drawTriangulationEdgesMapped(final WB_Triangulation2D tri, final List<? extends WB_Coord> points,
			final WB_Map map) {
		final int[] edges = tri.getEdges();
		for (int i = 0; i < edges.length; i += 2) {
			drawSegmentMapped(points.get(edges[i]), points.get(edges[i + 1]), map);
		}
	}

	public void drawTriangulationEdgesMapped(final WB_Triangulation2D tri, final WB_Coord[] points, final WB_Map map) {
		final int[] edges = tri.getEdges();
		for (int i = 0; i < edges.length; i += 2) {
			drawSegmentMapped(points[edges[i]], points[edges[i + 1]], map);
		}
	}

	public void drawTriangulationEdgesUnmapped(final WB_Triangulation2D tri, final List<? extends WB_Coord> points,
			final WB_Map map) {
		final int[] edges = tri.getEdges();
		for (int i = 0; i < edges.length; i += 2) {
			drawSegmentUnmapped(points.get(edges[i]), points.get(edges[i + 1]), map);
		}
	}

	public void drawTriangulationEdgesUnmapped(final WB_Triangulation2D tri, final WB_Coord[] points,
			final WB_Map map) {
		final int[] edges = tri.getEdges();
		for (int i = 0; i < edges.length; i += 2) {
			drawSegmentUnmapped(points[edges[i]], points[edges[i + 1]], map);
		}
	}

	public void drawTriangulationEmbedded2D(final WB_Triangulation2D tri, final List<? extends WB_Coord> points,
			final WB_Map2D map) {
		final int[] triangles = tri.getTriangles();
		home.beginShape(PConstants.TRIANGLES);
		for (int i = 0; i < triangles.length; i += 3) {
			vertexEmbedded2D(points.get(triangles[i]), map);
			vertexEmbedded2D(points.get(triangles[i + 1]), map);
			vertexEmbedded2D(points.get(triangles[i + 2]), map);
		}
		home.endShape();
	}

	public void drawTriangulationEmbedded2D(final WB_Triangulation2D tri, final WB_Coord[] points, final WB_Map2D map) {
		final int[] triangles = tri.getTriangles();
		home.beginShape(PConstants.TRIANGLES);
		for (int i = 0; i < triangles.length; i += 3) {
			vertexEmbedded2D(points[triangles[i]], map);
			vertexEmbedded2D(points[triangles[i + 1]], map);
			vertexEmbedded2D(points[triangles[i + 2]], map);
		}
		home.endShape();
	}

	public void drawTriangulationMapped(final WB_Triangulation2D tri, final List<? extends WB_Coord> points,
			final WB_Map map) {
		final int[] triangles = tri.getTriangles();
		home.beginShape(PConstants.TRIANGLES);
		for (int i = 0; i < triangles.length; i += 3) {
			vertexMapped(points.get(triangles[i]), map);
			vertexMapped(points.get(triangles[i + 1]), map);
			vertexMapped(points.get(triangles[i + 2]), map);
		}
		home.endShape();
	}

	public void drawTriangulationMapped(final WB_Triangulation2D tri, final WB_Coord[] points, final WB_Map map) {
		final int[] triangles = tri.getTriangles();
		home.beginShape(PConstants.TRIANGLES);
		for (int i = 0; i < triangles.length; i += 3) {
			vertexMapped(points[triangles[i]], map);
			vertexMapped(points[triangles[i + 1]], map);
			vertexMapped(points[triangles[i + 2]], map);
		}
		home.endShape();
	}

	public void drawTriangulationUnmapped(final WB_Triangulation2D tri, final List<? extends WB_Coord> points,
			final WB_Map map) {
		final int[] triangles = tri.getTriangles();
		home.beginShape(PConstants.TRIANGLES);
		for (int i = 0; i < triangles.length; i += 3) {
			vertexUnmapped(points.get(triangles[i]), map);
			vertexUnmapped(points.get(triangles[i + 1]), map);
			vertexUnmapped(points.get(triangles[i + 2]), map);
		}
		home.endShape();
	}

	public void drawTriangulationUnmapped(final WB_Triangulation2D tri, final WB_Coord[] points, final WB_Map map) {
		final int[] triangles = tri.getTriangles();
		home.beginShape(PConstants.TRIANGLES);
		for (int i = 0; i < triangles.length; i += 3) {
			vertexUnmapped(points[triangles[i]], map);
			vertexUnmapped(points[triangles[i + 1]], map);
			vertexUnmapped(points[triangles[i + 2]], map);
		}
		home.endShape();
	}

	public void drawVector(final WB_Coord p, final WB_Coord v, final double r) {
		home.pushMatrix();
		translate(p);
		home.line(0f, 0f, 0f, (float) (r * v.xd()), (float) (r * v.yd()), (float) (r * v.zd()));
		home.popMatrix();
	}

	public void drawVectorEmbedded2D(final WB_Coord p, final WB_Coord v, final double r, final WB_Map2D map) {
		drawSegmentEmbedded2D(p, WB_Point.addMul(p, r, v), map);
	}

	public void drawVectorMapped(final WB_Coord p, final WB_Coord v, final double r, final WB_Map map) {
		drawSegmentMapped(p, WB_Point.addMul(p, r, v), map);
	}

	public void drawVectorUnmapped(final WB_Coord p, final WB_Coord v, final double r, final WB_Map map) {
		drawSegmentUnmapped(p, WB_Point.addMul(p, r, v), map);
	}

	public void drawVertex(final long key, final HE_Mesh mesh, final double d) {
		final HE_Vertex v = mesh.getVertexWithKey(key);
		if (v != null && v.isVisible()) {
			home.pushMatrix();
			translate(v);
			home.box((float) d);
			home.popMatrix();
		}
	}

	public void drawVertex(final HE_Vertex v, final double d) {
		if (v != null && v.isVisible()) {
			home.pushMatrix();
			translate(v);
			home.box((float) d);
			home.popMatrix();
		}
	}

	public void drawVertexNormals(final HE_HalfedgeStructure mesh, final double d) {
		final Iterator<HE_Vertex> vItr = mesh.vItr();
		WB_Coord vn;
		HE_Vertex v;
		while (vItr.hasNext()) {
			v = vItr.next();
			vn = HE_MeshOp.getVertexNormal(v);
			drawVector(v, vn, d);
		}
	}

	public void drawVertices(final HE_HalfedgeStructure mesh) {
		HE_Vertex v;
		final Iterator<HE_Vertex> vItr = mesh.vItr();
		while (vItr.hasNext()) {
			v = vItr.next();
			if (v.isVisible()) {
				home.point(v.xf(), v.yf(), v.zf());
			}
		}
	}

	public void drawVertices(final HE_HalfedgeStructure mesh, final double d) {
		HE_Vertex v;
		final Iterator<HE_Vertex> vItr = mesh.vItr();
		while (vItr.hasNext()) {
			v = vItr.next();
			if (v.isVisible()) {
				home.pushMatrix();
				translate(v);
				home.box((float) d);
				home.popMatrix();
			}
		}
	}

	public PGraphicsOpenGL getHome() {
		return home;
	}

	class PickingRay {
		double iwidth = 0.0; // inverse width
		double iheight = 0.0; // inverse height
		PMatrix3D unprojection = new PMatrix3D();

		void getContext(final PGraphicsOpenGL home) {
			unprojection.set(home.projection);
			unprojection.apply(home.modelview);
			unprojection.invert();
			iwidth = 1.0 / home.width;
			iheight = 1.0 / home.height;
		}

		WB_Coord unproject(final double x, final double y, final double z) {
			// Normalized screen coordinates
			// x=mouseX to -1..1
			// y=mouseY to -1..1
			// z=0 if on near clipping plane, z=1 if on far clipping plane to
			// -1..1
			final double normScreenX = 2.0 * (iwidth * x - 0.5);
			final double normScreenY = 2.0 * (0.5 - iheight * y);
			final double normScreenZ = 2.0 * (constrain(z, 0, 1) - 0.5);
			// normScreenW=1.0
			// Homogeneous coordinates in world space x',y',z',w' = normalized
			// screen coordinates * inverse projection matrix
			// Only possible if w' is not zero
			final double wPrime = normScreenX * unprojection.m30 + normScreenY * unprojection.m31
					+ unprojection.m32 * normScreenZ + unprojection.m33;
			// Homogeneous coordinates to Cartesian coordinates: x=x'/w',
			// y=y'/w',z=z'/w'
			if (Math.abs(wPrime) < 1e-12) { // "Point in infinity"
				return null;
			}
			final double iw = 1.0 / wPrime;
			final double xPrime = normScreenX * unprojection.m00 + normScreenY * unprojection.m01
					+ normScreenZ * unprojection.m02 + unprojection.m03;
			final double yPrime = normScreenX * unprojection.m10 + normScreenY * unprojection.m11
					+ normScreenZ * unprojection.m12 + unprojection.m13;
			final double zPrime = normScreenX * unprojection.m20 + normScreenY * unprojection.m21
					+ normScreenZ * unprojection.m22 + unprojection.m23;
			return new WB_Point(xPrime * iw, yPrime * iw, zPrime * iw);
		}

		WB_Ray getPickingRay(final double x, final double y) {
			return geometryfactory.createRayThroughPoints(unproject(x, y, 0.0), unproject(x, y, 1.0));
		}

		double constrain(final double val, final double min, final double max) {
			if (val <= min) {
				return min;
			}
			if (val >= max) {
				return max;
			}
			return val;
		}
	}

	final PickingRay PRay = new PickingRay();

	public WB_Ray getPickingRay(final double x, final double y) {
		PRay.getContext(home);
		final WB_Ray ray = PRay.getPickingRay(x, y);
		return ray;
	}

	public HE_Face pickClosestFace(final HE_Mesh mesh, final double x, final double y) {
		final WB_Ray mouseRay3d = getPickingRay(x, y);
		final HE_FaceLineIntersection p = HE_MeshOp.getClosestIntersection(mesh, mouseRay3d);
		return p == null ? null : p.getFace();
	}

	public HE_Face pickClosestFace(final WB_AABBTree3D meshtree, final double x, final double y) {
		final WB_Ray mouseRay3d = getPickingRay(x, y);
		final HE_FaceLineIntersection p = HE_MeshOp.getClosestIntersection(meshtree, mouseRay3d);
		return p == null ? null : p.getFace();
	}

	public HE_Halfedge pickEdge(final HE_Mesh mesh, final double x, final double y) {
		final WB_Ray mouseRay3d = getPickingRay(x, y);
		final HE_FaceLineIntersection p = HE_MeshOp.getClosestIntersection(mesh, mouseRay3d);
		if (p == null) {
			return null;
		}
		final HE_Face f = p.getFace();
		final HE_FaceEdgeCirculator fec = f.feCrc();
		HE_Halfedge trial;
		HE_Halfedge closest = null;
		double d2 = 0;
		double d2min = Double.MAX_VALUE;
		while (fec.hasNext()) {
			trial = fec.next();
			d2 = WB_GeometryOp3D.getDistanceToSegment3D(p.getPoint(), trial.getStartVertex(), trial.getEndVertex());
			if (d2 < d2min) {
				d2min = d2;
				closest = trial;
			}
		}
		return closest;
	}

	public HE_Halfedge pickEdge(final WB_AABBTree3D tree, final double x, final double y) {
		final WB_Ray mouseRay3d = getPickingRay(x, y);
		final HE_FaceLineIntersection p = HE_MeshOp.getClosestIntersection(tree, mouseRay3d);
		if (p == null) {
			return null;
		}
		final HE_Face f = p.getFace();
		final HE_FaceEdgeCirculator fec = f.feCrc();
		HE_Halfedge trial;
		HE_Halfedge closest = null;
		double d2 = 0;
		double d2min = Double.MAX_VALUE;
		while (fec.hasNext()) {
			trial = fec.next();
			d2 = WB_GeometryOp3D.getDistanceToSegment3D(p.getPoint(), trial.getStartVertex(), trial.getEndVertex());
			if (d2 < d2min) {
				d2min = d2;
				closest = trial;
			}
		}
		return closest;
	}

	public List<HE_Face> pickFaces(final HE_Mesh mesh, final double x, final double y) {
		final WB_Ray mouseRay3d = getPickingRay(x, y);
		final List<HE_FaceLineIntersection> p = HE_MeshOp.getIntersection(mesh, mouseRay3d);
		final List<HE_Face> result = new ArrayList<>();
		for (final HE_FaceLineIntersection fi : p) {
			result.add(fi.getFace());
		}
		return result;
	}

	public List<HE_Face> pickFaces(final WB_AABBTree3D meshtree, final double x, final double y) {
		final WB_Ray mouseRay3d = getPickingRay(x, y);
		final List<HE_FaceLineIntersection> p = HE_MeshOp.getIntersection(meshtree, mouseRay3d);
		final List<HE_Face> result = new ArrayList<>();
		for (final HE_FaceLineIntersection fi : p) {
			result.add(fi.getFace());
		}
		return result;
	}

	public HE_Face pickFurthestFace(final HE_Mesh mesh, final double x, final double y) {
		final WB_Ray mouseRay3d = getPickingRay(x, y);
		final HE_FaceLineIntersection p = HE_MeshOp.getFurthestIntersection(mesh, mouseRay3d);
		return p == null ? null : p.getFace();
	}

	public HE_Face pickFurthestFace(final WB_AABBTree3D meshtree, final double x, final double y) {
		final WB_Ray mouseRay3d = getPickingRay(x, y);
		final HE_FaceLineIntersection p = HE_MeshOp.getFurthestIntersection(meshtree, mouseRay3d);
		return p == null ? null : p.getFace();
	}

	public HE_Vertex pickVertex(final HE_Mesh mesh, final double x, final double y) {
		final WB_Ray mouseRay3d = getPickingRay(x, y);
		final HE_FaceLineIntersection p = HE_MeshOp.getClosestIntersection(mesh, mouseRay3d);
		if (p == null) {
			return null;
		}
		final HE_Face f = p.getFace();
		final HE_FaceVertexCirculator fvc = f.fvCrc();
		HE_Vertex trial;
		HE_Vertex closest = null;
		double d2 = 0;
		double d2min = Double.MAX_VALUE;
		while (fvc.hasNext()) {
			trial = fvc.next();
			d2 = trial.getPosition().getSqDistance(p.getPoint());
			if (d2 < d2min) {
				d2min = d2;
				closest = trial;
			}
		}
		return closest;
	}

	public HE_Vertex pickVertex(final WB_AABBTree3D tree, final double x, final double y) {
		final WB_Ray mouseRay3d = getPickingRay(x, y);
		final HE_FaceLineIntersection p = HE_MeshOp.getClosestIntersection(tree, mouseRay3d);
		if (p == null) {
			return null;
		}
		final HE_Face f = p.getFace();
		final HE_FaceVertexCirculator fvc = f.fvCrc();
		HE_Vertex trial;
		HE_Vertex closest = null;
		double d2 = 0;
		double d2min = Double.MAX_VALUE;
		while (fvc.hasNext()) {
			trial = fvc.next();
			d2 = trial.getPosition().getSqDistance(p.getPoint());
			if (d2 < d2min) {
				d2min = d2;
				closest = trial;
			}
		}
		return closest;
	}

	public void setFaceColorFromTexture(final HE_Mesh mesh, final PImage texture) {
		final HE_FaceIterator fitr = mesh.fItr();
		HE_Face f;
		HE_Halfedge he;
		HE_TextureCoordinate uvw;
		while (fitr.hasNext()) {
			f = fitr.next();
			final HE_FaceHalfedgeInnerCirculator fhc = f.fheiCrc();
			final WB_Point p = new WB_Point();
			int id = 0;
			while (fhc.hasNext()) {
				he = fhc.next();
				uvw = he.getUVW();
				p.addSelf(uvw.ud(), uvw.vd(), 0);
				id++;
			}
			p.divSelf(id);
			f.setColor(getColorFromPImage(p.xd(), p.yd(), texture));
		}
	}

	public void setHalfedgeColorFromTexture(final HE_Mesh mesh, final PImage texture) {
		final HE_FaceIterator fitr = mesh.fItr();
		HE_Face f;
		HE_Halfedge he;
		HE_TextureCoordinate p;
		while (fitr.hasNext()) {
			f = fitr.next();
			final HE_FaceHalfedgeInnerCirculator fhec = f.fheiCrc();
			while (fhec.hasNext()) {
				he = fhec.next();
				p = he.getUVW();
				he.setColor(getColorFromPImage(p.ud(), p.vd(), texture));
			}
		}
	}

	public void setVertexColorFromTexture(final HE_Mesh mesh, final PImage texture) {
		final HE_VertexIterator vitr = mesh.vItr();
		HE_Vertex v;
		HE_TextureCoordinate p;
		while (vitr.hasNext()) {
			v = vitr.next();
			p = v.getHalfedge().getUVW();
			v.setColor(getColorFromPImage(p.ud(), p.vd(), texture));
		}
	}

	public void translate(final WB_Coord p) {
		home.translate(p.xf(), p.yf(), p.zf());
	}

	public void vertex(final WB_Coord p) {
		home.vertex(p.xf(), p.yf(), p.zf());
	}

	public void vertex(final double x, final double y, final double z) {
		home.vertex((float) x, (float) y, (float) z);
	}

	public void normal(final WB_Coord n) {
		home.normal(n.xf(), n.yf(), n.zf());
	}

	public void vertexEmbedded2D(final WB_Coord p, final WB_Map2D map) {
		final WB_Point q = new WB_Point();
		map.mapPoint3D(p, q);
		map.unmapPoint2D(q, q);
		vertex(q);
	}

	public void vertexMapped(final WB_Coord p, final WB_Map map) {
		final WB_Point q = new WB_Point();
		map.mapPoint3D(p, q);
		vertex(q);
	}

	public void vertexUnmapped(final WB_Coord p, final WB_Map map) {
		final WB_Point q = new WB_Point();
		map.mapPoint3D(p, q);
		vertex(q);
	}

	private void line(final WB_Coord p, final WB_Coord q) {
		home.beginShape(PConstants.LINES);
		vertex(p);
		vertex(q);
		home.endShape();
	}

	public void line(final double x1, final double y1, final double z1, final double x2, final double y2,
			final double z2) {
		home.beginShape(PConstants.LINES);
		vertex(x1, y1, z1);
		vertex(x2, y2, z2);
		home.endShape();
	}

	public void drawVerticesWithInternalLabel(final int label, final HE_HalfedgeStructure mesh, final double d) {
		final Iterator<HE_Vertex> vItr = mesh.vItr();
		HE_Vertex v;
		while (vItr.hasNext()) {
			v = vItr.next();
			if (v.isVisible() && v.getInternalLabel() == label) {
				drawVertex(v, d);
			}
		}
	}

	public void drawVerticesWithLabel(final int label, final HE_HalfedgeStructure mesh, final double d) {
		final Iterator<HE_Vertex> vItr = mesh.vItr();
		HE_Vertex v;
		while (vItr.hasNext()) {
			v = vItr.next();
			if (v.isVisible() && v.getLabel() == label) {
				drawVertex(v, d);
			}
		}
	}

	public void drawQuad(final WB_Quad quad) {
		home.beginShape();
		vertex(quad.getP1());
		vertex(quad.getP2());
		vertex(quad.getP3());
		vertex(quad.getP4());
		home.endShape(PConstants.CLOSE);
	}

	public void drawQuad(final Collection<? extends WB_Quad> quads) {
		final Iterator<? extends WB_Quad> qItr = quads.iterator();
		while (qItr.hasNext()) {
			drawQuad(qItr.next());
		}
	}

	public void drawPentagon(final WB_Pentagon pentagon) {
		home.beginShape();
		vertex(pentagon.getP1());
		vertex(pentagon.getP2());
		vertex(pentagon.getP3());
		vertex(pentagon.getP4());
		vertex(pentagon.getP5());
		home.endShape(PConstants.CLOSE);
	}

	public void drawPentagon(final Collection<? extends WB_Pentagon> pentagons) {
		final Iterator<? extends WB_Pentagon> pItr = pentagons.iterator();
		while (pItr.hasNext()) {
			drawPentagon(pItr.next());
		}
	}

	public void drawHexagon(final WB_Hexagon hexagon) {
		home.beginShape();
		vertex(hexagon.getP1());
		vertex(hexagon.getP2());
		vertex(hexagon.getP3());
		vertex(hexagon.getP4());
		vertex(hexagon.getP5());
		vertex(hexagon.getP6());
		home.endShape(PConstants.CLOSE);
	}

	public void drawHexagon(final Collection<? extends WB_Hexagon> hexagons) {
		final Iterator<? extends WB_Hexagon> hItr = hexagons.iterator();
		while (hItr.hasNext()) {
			drawHexagon(hItr.next());
		}
	}

	public void drawOctagon(final WB_Octagon octagon) {
		home.beginShape();
		vertex(octagon.getP1());
		vertex(octagon.getP2());
		vertex(octagon.getP3());
		vertex(octagon.getP4());
		vertex(octagon.getP5());
		vertex(octagon.getP6());
		vertex(octagon.getP7());
		vertex(octagon.getP8());
		home.endShape(PConstants.CLOSE);
	}

	public void drawOctagon(final Collection<? extends WB_Octagon> octagons) {
		final Iterator<? extends WB_Octagon> oItr = octagons.iterator();
		while (oItr.hasNext()) {
			drawOctagon(oItr.next());
		}
	}

	@Override
	public void drawGizmo(final double d) {
		home.pushStyle();
		home.noFill();
		home.strokeWeight(2.0f);
		home.stroke(255, 0, 0);
		home.line(0f, 0f, 0f, (float) d, 0, 0);
		home.stroke(0, 255, 0);
		home.line(0f, 0f, 0f, 0, (float) d, 0);
		home.stroke(0, 0, 255);
		home.line(0f, 0f, 0f, 0, 0, (float) d);
		home.popStyle();
	}

	public void drawGizmo(final double x, final double y, final double z) {
		home.pushStyle();
		home.noFill();
		home.strokeWeight(2.0f);
		home.stroke(255, 0, 0);
		home.line(0f, 0f, 0f, (float) x, 0, 0);
		home.stroke(0, 255, 0);
		home.line(0f, 0f, 0f, 0, (float) y, 0);
		home.stroke(0, 0, 255);
		home.line(0f, 0f, 0f, 0, 0, (float) z);
		home.popStyle();
	}

	public void drawGrid(final int color, final WB_Plane P, final double W, final double H, final int divW,
			final int divH) {
		home.pushStyle();
		home.noFill();
		home.stroke(color);
		final double dW = W / divW;
		final double dH = H / divH;
		for (int i = 0; i <= divW; i++) {
			drawSegment(P.extractPoint(-W / 2 + i * dW, -H / 2), P.extractPoint(-W / 2 + i * dW, H / 2));
		}
		for (int i = 0; i <= divH; i++) {
			drawSegment(P.extractPoint(-W / 2, i * dH - H / 2), P.extractPoint(W / 2, i * dH - H / 2));
		}
		home.popStyle();
	}

	public void drawTriangle(final WB_TriangleFactory triangleGenerator) {
		final int[] tri = triangleGenerator.getTriangles();
		final WB_CoordCollection points = triangleGenerator.getPoints();
		for (int i = 0; i < tri.length; i += 3) {
			home.beginShape(PConstants.TRIANGLES);
			vertex(points.get(tri[i]));
			vertex(points.get(tri[i + 1]));
			vertex(points.get(tri[i + 2]));
			home.endShape();
		}
	}

	public void drawBinaryGridOutline3D(final WB_BinaryGrid3D grid) {
		home.pushMatrix();
		translate(grid.getMin());
		drawXEdges(grid);
		drawYEdges(grid);
		drawZEdges(grid);
		home.popMatrix();
	}

	public void drawBinaryGrid3D(final WB_BinaryGrid3D grid) {
		home.pushMatrix();
		translate(grid.getMin());
		drawXFaces(grid);
		drawYFaces(grid);
		drawZFaces(grid);
		home.popMatrix();
	}

	public void drawBinaryGrid3D(final WB_BinaryGrid3D grid, final int cx, final int cy, final int cz) {
		home.pushStyle();
		home.pushMatrix();
		translate(grid.getMin());
		home.fill(cx);
		drawXFaces(grid);
		home.fill(cy);
		drawYFaces(grid);
		home.fill(cz);
		drawZFaces(grid);
		home.popMatrix();
		home.popStyle();
	}

	public void drawBinaryGridOutline3D(final WB_BinaryGrid3D grid, final int cx, final int cy, final int cz) {
		home.pushStyle();
		home.pushMatrix();
		translate(grid.getMin());
		home.stroke(cx);
		drawXEdges(grid);
		home.stroke(cy);
		drawYEdges(grid);
		home.stroke(cz);
		drawZEdges(grid);
		home.popMatrix();
		home.popStyle();
	}

	void drawXEdges(final WB_BinaryGrid3D grid) {
		int val00, valm0, valmm, val0m, sum;
		double x, y, z;
		for (int i = grid.lx(); i < grid.ux(); i++) {
			x = i * grid.getDX();
			for (int j = grid.ly(); j <= grid.uy(); j++) {
				y = j * grid.getDY();
				for (int k = grid.lz(); k <= grid.uz(); k++) {
					z = k * grid.getDZ();
					val00 = grid.get(i, j, k) ? 1 : 0;
					valm0 = grid.get(i, j - 1, k) ? 1 : 0;
					valmm = grid.get(i, j - 1, k - 1) ? 1 : 0;
					val0m = grid.get(i, j, k - 1) ? 1 : 0;
					sum = val00 + valm0 + valmm + val0m;
					if (sum == 1 || sum == 3) {
						line(x, y, z, x + grid.getDX(), y, z);
					}
					if (sum == 2) {
						if (val00 + valmm == 2 || val0m + valm0 == 2) {
							line(x, y, z, x + grid.getDX(), y, z);
						}
					}
				}
			}
		}
	}

	void drawXFaces(final WB_BinaryGrid3D grid) {
		int val0, valm, sum;
		double x, y, z;
		for (int i = grid.lx(); i <= grid.ux(); i++) {
			x = i * grid.getDX();
			for (int j = grid.ly(); j < grid.uy(); j++) {
				y = j * grid.getDY();
				for (int k = grid.lz(); k < grid.uz(); k++) {
					z = k * grid.getDZ();
					val0 = grid.get(i, j, k) ? 1 : 0;
					valm = grid.get(i - 1, j, k) ? 1 : 0;
					sum = val0 + valm;
					if (sum == 1) {
						home.beginShape();
						vertex(x, y, z);
						vertex(x, y + grid.getDY(), z);
						vertex(x, y + grid.getDY(), z + grid.getDZ());
						vertex(x, y, z + grid.getDZ());
						home.endShape();
					}
				}
			}
		}
	}

	void drawYEdges(final WB_BinaryGrid3D grid) {
		int val00, valm0, valmm, val0m, sum;
		double x, y, z;
		for (int j = grid.ly(); j < grid.uy(); j++) {
			y = j * grid.getDY();
			for (int i = grid.lx(); i <= grid.ux(); i++) {
				x = i * grid.getDX();
				for (int k = grid.lz(); k <= grid.uz(); k++) {
					z = k * grid.getDZ();
					val00 = grid.get(i, j, k) ? 1 : 0;
					valm0 = grid.get(i - 1, j, k) ? 1 : 0;
					valmm = grid.get(i - 1, j, k - 1) ? 1 : 0;
					val0m = grid.get(i, j, k - 1) ? 1 : 0;
					sum = val00 + valm0 + valmm + val0m;
					if (sum == 1 || sum == 3) {
						line(x, y, z, x, y + grid.getDY(), z);
					}
					if (sum == 2) {
						if (val00 + valmm == 2 || val0m + valm0 == 2) {
							line(x, y, z, x, y + grid.getDY(), z);
						}
					}
				}
			}
		}
	}

	void drawYFaces(final WB_BinaryGrid3D grid) {
		int val0, valm, sum;
		double x, y, z;
		for (int i = grid.lx(); i < grid.ux(); i++) {
			x = i * grid.getDX();
			for (int j = grid.ly(); j <= grid.uy(); j++) {
				y = j * grid.getDY();
				for (int k = grid.lz(); k < grid.uz(); k++) {
					z = k * grid.getDZ();
					val0 = grid.get(i, j, k) ? 1 : 0;
					valm = grid.get(i, j - 1, k) ? 1 : 0;
					sum = val0 + valm;
					if (sum == 1) {
						home.beginShape();
						vertex(x, y, z);
						vertex(x + grid.getDX(), y, z);
						vertex(x + grid.getDX(), y, z + grid.getDZ());
						vertex(x, y, z + grid.getDZ());
						home.endShape();
					}
				}
			}
		}
	}

	void drawZEdges(final WB_BinaryGrid3D grid) {
		int val00, valm0, valmm, val0m, sum;
		double x, y, z;
		for (int k = grid.lz(); k < grid.uz(); k++) {
			z = k * grid.getDZ();
			for (int j = grid.ly(); j <= grid.uy(); j++) {
				y = j * grid.getDY();
				for (int i = grid.lx(); i <= grid.ux(); i++) {
					x = i * grid.getDX();
					val00 = grid.get(i, j, k) ? 1 : 0;
					valm0 = grid.get(i - 1, j, k) ? 1 : 0;
					valmm = grid.get(i - 1, j - 1, k) ? 1 : 0;
					val0m = grid.get(i, j - 1, k) ? 1 : 0;
					sum = val00 + valm0 + valmm + val0m;
					if (sum == 1 || sum == 3) {
						line(x, y, z, x, y, z + grid.getDZ());
					}
					if (sum == 2) {
						if (val00 + valmm == 2 || val0m + valm0 == 2) {
							line(x, y, z, x, y, z + grid.getDZ());
						}
					}
				}
			}
		}
	}

	void drawZFaces(final WB_BinaryGrid3D grid) {
		int val0, valm, sum;
		double x, y, z;
		for (int i = grid.lx(); i < grid.ux(); i++) {
			x = i * grid.getDX();
			for (int j = grid.ly(); j < grid.uy(); j++) {
				y = j * grid.getDY();
				for (int k = grid.lz(); k <= grid.uz(); k++) {
					z = k * grid.getDZ();
					val0 = grid.get(i, j, k) ? 1 : 0;
					valm = grid.get(i, j, k - 1) ? 1 : 0;
					sum = val0 + valm;
					if (sum == 1) {
						home.beginShape();
						vertex(x, y, z);
						vertex(x + grid.getDX(), y, z);
						vertex(x + grid.getDX(), y + grid.getDY(), z);
						vertex(x, y + grid.getDY(), z);
						home.endShape();
					}
				}
			}
		}
	}

	public void drawDanzer3D(final WB_Danzer3D danzer) {
		for (final WB_DanzerTile3D tile : danzer.getTiles()) {
			drawTetrahedron(tile);
		}
	}

	public void drawUnpairedHalfedges(final HE_HalfedgeStructure mesh) {
		HE_Halfedge he;
		final Iterator<HE_Halfedge> heItr = mesh.heItr();
		home.pushStyle();
		while (heItr.hasNext()) {
			he = heItr.next();
			if (he.getPair() == null) {
				home.stroke(255, 0, 0);
				home.line(he.getVertex().xf(), he.getVertex().yf(), he.getVertex().zf(),
						he.getNextInFace().getVertex().xf(), he.getNextInFace().getVertex().yf(),
						he.getNextInFace().getVertex().zf());
			}
		}
		home.popStyle();
	}

	private boolean flipTextureV = false;

	public WB_Render3D flipTextureV(final boolean b) {
		flipTextureV = b;
		return this;
	}

	float uOffset = 0.0f;
	float uRepeat = 1.0f;
	float vOffset = 0.0f;
	float vRepeat = 1.0f;

	public WB_Render3D setUOffset(final float value) {
		uOffset = value;
		return this;
	}

	public WB_Render3D setVOffset(final float value) {
		vOffset = value;
		return this;
	}

	public WB_Render3D setURepeat(final float value) {
		uRepeat = value;
		return this;
	}

	public WB_Render3D setVRepeat(final float value) {
		vRepeat = value;
		return this;
	}

	private float getU(final HE_TextureCoordinate uvw) {
		return uOffset + uRepeat * uvw.uf();
	}

	private float getV(final HE_TextureCoordinate uvw) {
		return vOffset + vRepeat * (flipTextureV ? 1.0f - uvw.vf() : uvw.vf());
	}

	public static int getColorFromPImage(final double u, final double v, final PImage texture) {
		return texture.get(Math.max(0, Math.min((int) (u * texture.width), texture.width - 1)),
				Math.max(0, Math.min((int) (v * texture.height), texture.height - 1)));
	}
}
