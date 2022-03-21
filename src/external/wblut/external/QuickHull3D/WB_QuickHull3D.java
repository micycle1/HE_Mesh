// 
// Decompiled by Procyon v0.5.36
// 

package wblut.external.QuickHull3D;

import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Vector;

import org.eclipse.collections.impl.list.mutable.FastList;

import wblut.geom.WB_Coord;
import wblut.geom.WB_GeometryOp;
import wblut.geom.WB_Point;
import wblut.geom.WB_Vector;
import wblut.math.WB_Epsilon;

public class WB_QuickHull3D {
	protected double charLength;
	protected List<WB_Coord> points;
	protected Vertex[] pointBuffer;
	protected int[] vertexPointIndices;
	private final Face[] discardedFaces;
	private final Vertex[] maxVtxs;
	private final Vertex[] minVtxs;
	protected Vector<Face> faces;
	protected Vector<HalfEdge> horizon;
	private final FastList<Face> newFaces;
	private final VertexList unclaimed;
	private final VertexList claimed;
	private int numVertices;
	private int numPoints;
	private double tolerance;
	private static final int NONCONVEX_WRT_LARGER_FACE = 1;
	private static final int NONCONVEX = 2;

	public WB_QuickHull3D(final Collection<? extends WB_Coord> points) throws IllegalArgumentException {
		this.pointBuffer = new Vertex[0];
		this.vertexPointIndices = new int[0];
		this.discardedFaces = new Face[3];
		this.maxVtxs = new Vertex[3];
		this.minVtxs = new Vertex[3];
		this.faces = new Vector<>(16);
		this.horizon = new Vector<>(16);
		this.newFaces = new FastList();
		this.unclaimed = new VertexList();
		this.claimed = new VertexList();
		this.points = new FastList();
		for (final WB_Coord point : points) {
			this.points.add(new WB_Point(point));
		}
		this.build(false);
	}

	public WB_QuickHull3D(final Collection<? extends WB_Coord> points, final boolean triangulate) throws IllegalArgumentException {
		this.pointBuffer = new Vertex[0];
		this.vertexPointIndices = new int[0];
		this.discardedFaces = new Face[3];
		this.maxVtxs = new Vertex[3];
		this.minVtxs = new Vertex[3];
		this.faces = new Vector<>(16);
		this.horizon = new Vector<>(16);
		this.newFaces = new FastList();
		this.unclaimed = new VertexList();
		this.claimed = new VertexList();
		this.points = new FastList();
		for (final WB_Coord point : points) {
			this.points.add(new WB_Point(point));
		}
		this.build(triangulate);
	}

	public WB_QuickHull3D(final WB_Coord[] points) throws IllegalArgumentException {
		this.pointBuffer = new Vertex[0];
		this.vertexPointIndices = new int[0];
		this.discardedFaces = new Face[3];
		this.maxVtxs = new Vertex[3];
		this.minVtxs = new Vertex[3];
		this.faces = new Vector<>(16);
		this.horizon = new Vector<>(16);
		this.newFaces = new FastList();
		this.unclaimed = new VertexList();
		this.claimed = new VertexList();
		this.points = new FastList();
		for (final WB_Coord point : points) {
			this.points.add(new WB_Point(point));
		}
		this.build(false);
	}

	public WB_QuickHull3D(final WB_Coord[] points, final boolean triangulate) throws IllegalArgumentException {
		this.pointBuffer = new Vertex[0];
		this.vertexPointIndices = new int[0];
		this.discardedFaces = new Face[3];
		this.maxVtxs = new Vertex[3];
		this.minVtxs = new Vertex[3];
		this.faces = new Vector<>(16);
		this.horizon = new Vector<>(16);
		this.newFaces = new FastList();
		this.unclaimed = new VertexList();
		this.claimed = new VertexList();
		this.points = new FastList();
		for (final WB_Coord point : points) {
			this.points.add(new WB_Point(point));
		}
		this.build(triangulate);
	}

	public void build(final boolean triangulate) throws IllegalArgumentException {
		final int nump = this.points.size();
		if (nump < 4) {
			throw new IllegalArgumentException("Less than four input points specified");
		}
		this.initBuffers(nump);
		this.setPoints(this.points);
		this.buildHull();
		if (triangulate) {
			this.triangulate();
		}
	}

	protected void initBuffers(final int nump) {
		this.pointBuffer = new Vertex[nump];
		this.vertexPointIndices = new int[nump];
		this.faces.clear();
		this.claimed.clear();
		this.numPoints = nump;
	}

	private void addPointToFace(final Vertex vtx, final Face face) {
		vtx.face = face;
		if (face.outside == null) {
			this.claimed.add(vtx);
		} else {
			this.claimed.insertBefore(vtx, face.outside);
		}
		face.outside = vtx;
	}

	private void removePointFromFace(final Vertex vtx, final Face face) {
		if (vtx == face.outside) {
			if (vtx.next != null && vtx.next.face == face) {
				face.outside = vtx.next;
			} else {
				face.outside = null;
			}
		}
		this.claimed.delete(vtx);
	}

	private Vertex removeAllPointsFromFace(final Face face) {
		if (face.outside != null) {
			Vertex end;
			for (end = face.outside; end.next != null && end.next.face == face; end = end.next) {
			}
			this.claimed.delete(face.outside, end);
			end.next = null;
			return face.outside;
		}
		return null;
	}

	public void triangulate() {
		final double minArea = this.charLength * 2.220446049250313E-13;
		this.newFaces.clear();
		for (final Face face : this.faces) {
			if (face.mark == 1) {
				face.triangulate(this.newFaces, minArea);
			}
		}
		for (final Face face2 : this.newFaces) {
			this.faces.add(face2);
		}
	}

	protected void setPoints(final List<WB_Coord> points) {
		int i = 0;
		for (final WB_Coord p : points) {
			this.pointBuffer[i] = new Vertex(p, i);
			++i;
		}
	}

	protected void computeMaxAndMin() {
		final WB_Vector max = new WB_Vector();
		final WB_Vector min = new WB_Vector();
		for (int i = 0; i < 3; ++i) {
			this.maxVtxs[i] = this.pointBuffer[0];
			this.minVtxs[i] = this.pointBuffer[0];
		}
		double xm = this.pointBuffer[0].pos.xd();
		double ym = this.pointBuffer[0].pos.yd();
		double zm = this.pointBuffer[0].pos.zd();
		double xM = this.pointBuffer[0].pos.xd();
		double yM = this.pointBuffer[0].pos.yd();
		double zM = this.pointBuffer[0].pos.zd();
		for (int j = 1; j < this.numPoints; ++j) {
			final Vertex v = this.pointBuffer[j];
			final WB_Vector pnt = v.pos;
			if (pnt.xd() > max.xd()) {
				xM = pnt.xd();
				this.maxVtxs[0] = v;
			} else if (pnt.xd() < min.xd()) {
				xm = pnt.xd();
				this.minVtxs[0] = v;
			}
			if (pnt.yd() > max.yd()) {
				yM = pnt.yd();
				this.maxVtxs[1] = v;
			} else if (pnt.yd() < min.yd()) {
				ym = pnt.yd();
				this.minVtxs[1] = v;
			}
			if (pnt.zd() > max.zd()) {
				zM = pnt.zd();
				this.maxVtxs[2] = v;
			} else if (pnt.zd() < min.zd()) {
				zm = pnt.zd();
				this.minVtxs[2] = v;
			}
			max.set(xM, yM, zM);
			min.set(xm, ym, zm);
		}
		this.charLength = Math.max(max.xd() - min.xd(), max.yd() - min.yd());
		this.charLength = Math.max(max.zd() - min.zd(), this.charLength);
		this.tolerance = WB_Epsilon.EPSILON;
	}

	protected void createInitialSimplex() throws IllegalArgumentException {
		double max = 0.0;
		int imax = 0;
		for (int i = 0; i < 3; ++i) {
			final double diff = this.maxVtxs[i].pos.getd(i) - this.minVtxs[i].pos.getd(i);
			if (diff > max) {
				max = diff;
				imax = i;
			}
		}
		if (max <= this.tolerance) {
			throw new IllegalArgumentException("Input points appear to be coincident");
		}
		final Vertex[] vtx = { this.maxVtxs[imax], this.minVtxs[imax], null, null };
		final WB_Vector u01 = new WB_Vector(vtx[1].pos.xd(), vtx[1].pos.yd(), vtx[1].pos.zd());
		final WB_Vector diff2 = new WB_Vector();
		final WB_Vector nrml = new WB_Vector();
		final WB_Vector xprod = new WB_Vector();
		double maxSqr = 0.0;
		u01.subSelf(vtx[0].pos);
		u01.normalizeSelf();
		for (int j = 0; j < this.numPoints; ++j) {
			final Vertex p = this.pointBuffer[j];
			diff2.set(p.pos.xd() - vtx[0].pos.xd(), p.pos.yd() - vtx[0].pos.yd(), p.pos.zd() - vtx[0].pos.zd());
			xprod.set(u01.yd() * diff2.zd() - u01.zd() * diff2.yd(), u01.zd() * diff2.xd() - u01.xd() * diff2.zd(),
					u01.xd() * diff2.yd() - u01.yd() * diff2.xd());
			final double lenSqr = xprod.getSqLength3D();
			if (lenSqr > maxSqr && p != vtx[0] && p != vtx[1]) {
				maxSqr = lenSqr;
				vtx[2] = p;
				nrml.set(xprod);
			}
		}
		if (Math.sqrt(maxSqr) <= 100.0 * this.tolerance) {
			throw new IllegalArgumentException("Input points appear to be colinear");
		}
		nrml.normalizeSelf();
		double maxDist = 0.0;
		final double d0 = nrml.dot(vtx[2].pos);
		for (int k = 0; k < this.numPoints; ++k) {
			final Vertex p2 = this.pointBuffer[k];
			final double dist = Math.abs(nrml.dot(p2.pos) - d0);
			if (dist > maxDist && p2 != vtx[0] && p2 != vtx[1] && p2 != vtx[2]) {
				maxDist = dist;
				vtx[3] = p2;
			}
		}
		if (Math.abs(maxDist) <= 100.0 * this.tolerance) {
			throw new IllegalArgumentException("Input points appear to be coplanar");
		}
		final Face[] tris = new Face[4];
		if (nrml.dot(vtx[3].pos) - d0 < 0.0) {
			tris[0] = Face.createTriangle(vtx[0], vtx[1], vtx[2]);
			tris[1] = Face.createTriangle(vtx[3], vtx[1], vtx[0]);
			tris[2] = Face.createTriangle(vtx[3], vtx[2], vtx[1]);
			tris[3] = Face.createTriangle(vtx[3], vtx[0], vtx[2]);
			for (int l = 0; l < 3; ++l) {
				final int m = (l + 1) % 3;
				tris[l + 1].getEdge(1).setOpposite(tris[m + 1].getEdge(0));
				tris[l + 1].getEdge(2).setOpposite(tris[0].getEdge(m));
			}
		} else {
			tris[0] = Face.createTriangle(vtx[0], vtx[2], vtx[1]);
			tris[1] = Face.createTriangle(vtx[3], vtx[0], vtx[1]);
			tris[2] = Face.createTriangle(vtx[3], vtx[1], vtx[2]);
			tris[3] = Face.createTriangle(vtx[3], vtx[2], vtx[0]);
			for (int l = 0; l < 3; ++l) {
				final int m = (l + 1) % 3;
				tris[l + 1].getEdge(0).setOpposite(tris[m + 1].getEdge(1));
				tris[l + 1].getEdge(2).setOpposite(tris[0].getEdge((3 - l) % 3));
			}
		}
		for (int l = 0; l < 4; ++l) {
			this.faces.add(tris[l]);
		}
		for (int l = 0; l < this.numPoints; ++l) {
			final Vertex v = this.pointBuffer[l];
			if (v != vtx[0] && v != vtx[1] && v != vtx[2]) {
				if (v != vtx[3]) {
					maxDist = this.tolerance;
					Face maxFace = null;
					for (int k2 = 0; k2 < 4; ++k2) {
						final double dist2 = tris[k2].distanceToPlane(v.pos);
						if (dist2 > maxDist) {
							maxFace = tris[k2];
							maxDist = dist2;
						}
					}
					if (maxFace != null) {
						this.addPointToFace(v, maxFace);
					}
				}
			}
		}
	}

	public int getNumVertices() {
		return this.numVertices;
	}

	public List<WB_Coord> getVertices() {
		final List<WB_Coord> vtxs = new FastList();
		for (int i = 0; i < this.numVertices; ++i) {
			vtxs.add(this.points.get(this.vertexPointIndices[i]));
		}
		return vtxs;
	}

	public int[] getVertexPointIndices() {
		final int[] indices = new int[this.numVertices];
		for (int i = 0; i < this.numVertices; ++i) {
			indices[i] = this.vertexPointIndices[i];
		}
		return indices;
	}

	public int getNumFaces() {
		return this.faces.size();
	}

	public int[][] getFaces() {
		final int[][] allFaces = new int[this.faces.size()][];
		int k = 0;
		for (final Face face : this.faces) {
			this.getFaceIndices(allFaces[k] = new int[face.numVertices()], face);
			++k;
		}
		return allFaces;
	}

	private void getFaceIndices(final int[] indices, final Face face) {
		HalfEdge hedge = face.he0;
		int k = 0;
		do {
			final int idx = hedge.head().hullindex;
			indices[k++] = idx;
			hedge = hedge.next;
		} while (hedge != face.he0);
	}

	protected void resolveUnclaimedPoints(final List<Face> newFaces) {
		Vertex vtx;
		for (Vertex vtxNext = vtx = this.unclaimed.first(); vtx != null; vtx = vtxNext) {
			vtxNext = vtx.next;
			double maxDist = this.tolerance;
			Face maxFace = null;
			for (final Face newFace : newFaces) {
				if (newFace.mark == 1) {
					final double dist = newFace.distanceToPlane(vtx);
					if (dist > maxDist) {
						maxDist = dist;
						maxFace = newFace;
					}
					if (maxDist > 1000.0 * this.tolerance) {
						break;
					}
					continue;
				}
			}
			if (maxFace != null) {
				this.addPointToFace(vtx, maxFace);
			}
		}
	}

	protected void deleteFacePoints(final Face face, final Face absorbingFace) {
		final Vertex faceVtxs = this.removeAllPointsFromFace(face);
		if (faceVtxs != null) {
			if (absorbingFace == null) {
				this.unclaimed.addAll(faceVtxs);
			} else {
				Vertex vtx;
				for (Vertex vtxNext = vtx = faceVtxs; vtx != null; vtx = vtxNext) {
					vtxNext = vtx.next;
					final double dist = absorbingFace.distanceToPlane(vtx);
					if (dist > this.tolerance) {
						this.addPointToFace(vtx, absorbingFace);
					} else {
						this.unclaimed.add(vtx);
					}
				}
			}
		}
	}

	protected double oppFaceDistance(final HalfEdge he) {
		return he.face.distanceToPlane(he.opposite.face.getCentroid());
	}

	private boolean doAdjacentMerge(final Face face, final int mergeType) {
		HalfEdge hedge = face.he0;
		boolean convex = true;
		do {
			final Face oppFace = hedge.oppositeFace();
			boolean merge = false;
			if (mergeType == 2) {
				if (this.oppFaceDistance(hedge) > -this.tolerance || this.oppFaceDistance(hedge.opposite) > -this.tolerance) {
					merge = true;
				}
			} else if (face.area > oppFace.area) {
				if (this.oppFaceDistance(hedge) > -this.tolerance) {
					merge = true;
				} else if (this.oppFaceDistance(hedge.opposite) > -this.tolerance) {
					convex = false;
				}
			} else if (this.oppFaceDistance(hedge.opposite) > -this.tolerance) {
				merge = true;
			} else if (this.oppFaceDistance(hedge) > -this.tolerance) {
				convex = false;
			}
			if (merge) {
				for (int numd = face.mergeAdjacentFace(hedge, this.discardedFaces), i = 0; i < numd; ++i) {
					this.deleteFacePoints(this.discardedFaces[i], face);
				}
				return true;
			}
			hedge = hedge.next;
		} while (hedge != face.he0);
		if (!convex) {
			face.mark = 2;
		}
		return false;
	}

	protected void calculateHorizon(final double eyePntx, final double eyePnty, final double eyePntz, HalfEdge edge0, final Face face,
			final Vector<HalfEdge> horizon) {
		this.deleteFacePoints(face, null);
		face.mark = 3;
		HalfEdge edge;
		if (edge0 == null) {
			edge0 = (edge = face.getEdge(0));
		} else {
			edge = edge0.getNext();
		}
		do {
			final Face oppFace = edge.oppositeFace();
			if (oppFace.mark == 1) {
				if (oppFace.distanceToPlane(eyePntx, eyePnty, eyePntz) > this.tolerance) {
					this.calculateHorizon(eyePntx, eyePnty, eyePntz, edge.getOpposite(), oppFace, horizon);
				} else {
					horizon.add(edge);
				}
			}
			edge = edge.getNext();
		} while (edge != edge0);
	}

	private HalfEdge addAdjoiningFace(final Vertex eyeVtx, final HalfEdge he) {
		final Face face = Face.createTriangle(eyeVtx, he.tail(), he.head());
		this.faces.add(face);
		face.getEdge(-1).setOpposite(he.getOpposite());
		return face.getEdge(0);
	}

	protected void addNewFaces(final List<Face> newFaces, final Vertex eyeVtx, final Vector<HalfEdge> horizon) {
		newFaces.clear();
		HalfEdge hedgeSidePrev = null;
		HalfEdge hedgeSideBegin = null;
		for (final HalfEdge horizonHe : horizon) {
			final HalfEdge hedgeSide = this.addAdjoiningFace(eyeVtx, horizonHe);
			if (hedgeSidePrev != null) {
				hedgeSide.next.setOpposite(hedgeSidePrev);
			} else {
				hedgeSideBegin = hedgeSide;
			}
			newFaces.add(hedgeSide.getFace());
			hedgeSidePrev = hedgeSide;
		}
		hedgeSideBegin.next.setOpposite(hedgeSidePrev);
	}

	protected Vertex nextPointToAdd() {
		if (!this.claimed.isEmpty()) {
			final Face eyeFace = this.claimed.first().face;
			Vertex eyeVtx = null;
			double maxDist = 0.0;
			for (Vertex vtx = eyeFace.outside; vtx != null && vtx.face == eyeFace; vtx = vtx.next) {
				final double dist = eyeFace.distanceToPlane(vtx);
				if (dist > maxDist) {
					maxDist = dist;
					eyeVtx = vtx;
				}
			}
			return eyeVtx;
		}
		return null;
	}

	protected void addPointToHull(final Vertex eyeVtx) {
		this.horizon.clear();
		this.unclaimed.clear();
		this.removePointFromFace(eyeVtx, eyeVtx.face);
		this.calculateHorizon(eyeVtx.pos.xd(), eyeVtx.pos.yd(), eyeVtx.pos.zd(), null, eyeVtx.face, this.horizon);
		this.newFaces.clear();
		this.addNewFaces(this.newFaces, eyeVtx, this.horizon);
		for (final Face face : this.newFaces) {
			if (face.mark == 1) {
				while (this.doAdjacentMerge(face, 1)) {
				}
			}
		}
		for (final Face face : this.newFaces) {
			if (face.mark == 2) {
				face.mark = 1;
				while (this.doAdjacentMerge(face, 2)) {
				}
			}
		}
		this.resolveUnclaimedPoints(this.newFaces);
	}

	protected void buildHull() {
		this.computeMaxAndMin();
		this.createInitialSimplex();
		Vertex eyeVtx;
		while ((eyeVtx = this.nextPointToAdd()) != null) {
			this.addPointToHull(eyeVtx);
		}
		this.reindexFacesAndVertices();
	}

	private void markFaceVertices(final Face face) {
		HalfEdge he2;
		final HalfEdge he0 = he2 = face.getFirstEdge();
		do {
			he2.head().flag = true;
			he2.head().hullindex = 0;
			he2 = he2.next;
		} while (he2 != he0);
	}

	protected void reindexFacesAndVertices() {
		for (int i = 1; i < this.numPoints; ++i) {
			final Vertex p = this.pointBuffer[i];
			p.flag = false;
			p.hullindex = -1;
		}
		final Iterator<Face> it = this.faces.iterator();
		while (it.hasNext()) {
			final Face face = it.next();
			if (face.mark != 1) {
				it.remove();
			} else {
				this.markFaceVertices(face);
			}
		}
		this.numVertices = 0;
		for (int i = 0; i < this.numPoints; ++i) {
			final Vertex vtx = this.pointBuffer[i];
			if (vtx.flag) {
				this.vertexPointIndices[this.numVertices] = i;
				vtx.hullindex = this.numVertices++;
			}
		}
	}

	static class Vertex {
		WB_Vector pos;
		int hullindex;
		Vertex prev;
		Vertex next;
		Face face;
		boolean flag;

		protected Vertex() {
			this.flag = false;
		}

		protected Vertex(final WB_Coord p, final int idx) {
			this.pos = new WB_Vector(p);
			this.hullindex = idx;
			this.flag = false;
		}
	}

	static class Face {
		HalfEdge he0;
		private final WB_Vector normal;
		double area;
		private WB_Vector centroid;
		double planeOffset;
		int index;
		int numVerts;
		Face next;
		static final int VISIBLE = 1;
		static final int NON_CONVEX = 2;
		static final int DELETED = 3;
		int mark;
		Vertex outside;

		protected WB_Vector computeCentroid() {
			this.centroid = new WB_Vector();
			HalfEdge he = this.he0;
			do {
				final Vertex v = he.head();
				this.centroid.addSelf(v.pos);
				he = he.next;
			} while (he != this.he0);
			return this.centroid.divSelf(this.numVerts);
		}

		protected void computeNormal(final WB_Vector normal, final double minArea) {
			this.computeNormal(normal);
			if (this.area < minArea) {
				System.out.println("area=" + this.area);
				HalfEdge hedgeMax = null;
				double lenSqrMax = 0.0;
				HalfEdge hedge = this.he0;
				do {
					final double lenSqr = hedge.lengthSquared();
					if (lenSqr > lenSqrMax) {
						hedgeMax = hedge;
						lenSqrMax = lenSqr;
					}
					hedge = hedge.next;
				} while (hedge != this.he0);
				final Vertex p2 = hedgeMax.head();
				final Vertex p3 = hedgeMax.tail();
				final double lenMax = Math.sqrt(lenSqrMax);
				final WB_Vector u = p2.pos.sub(p3.pos).divSelf(lenMax);
				final double dot = normal.dot(u);
				normal.addMul(-dot, u);
				normal.normalizeSelf();
			}
		}

		protected void computeNormal(final WB_Vector normal) {
			HalfEdge he1 = this.he0.next;
			HalfEdge he2 = he1.next;
			final Vertex p0 = this.he0.head();
			Vertex p2 = he1.head();
			double d2x = p2.pos.xd() - p0.pos.xd();
			double d2y = p2.pos.yd() - p0.pos.yd();
			double d2z = p2.pos.zd() - p0.pos.zd();
			normal.set(0.0, 0.0, 0.0);
			this.numVerts = 2;
			while (he2 != this.he0) {
				final double d1x = d2x;
				final double d1y = d2y;
				final double d1z = d2z;
				p2 = he2.head();
				d2x = p2.pos.xd() - p0.pos.xd();
				d2y = p2.pos.yd() - p0.pos.yd();
				d2z = p2.pos.zd() - p0.pos.zd();
				normal.addSelf(d1y * d2z - d1z * d2y, d1z * d2x - d1x * d2z, d1x * d2y - d1y * d2x);
				he1 = he2;
				he2 = he2.next;
				++this.numVerts;
			}
			this.area = normal.getLength3D();
			normal.normalizeSelf();
		}

		private void computeNormalAndCentroid() {
			this.computeNormal(this.normal);
			this.centroid = this.computeCentroid();
			this.planeOffset = this.normal.dot(this.centroid);
			int numv = 0;
			HalfEdge he = this.he0;
			do {
				++numv;
				he = he.next;
			} while (he != this.he0);
			if (numv != this.numVerts) {
				throw new InternalErrorException("face " + this.getVertexString() + " numVerts=" + this.numVerts + " should be " + numv);
			}
		}

		private void computeNormalAndCentroid(final double minArea) {
			this.computeNormal(this.normal, minArea);
			this.centroid = this.computeCentroid();
			this.planeOffset = this.normal.dot(this.centroid);
		}

		protected static Face createTriangle(final Vertex v0, final Vertex v1, final Vertex v2) {
			return createTriangle(v0, v1, v2, 0.0);
		}

		protected static Face createTriangle(final Vertex v0, final Vertex v1, final Vertex v2, final double minArea) {
			final Face face = new Face();
			final HalfEdge he0 = new HalfEdge(v0, face);
			final HalfEdge he2 = new HalfEdge(v1, face);
			final HalfEdge he3 = new HalfEdge(v2, face);
			he0.prev = he3;
			he0.next = he2;
			he2.prev = he0;
			he2.next = he3;
			he3.prev = he2;
			he3.next = he0;
			face.he0 = he0;
			face.computeNormalAndCentroid(minArea);
			return face;
		}

		protected static Face create(final FastList<Vertex> vtxArray, final int[] indices) {
			final Face face = new Face();
			HalfEdge hePrev = null;
			for (final int indice : indices) {
				final HalfEdge he = new HalfEdge(vtxArray.get(indice), face);
				if (hePrev != null) {
					he.setPrev(hePrev);
					hePrev.setNext(he);
				} else {
					face.he0 = he;
				}
				hePrev = he;
			}
			face.he0.setPrev(hePrev);
			hePrev.setNext(face.he0);
			face.computeNormalAndCentroid();
			return face;
		}

		protected Face() {
			this.mark = 1;
			this.normal = new WB_Vector();
			this.centroid = new WB_Vector();
			this.mark = 1;
		}

		protected HalfEdge getEdge(int i) {
			HalfEdge he = this.he0;
			while (i > 0) {
				he = he.next;
				--i;
			}
			while (i < 0) {
				he = he.prev;
				++i;
			}
			return he;
		}

		protected HalfEdge getFirstEdge() {
			return this.he0;
		}

		protected HalfEdge findEdge(final Vertex vt, final Vertex vh) {
			HalfEdge he = this.he0;
			while (he.head() != vh || he.tail() != vt) {
				he = he.next;
				if (he == this.he0) {
					return null;
				}
			}
			return he;
		}

		protected double distanceToPlane(final WB_Coord p) {
			return this.normal.xd() * p.xd() + this.normal.yd() * p.yd() + this.normal.zd() * p.zd() - this.planeOffset;
		}

		protected double distanceToPlane(final double x, final double y, final double z) {
			return this.normal.xd() * x + this.normal.yd() * y + this.normal.zd() * z - this.planeOffset;
		}

		protected double distanceToPlane(final double[] p) {
			return this.normal.xd() * p[0] + this.normal.yd() * p[1] + this.normal.zd() * p[2] - this.planeOffset;
		}

		protected double distanceToPlane(final Vertex v) {
			return this.normal.xd() * v.pos.xd() + this.normal.yd() * v.pos.yd() + this.normal.zd() * v.pos.zd() - this.planeOffset;
		}

		protected WB_Vector getNormal() {
			return this.normal;
		}

		protected WB_Vector getCentroid() {
			return this.centroid;
		}

		protected int numVertices() {
			return this.numVerts;
		}

		protected String getVertexString() {
			String s = null;
			HalfEdge he = this.he0;
			do {
				if (s == null) {
					s = new StringBuilder().append(he.head().hullindex).toString();
				} else {
					s = String.valueOf(s) + " " + he.head().hullindex;
				}
				he = he.next;
			} while (he != this.he0);
			return s;
		}

		protected void getVertexIndices(final int[] idxs) {
			HalfEdge he = this.he0;
			int i = 0;
			do {
				idxs[i++] = he.head().hullindex;
				he = he.next;
			} while (he != this.he0);
		}

		private Face connectHalfEdges(final HalfEdge hedgePrev, final HalfEdge hedge) {
			Face discardedFace = null;
			if (hedgePrev.oppositeFace() == hedge.oppositeFace()) {
				final Face oppFace = hedge.oppositeFace();
				if (hedgePrev == this.he0) {
					this.he0 = hedge;
				}
				HalfEdge hedgeOpp;
				if (oppFace.numVertices() == 3) {
					hedgeOpp = hedge.getOpposite().prev.getOpposite();
					oppFace.mark = 3;
					discardedFace = oppFace;
				} else {
					hedgeOpp = hedge.getOpposite().next;
					if (oppFace.he0 == hedgeOpp.prev) {
						oppFace.he0 = hedgeOpp;
					}
					hedgeOpp.prev = hedgeOpp.prev.prev;
					hedgeOpp.prev.next = hedgeOpp;
				}
				hedge.prev = hedgePrev.prev;
				hedge.prev.next = hedge;
				hedge.opposite = hedgeOpp;
				hedgeOpp.opposite = hedge;
				oppFace.computeNormalAndCentroid();
			} else {
				hedgePrev.next = hedge;
				hedge.prev = hedgePrev;
			}
			return discardedFace;
		}

		void checkConsistency() {
			HalfEdge hedge = this.he0;
			double maxd = 0.0;
			int numv = 0;
			if (this.numVerts < 3) {
				throw new InternalErrorException("degenerate face: " + this.getVertexString());
			}
			do {
				final HalfEdge hedgeOpp = hedge.getOpposite();
				if (hedgeOpp == null) {
					throw new InternalErrorException(
							"face " + this.getVertexString() + ": " + "unreflected half edge " + hedge.getVertexString());
				}
				if (hedgeOpp.getOpposite() != hedge) {
					throw new InternalErrorException("face " + this.getVertexString() + ": " + "opposite half edge "
							+ hedgeOpp.getVertexString() + " has opposite " + hedgeOpp.getOpposite().getVertexString());
				}
				if (hedgeOpp.head() != hedge.tail() || hedge.head() != hedgeOpp.tail()) {
					throw new InternalErrorException("face " + this.getVertexString() + ": " + "half edge " + hedge.getVertexString()
							+ " reflected by " + hedgeOpp.getVertexString());
				}
				final Face oppFace = hedgeOpp.face;
				if (oppFace == null) {
					throw new InternalErrorException(
							"face " + this.getVertexString() + ": " + "no face on half edge " + hedgeOpp.getVertexString());
				}
				if (oppFace.mark == 3) {
					throw new InternalErrorException(
							"face " + this.getVertexString() + ": " + "opposite face " + oppFace.getVertexString() + " not on hull");
				}
				final double d = Math.abs(this.distanceToPlane(hedge.head()));
				if (d > maxd) {
					maxd = d;
				}
				++numv;
				hedge = hedge.next;
			} while (hedge != this.he0);
			if (numv != this.numVerts) {
				throw new InternalErrorException("face " + this.getVertexString() + " numVerts=" + this.numVerts + " should be " + numv);
			}
		}

		protected int mergeAdjacentFace(final HalfEdge hedgeAdj, final Face[] discarded) {
			final Face oppFace = hedgeAdj.oppositeFace();
			int numDiscarded = 0;
			discarded[numDiscarded++] = oppFace;
			oppFace.mark = 3;
			final HalfEdge hedgeOpp = hedgeAdj.getOpposite();
			HalfEdge hedgeAdjPrev = hedgeAdj.prev;
			HalfEdge hedgeAdjNext = hedgeAdj.next;
			HalfEdge hedgeOppPrev = hedgeOpp.prev;
			HalfEdge hedgeOppNext;
			for (hedgeOppNext = hedgeOpp.next; hedgeAdjPrev
					.oppositeFace() == oppFace; hedgeAdjPrev = hedgeAdjPrev.prev, hedgeOppNext = hedgeOppNext.next) {
			}
			while (hedgeAdjNext.oppositeFace() == oppFace) {
				hedgeOppPrev = hedgeOppPrev.prev;
				hedgeAdjNext = hedgeAdjNext.next;
			}
			for (HalfEdge hedge = hedgeOppNext; hedge != hedgeOppPrev.next; hedge = hedge.next) {
				hedge.face = this;
			}
			if (hedgeAdj == this.he0) {
				this.he0 = hedgeAdjNext;
			}
			Face discardedFace = this.connectHalfEdges(hedgeOppPrev, hedgeAdjNext);
			if (discardedFace != null) {
				discarded[numDiscarded++] = discardedFace;
			}
			discardedFace = this.connectHalfEdges(hedgeAdjPrev, hedgeOppNext);
			if (discardedFace != null) {
				discarded[numDiscarded++] = discardedFace;
			}
			this.computeNormalAndCentroid();
			this.checkConsistency();
			return numDiscarded;
		}

		protected void triangulate(final List<Face> newFaces, final double minArea) {
			if (this.numVertices() < 4) {
				return;
			}
			final Vertex v0 = this.he0.head();
			HalfEdge hedge = this.he0.next;
			HalfEdge oppPrev = hedge.opposite;
			Face face0 = null;
			for (hedge = hedge.next; hedge != this.he0.prev; hedge = hedge.next) {
				final Face face2 = createTriangle(v0, hedge.prev.head(), hedge.head(), minArea);
				face2.he0.next.setOpposite(oppPrev);
				face2.he0.prev.setOpposite(hedge.opposite);
				oppPrev = face2.he0;
				newFaces.add(face2);
				if (face0 == null) {
					face0 = face2;
				}
			}
			hedge = new HalfEdge(this.he0.prev.prev.head(), this);
			hedge.setOpposite(oppPrev);
			hedge.prev = this.he0;
			hedge.prev.next = hedge;
			hedge.next = this.he0.prev;
			hedge.next.prev = hedge;
			this.computeNormalAndCentroid(minArea);
			this.checkConsistency();
			for (Face face2 = face0; face2 != null; face2 = face2.next) {
				face2.checkConsistency();
			}
		}
	}

	static class HalfEdge {
		Vertex vertex;
		Face face;
		HalfEdge next;
		HalfEdge prev;
		HalfEdge opposite;

		protected HalfEdge(final Vertex v, final Face f) {
			this.vertex = v;
			this.face = f;
		}

		protected HalfEdge() {
		}

		protected void setNext(final HalfEdge edge) {
			this.next = edge;
		}

		protected HalfEdge getNext() {
			return this.next;
		}

		protected void setPrev(final HalfEdge edge) {
			this.prev = edge;
		}

		protected HalfEdge getPrev() {
			return this.prev;
		}

		protected Face getFace() {
			return this.face;
		}

		protected HalfEdge getOpposite() {
			return this.opposite;
		}

		protected void setOpposite(final HalfEdge edge) {
			this.opposite = edge;
			edge.opposite = this;
		}

		protected Vertex head() {
			return this.vertex;
		}

		protected Vertex tail() {
			return (this.prev != null) ? this.prev.vertex : null;
		}

		protected Face oppositeFace() {
			return (this.opposite != null) ? this.opposite.face : null;
		}

		protected String getVertexString() {
			if (this.tail() != null) {
				return this.tail().hullindex + "-" + this.head().hullindex;
			}
			return "?-" + this.head().hullindex;
		}

		protected double length() {
			if (this.tail() != null) {
				return WB_GeometryOp.getDistance3D(this.head().pos, this.tail().pos);
			}
			return -1.0;
		}

		protected double lengthSquared() {
			if (this.tail() != null) {
				return WB_GeometryOp.getSqDistance3D(this.head().pos, this.tail().pos);
			}
			return -1.0;
		}
	}

	static class VertexList {
		private Vertex head;
		private Vertex tail;

		protected void clear() {
			final Vertex vertex = null;
			this.tail = vertex;
			this.head = vertex;
		}

		protected void add(final Vertex vtx) {
			if (this.head == null) {
				this.head = vtx;
			} else {
				this.tail.next = vtx;
			}
			vtx.prev = this.tail;
			vtx.next = null;
			this.tail = vtx;
		}

		protected void addAll(Vertex vtx) {
			if (this.head == null) {
				this.head = vtx;
			} else {
				this.tail.next = vtx;
			}
			vtx.prev = this.tail;
			while (vtx.next != null) {
				vtx = vtx.next;
			}
			this.tail = vtx;
		}

		protected void delete(final Vertex vtx) {
			if (vtx.prev == null) {
				this.head = vtx.next;
			} else {
				vtx.prev.next = vtx.next;
			}
			if (vtx.next == null) {
				this.tail = vtx.prev;
			} else {
				vtx.next.prev = vtx.prev;
			}
		}

		protected void delete(final Vertex vtx1, final Vertex vtx2) {
			if (vtx1.prev == null) {
				this.head = vtx2.next;
			} else {
				vtx1.prev.next = vtx2.next;
			}
			if (vtx2.next == null) {
				this.tail = vtx1.prev;
			} else {
				vtx2.next.prev = vtx1.prev;
			}
		}

		protected void insertBefore(final Vertex vtx, final Vertex next) {
			vtx.prev = next.prev;
			if (next.prev == null) {
				this.head = vtx;
			} else {
				next.prev.next = vtx;
			}
			vtx.next = next;
			next.prev = vtx;
		}

		protected Vertex first() {
			return this.head;
		}

		protected boolean isEmpty() {
			return this.head == null;
		}
	}

	protected static class InternalErrorException extends RuntimeException {
		private static final long serialVersionUID = 4257980575065868777L;

		protected InternalErrorException(final String msg) {
			super(msg);
		}
	}
}
