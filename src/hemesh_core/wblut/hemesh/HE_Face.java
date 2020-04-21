package wblut.hemesh;

import java.util.List;

import wblut.geom.WB_AABB;
import wblut.geom.WB_Classification;
import wblut.geom.WB_Coord;
import wblut.geom.WB_CoordCollection;
import wblut.geom.WB_CoordinateSystem;
import wblut.geom.WB_GeometryOp3D;
import wblut.geom.WB_JTS;
import wblut.geom.WB_List;
import wblut.geom.WB_Plane;
import wblut.geom.WB_Point;
import wblut.geom.WB_Polygon;
import wblut.geom.WB_Triangle;
import wblut.geom.WB_TriangleFactory;
import wblut.geom.WB_Vector;
import wblut.math.WB_Epsilon;

public class HE_Face extends HE_MeshElement implements Comparable<HE_Face>, WB_TriangleFactory {
	private HE_Halfedge _halfedge;
	private int textureId;
	private int[] triangles;

	public HE_Face() {
		super();
		triangles = null;
	}

	public HE_Halfedge getHalfedge() {
		return _halfedge;
	}

	protected void setHalfedge(final HE_Halfedge halfedge) {
		_halfedge = halfedge;
	}

	protected void clearHalfedge() {
		_halfedge = null;
	}

	public HE_Halfedge getHalfedge(final HE_Vertex v) {
		HE_Halfedge he = _halfedge;
		if (he == null) {
			return null;
		}
		do {
			if (he.getVertex() == v) {
				return he;
			}
			he = he.getNextInFace();
		} while (he != _halfedge);
		return null;
	}

	public HE_Halfedge getHalfedge(final HE_Face f) {
		if (getHalfedge() == null || f == null) {
			return null;
		}
		HE_Halfedge he = getHalfedge();
		do {
			if (he.getPair() != null && he.getPair().getFace() != null && he.getPair().getFace() == f) {
				return he;
			}
			he = he.getNextInFace();
		} while (he != getHalfedge());
		return null;
	}

	public WB_Point getFaceCenter() {
		return HE_MeshOp.getFaceCenter(this);
	}

	public WB_Vector getFaceNormal() {
		return HE_MeshOp.getFaceNormal(this);
	}

	public double getFaceArea() {
		return HE_MeshOp.getFaceArea(this);
	}

	public int getFaceDegree() {
		int result = 0;
		if (_halfedge == null) {
			return 0;
		}
		HE_Halfedge he = _halfedge;
		do {
			result++;
			he = he.getNextInFace();
		} while (he != _halfedge);
		return result;
	}

	public HE_FaceEdgeCirculator feCrc() {
		return new HE_FaceEdgeCirculator(this);
	}

	public HE_FaceFaceCirculator ffCrc() {
		return new HE_FaceFaceCirculator(this);
	}

	public HE_FaceVertexCirculator fvCrc() {
		return new HE_FaceVertexCirculator(this);
	}

	public HE_FaceHalfedgeInnerCirculator fheiCrc() {
		return new HE_FaceHalfedgeInnerCirculator(this);
	}

	public HE_FaceHalfedgeOuterCirculator fheoCrc() {
		return new HE_FaceHalfedgeOuterCirculator(this);
	}

	public HE_FaceEdgeRevCirculator feRevCrc() {
		return new HE_FaceEdgeRevCirculator(this);
	}

	public HE_FaceFaceRevCirculator ffRevCrc() {
		return new HE_FaceFaceRevCirculator(this);
	}

	public HE_FaceHalfedgeInnerRevCirculator fheiRevCrc() {
		return new HE_FaceHalfedgeInnerRevCirculator(this);
	}

	public HE_FaceHalfedgeOuterRevCirculator fheoRevCrc() {
		return new HE_FaceHalfedgeOuterRevCirculator(this);
	}

	public HE_FaceVertexRevCirculator fvRevCrc() {
		return new HE_FaceVertexRevCirculator(this);
	}

	public List<HE_Vertex> getFaceVertices() {
		final HE_VertexList fv = new HE_VertexList();
		if (_halfedge == null) {
			return fv;
		}
		HE_Halfedge he = _halfedge;
		do {
			fv.add(he.getVertex());
			he = he.getNextInFace();
		} while (he != _halfedge);
		return fv.asUnmodifiable();
	}

	public List<HE_Vertex> getUniqueFaceVertices() {
		final HE_VertexList fv = new HE_VertexList();
		if (_halfedge == null) {
			return fv;
		}
		HE_Halfedge he = _halfedge;
		do {
			if (!fv.contains(he.getVertex())) {
				fv.add(he.getVertex());
			}
			he = he.getNextInFace();
		} while (he != _halfedge);
		return fv.asUnmodifiable();
	}

	public List<HE_Halfedge> getFaceHalfedges() {
		final HE_HalfedgeList fhe = new HE_HalfedgeList();
		if (_halfedge == null) {
			return fhe;
		}
		HE_Halfedge he = _halfedge;
		do {
			if (!fhe.contains(he)) {
				fhe.add(he);
			}
			he = he.getNextInFace();
		} while (he != _halfedge);
		return fhe.asUnmodifiable();
	}

	public List<HE_Halfedge> getFaceHalfedgesTwoSided() {
		final HE_HalfedgeList fhe = new HE_HalfedgeList();
		if (_halfedge == null) {
			return fhe;
		}
		HE_Halfedge he = _halfedge;
		do {
			if (!fhe.contains(he)) {
				fhe.add(he);
				if (he.getPair() != null) {
					if (!fhe.contains(he.getPair())) {
						fhe.add(he.getPair());
					}
				}
			}
			he = he.getNextInFace();
		} while (he != _halfedge);
		return fhe.asUnmodifiable();
	}

	public List<HE_Halfedge> getFaceEdges() {
		final HE_HalfedgeList fe = new HE_HalfedgeList();
		if (_halfedge == null) {
			return fe;
		}
		HE_Halfedge he = _halfedge;
		do {
			if (he.isEdge()) {
				if (!fe.contains(he)) {
					fe.add(he);
				}
			} else {
				if (!fe.contains(he.getPair())) {
					fe.add(he.getPair());
				}
			}
			he = he.getNextInFace();
		} while (he != _halfedge);
		return fe.asUnmodifiable();
	}

	public List<HE_Face> getNeighborFaces() {
		final HE_FaceList ff = new HE_FaceList();
		if (getHalfedge() == null) {
			return ff;
		}
		HE_Halfedge he = getHalfedge();
		do {
			final HE_Halfedge hep = he.getPair();
			if (hep != null && hep.getFace() != null) {
				if (hep.getFace() != this) {
					if (!ff.contains(hep.getFace())) {
						ff.add(hep.getFace());
					}
				}
			}
			he = he.getNextInFace();
		} while (he != getHalfedge());
		return ff.asUnmodifiable();
	}

	public List<HE_TextureCoordinate> getFaceUVWs() {
		final WB_List<HE_TextureCoordinate> fv = new WB_List<>();
		if (_halfedge == null) {
			return fv;
		}
		HE_Halfedge he = _halfedge;
		do {
			fv.add(he.getUVW());
			he = he.getNextInFace();
		} while (he != _halfedge);
		return fv.asUnmodifiable();
	}

	public void move(final WB_Coord c) {
		HE_Halfedge he = _halfedge;
		do {
			he.getVertex().getPosition().addSelf(c);
			he = he.getNextInFace();
		} while (he != _halfedge);
	}

	@Override
	public int compareTo(final HE_Face f) {
		if (f.getHalfedge() == null) {
			if (getHalfedge() == null) {
				return 0;
			} else {
				return 1;
			}
		} else if (getHalfedge() == null) {
			return -1;
		}
		return getHalfedge().compareTo(f.getHalfedge());
	}

	@Override
	public int[] getTriangles() {
		return getTriangles(true);
	}

	public int[] getTriangles(final boolean optimize) {
		if (triangles != null) {
			return triangles;
		}
		final int fo = getFaceDegree();
		if (fo < 3) {
			return triangles = new int[] { 0, 0, 0 };
		} else if (fo == 3) {
			return triangles = new int[] { 0, 1, 2 };
		} else if (isDegenerate()) {
			triangles = new int[3 * (fo - 2)];
			for (int i = 0; i < fo - 2; i++) {
				triangles[3 * i] = 0;
				triangles[3 * i + 1] = i + 1;
				triangles[3 * i + 2] = i + 2;
			}
			return triangles;
		} else if (fo == 4) {
			final WB_Point[] points = new WB_Point[4];
			int i = 0;
			HE_Halfedge he = _halfedge;
			do {
				points[i] = new WB_Point(he.getVertex().xd(), he.getVertex().yd(), he.getVertex().zd());
				he = he.getNextInFace();
				i++;
			} while (he != _halfedge);
			return triangles = WB_JTS.PolygonTriangulatorJTS.triangulateQuad(points[0], points[1], points[2],
					points[3]);
		}
		return triangles = new WB_JTS.PolygonTriangulatorJTS()
				.triangulatePolygon2D(HE_MeshOp.getOrthoPolygon(this), optimize).getTriangles();
	}

	public WB_AABB getAABB() {
		return HE_MeshOp.getAABB(this);
	}

	@Override
	public String toString() {
		String s = "HE_Face key: " + getKey() + ". Connects " + getFaceDegree() + " vertices: ";
		HE_Halfedge he = getHalfedge();
		for (int i = 0; i < getFaceDegree() - 1; i++) {
			s += he.getVertex().getKey() + "-";
			he = he.getNextInFace();
		}
		s += he.getVertex().getKey() + "." + " (" + getLabel() + "," + getInternalLabel() + ")";
		return s;
	}

	public boolean isPlanar() {
		final WB_Plane P = HE_MeshOp.getPlane(this);
		HE_Halfedge he = getHalfedge();
		do {
			if (!WB_Epsilon.isZero(WB_GeometryOp3D.getDistance3D(he.getVertex(), P))) {
				return false;
			}
			he = he.getNextInFace();
		} while (he != getHalfedge());
		return true;
	}

	public boolean isBoundary() {
		HE_Halfedge he = _halfedge;
		do {
			if (he.getPair().getFace() == null) {
				return true;
			}
			he = he.getNextInFace();
		} while (he != _halfedge);
		return false;
	}

	public boolean isDegenerate() {
		return WB_Vector.getLength3D(HE_MeshOp.getFaceNormal(this)) < 0.5;
	}

	public void copyProperties(final HE_Face el) {
		super.copyProperties(el);
		textureId = el.textureId;
	}

	@Override
	public void clear() {
		_halfedge = null;
	}

	public int getTextureId() {
		return textureId;
	}

	public void setTextureId(final int i) {
		textureId = i;
	}

	public boolean isNeighbor(final HE_Face f) {
		if (getHalfedge() == null) {
			return false;
		}
		HE_Halfedge he = getHalfedge();
		do {
			if (he.getPair() != null && he.getPair().getFace() != null && he.getPair().getFace() == f) {
				return true;
			}
			he = he.getNextInFace();
		} while (he != getHalfedge());
		return false;
	}

	@Override
	public WB_CoordCollection getPoints() {
		return WB_CoordCollection.getCollection(getFaceVertices());
	}

	@Override
	public void clearPrecomputed() {
		triangles = null;
	}

	public WB_CoordinateSystem getFaceCS() {
		return HE_MeshOp.getFaceCS(this);
	}

	public WB_Vector getFaceNormalNotNormalized() {
		return HE_MeshOp.getFaceNormalNotNormalized(this);
	}

	public WB_Classification getFaceType() {
		return HE_MeshOp.getFaceType(this);
	}

	public WB_Vector getNormalOffsetFaceCenter(final double d) {
		return HE_MeshOp.getNormalOffsetFaceCenter(this, d);
	}

	public WB_Plane getNormalOffsetPlane(final double d) {
		return HE_MeshOp.getNormalOffsetPlane(this, d);
	}

	public WB_Polygon getOrthoPolygon() {
		return HE_MeshOp.getOrthoPolygon(this);
	}

	public WB_Polygon getPlanarPolygon() {
		return HE_MeshOp.getPlanarPolygon(this);
	}

	public WB_Plane getPlane() {
		return HE_MeshOp.getPlane(this);
	}

	public WB_Polygon getPolygon() {
		return HE_MeshOp.getPolygon(this);
	}

	public WB_Triangle getTriangle() {
		return HE_MeshOp.getTriangle(this);
	}
}
