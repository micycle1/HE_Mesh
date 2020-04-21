package wblut.hemesh;

import java.util.Iterator;
import java.util.List;

import wblut.geom.WB_AABBTree3D;
import wblut.geom.WB_GeometryOp3D;
import wblut.geom.WB_List;
import wblut.geom.WB_Point;
import wblut.geom.WB_PolyLine;
import wblut.geom.WB_Segment;
import wblut.geom.WB_Vector;
import wblut.hemesh.HE_MeshOp.HE_FaceLineIntersection;

public class HEM_PolyLineInversion extends HEM_Modifier {
	private WB_PolyLine polyLine;
	private double r, r2;
	private double icutoff;
	private boolean linear;

	public HEM_PolyLineInversion() {
		super();
		icutoff = 0.0001;
		linear = false;
	}

	public HEM_PolyLineInversion(final WB_PolyLine poly, final double r) {
		super();
		polyLine = poly;
		this.r = r;
		r2 = r * r;
		icutoff = 0.0001;
		linear = false;
	}

	public HEM_PolyLineInversion setPolyLine(final WB_PolyLine poly) {
		polyLine = poly;
		return this;
	}

	public HEM_PolyLineInversion setRadius(final double r) {
		this.r = r;
		r2 = r * r;
		return this;
	}

	public HEM_PolyLineInversion setCutoff(final double cutoff) {
		icutoff = 1.0 / cutoff;
		return this;
	}

	public HEM_PolyLineInversion setLinear(final boolean b) {
		linear = b;
		return this;
	}

	@Override
	protected HE_Mesh applySelf(final HE_Mesh mesh) {
		if (polyLine == null) {
			return mesh;
		}
		if (r == 0) {
			return mesh;
		}
		final WB_AABBTree3D tree = new WB_AABBTree3D(mesh, 10);
		final List<HE_FaceLineIntersection> intersections = new WB_List<>();
		for (int i = 0; i < polyLine.getNumberSegments(); i++) {
			final WB_Segment S = polyLine.getSegment(i);
			intersections.addAll(HE_MeshOp.getIntersection(tree, S));
		}
		for (final HE_FaceLineIntersection fi : intersections) {
			if (mesh.contains(fi.face)) {
				mesh.deleteFace(fi.face);
			}
		}
		mesh.removeUnconnectedElements();
		HE_MeshOp.capHalfedges(mesh);
		final Iterator<HE_Vertex> vItr = mesh.vItr();
		HE_Vertex v;
		WB_Vector d;
		WB_Point surf;
		WB_Point q;
		double ri, rf;
		while (vItr.hasNext()) {
			v = vItr.next();
			q = WB_GeometryOp3D.getClosestPoint3D(v, polyLine);
			if (linear) {
				d = WB_Vector.subToVector3D(v, q);
				d.normalizeSelf();
				surf = q.addMulSelf(r, d);
				d = surf.subToVector3D(v).mulSelf(2);
				v.getPosition().addSelf(d);
			} else {
				d = WB_Vector.subToVector3D(v, q);
				ri = d.normalizeSelf();
				rf = r2 * Math.max(icutoff, 1.0 / ri);
				v.set(q);
				v.getPosition().addMulSelf(rf, d);
			}
		}
		return mesh;
	}

	@Override
	protected HE_Mesh applySelf(final HE_Selection selection) {
		if (polyLine == null) {
			return selection.getParent();
		}
		if (r == 0) {
			return selection.getParent();
		}
		final WB_AABBTree3D tree = new WB_AABBTree3D(selection.getParent(), 10);
		final List<HE_FaceLineIntersection> intersections = new WB_List<>();
		for (int i = 0; i < polyLine.getNumberSegments(); i++) {
			final WB_Segment S = polyLine.getSegment(i);
			intersections.addAll(HE_MeshOp.getIntersection(tree, S));
		}
		for (final HE_FaceLineIntersection fi : intersections) {
			if (selection.contains(fi.face)) {
				if (selection.getParent().contains(fi.face)) {
					selection.getParent().deleteFace(fi.face);
				}
			}
		}
		selection.getParent().removeUnconnectedElements();
		HE_MeshOp.capHalfedges(selection.getParent());
		selection.cleanSelection();
		final Iterator<HE_Vertex> vItr = selection.vItr();
		HE_Vertex v;
		WB_Vector d;
		WB_Point surf;
		WB_Point q;
		double ri, rf;
		while (vItr.hasNext()) {
			v = vItr.next();
			q = WB_GeometryOp3D.getClosestPoint3D(v, polyLine);
			if (linear) {
				d = WB_Vector.subToVector3D(v, q);
				d.normalizeSelf();
				surf = q.addMulSelf(r, d);
				d = surf.subToVector3D(v).mulSelf(2);
				v.getPosition().addSelf(d);
			} else {
				d = WB_Vector.subToVector3D(v, q);
				ri = d.normalizeSelf();
				rf = r2 * Math.max(icutoff, 1.0 / ri);
				v.set(q);
				v.getPosition().addMulSelf(rf, d);
			}
		}
		return selection.getParent();
	}
}
