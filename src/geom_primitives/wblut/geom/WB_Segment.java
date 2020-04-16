package wblut.geom;

import java.util.List;

import org.eclipse.collections.impl.list.mutable.FastList;

import wblut.math.WB_Math;

public class WB_Segment extends WB_Line implements Comparable<WB_Segment> {
	protected double length;
	protected WB_Point endpoint;
	private final WB_GeometryFactory3D geometryfactory = new WB_GeometryFactory3D();

	public WB_Segment() {
		super();
		endpoint = new WB_Point();
		length = 0;
	}

	public static WB_Segment X() {
		return new WB_Segment(0, 0, 0, 1, 0, 0);
	}

	public static WB_Segment Y() {
		return new WB_Segment(0, 0, 0, 0, 1, 0);
	}

	public static WB_Segment Z() {
		return new WB_Segment(0, 0, 0, 0, 0, 1);
	}

	public WB_Segment(final WB_Coord o, final WB_Coord d, final double l) {
		super(o, d);
		length = l;
		endpoint = new WB_Point(direction).mulSelf(l).addSelf(origin);
	}

	public WB_Segment(final WB_Coord p1, final WB_Coord p2) {
		super(p1, new WB_Vector(p1, p2));
		endpoint = new WB_Point(p2);
		length = Math.sqrt(WB_GeometryOp3D.getSqDistance3D(p1, p2));
	}

	public WB_Segment(final double p1x, final double p1y, final double p1z, final double p2x, final double p2y,
			final double p2z) {
		super(new WB_Point(p1x, p1y, p1z), new WB_Vector(p2x - p1x, p2y - p1y, p2z - p1z));
		endpoint = new WB_Point(p2x, p2y, p2z);
		length = Math.sqrt(WB_GeometryOp3D.getSqDistance3D(origin, endpoint));
	}

	@Override
	public WB_Point getParametricPoint(final double t) {
		final WB_Point result = new WB_Point(direction);
		result.scaleSelf(WB_Math.clamp(t, 0, 1) * length);
		result.addSelf(origin);
		return result;
	}

	@Override
	public void getParametricPointInto(final double t, final WB_MutableCoord result) {
		result.set(new WB_Vector(direction).mulSelf(WB_Math.clamp(t, 0, 1) * length).addSelf(origin));
	}

	public WB_Coord getEndpoint() {
		return endpoint;
	}

	public WB_Point getCenter() {
		return new WB_Point(endpoint).addSelf(origin).mulSelf(0.5);
	}

	public double getLength() {
		return length;
	}

	public static List<WB_Segment> negate(final List<WB_Segment> segs) {
		final List<WB_Segment> neg = new FastList<>();
		for (int i = 0; i < segs.size(); i++) {
			neg.add(segs.get(i).negate());
		}
		return neg;
	}

	public WB_Segment negate() {
		return new WB_Segment(endpoint, origin);
	}

	public void reverse() {
		final WB_Coord tmp = origin;
		origin = new WB_Point(endpoint);
		endpoint = new WB_Point(tmp);
		direction.mulSelf(-1);
	}

	@Override
	public WB_Point getPointOnCurve(final double u) {
		if (u < 0 || u > 1) {
			return null;
		}
		return this.getParametricPoint(u);
	}

	@Override
	public double getLowerU() {
		return 0;
	}

	@Override
	public double getUpperU() {
		return 1;
	}

	@Override
	public WB_Vector getDirectionOnCurve(final double u) {
		return new WB_Vector(direction);
	}

	@Override
	public WB_Vector getDerivative(final double u) {
		return new WB_Vector(direction);
	}

	@Override
	public boolean equals(final Object o) {
		if (o == this) {
			return true;
		}
		if (!(o instanceof WB_Segment)) {
			return false;
		}
		return origin.equals(((WB_Segment) o).getOrigin()) && endpoint.equals(((WB_Segment) o).getEndpoint());
	}

	@Override
	public int hashCode() {
		return 31 * origin.hashCode() + endpoint.hashCode();
	}

	@Override
	public WB_Segment apply2D(final WB_Transform2D T) {
		return geometryfactory.createSegment(new WB_Point(origin).applyAsPoint2D(T),
				new WB_Point(endpoint).applyAsPoint2D(T));
	}

	@Override
	public WB_Segment apply2DSelf(final WB_Transform2D T) {
		origin.apply2DSelf(T);
		endpoint.apply2DSelf(T);
		direction.apply2DSelf(T);
		length = Math.sqrt(WB_GeometryOp3D.getSqDistance3D(origin, endpoint));
		return this;
	}

	@Override
	public WB_Segment apply(final WB_Transform3D T) {
		return geometryfactory.createSegment(new WB_Point(origin).applyAsPoint(T),
				new WB_Point(endpoint).applyAsPoint(T));
	}

	@Override
	public WB_Segment applySelf(final WB_Transform3D T) {
		origin.applySelf(T);
		endpoint.applySelf(T);
		direction.applySelf(T);
		length = Math.sqrt(WB_GeometryOp3D.getSqDistance3D(origin, endpoint));
		return this;
	}

	@Override
	public int compareTo(final WB_Segment seg) {
		final WB_Coord a = getCenter();
		final WB_Coord b = seg.getCenter();
		return a.compareTo(b);
	}
}