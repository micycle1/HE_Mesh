package wblut.geom;

public class WB_Tetrahedron implements WB_Transformable3D {
	WB_Point p1;
	WB_Point p2;
	WB_Point p3;
	WB_Point p4;

	protected WB_Tetrahedron() {
	}

	private final WB_GeometryFactory3D geometryfactory = new WB_GeometryFactory3D();

	public WB_Tetrahedron(final WB_Coord p1, final WB_Coord p2, final WB_Coord p3, final WB_Coord p4) {
		this.p1 = geometryfactory.createPoint(p1);
		this.p2 = geometryfactory.createPoint(p2);
		this.p3 = geometryfactory.createPoint(p3);
		this.p4 = geometryfactory.createPoint(p4);
	}

	public WB_Coord p1() {
		return p1;
	}

	public WB_Coord p2() {
		return p2;
	}

	public WB_Coord p3() {
		return p3;
	}

	public WB_Coord p4() {
		return p4;
	}

	public WB_Coord getPoint(final int i) {
		if (i == 0) {
			return p1;
		} else if (i == 1) {
			return p2;
		} else if (i == 2) {
			return p3;
		} else if (i == 3) {
			return p4;
		}
		return null;
	}

	public WB_Point getCenter() {
		return geometryfactory.createMidpoint(p1, p2, p3, p3);
	}

	public double getVolume() {
		final WB_Vector a = geometryfactory.createVectorFromTo(p1, p4);
		final WB_Vector b = geometryfactory.createVectorFromTo(p2, p4);
		final WB_Vector c = geometryfactory.createVectorFromTo(p3, p4);
		return Math.abs(a.dot(b.crossSelf(c))) / 6.0;
	}

	public double getCircumradius() {
		final WB_Vector a = geometryfactory.createVectorFromTo(p1, p4);
		final WB_Vector b = geometryfactory.createVectorFromTo(p2, p4);
		final WB_Vector c = geometryfactory.createVectorFromTo(p3, p4);
		final WB_Vector O = b.cross(c).mulSelf(a.dot(a));
		O.addSelf(c.cross(a).mulSelf(b.dot(b)));
		O.addSelf(a.cross(b).mulSelf(c.dot(c)));
		O.mulSelf(1.0 / (2 * a.dot(b.crossSelf(c))));
		return O.getLength();
	}

	public WB_Point getCircumcenter() {
		final WB_Vector a = geometryfactory.createVectorFromTo(p1, p4);
		final WB_Vector b = geometryfactory.createVectorFromTo(p2, p4);
		final WB_Vector c = geometryfactory.createVectorFromTo(p3, p4);
		final WB_Vector O = b.cross(c).mulSelf(a.dot(a));
		O.addSelf(c.cross(a).mulSelf(b.dot(b)));
		O.addSelf(a.cross(b).mulSelf(c.dot(c)));
		O.mulSelf(1.0 / (2 * a.dot(b.crossSelf(c))));
		return p4.add(O);
	}

	public WB_Sphere getCircumsphere() {
		final WB_Vector a = geometryfactory.createVectorFromTo(p1, p4);
		final WB_Vector b = geometryfactory.createVectorFromTo(p2, p4);
		final WB_Vector c = geometryfactory.createVectorFromTo(p3, p4);
		final WB_Vector O = b.cross(c).mulSelf(a.dot(a));
		O.addSelf(c.cross(a).mulSelf(b.dot(b)));
		O.addSelf(a.cross(b).mulSelf(c.dot(c)));
		O.mulSelf(1.0 / (2 * a.dot(b.crossSelf(c))));
		return geometryfactory.createSphereWithRadius(p4.add(O), O.getLength());
	}

	public double getInradius() {
		final WB_Vector a = geometryfactory.createVectorFromTo(p1, p4);
		final WB_Vector b = geometryfactory.createVectorFromTo(p2, p4);
		final WB_Vector c = geometryfactory.createVectorFromTo(p3, p4);
		final WB_Vector bXc = b.cross(c);
		final double sixV = Math.abs(a.dot(bXc));
		c.crossSelf(a);
		a.crossSelf(b);
		final double denom = bXc.getLength() + c.getLength() + a.getLength() + bXc.addMulSelf(2, a).getLength();
		return sixV / denom;
	}

	public WB_Point getIncenter() {
		final WB_Vector a = geometryfactory.createVectorFromTo(p1, p4);
		final WB_Vector b = geometryfactory.createVectorFromTo(p2, p4);
		final WB_Vector c = geometryfactory.createVectorFromTo(p3, p4);
		final WB_Vector bXc = b.cross(c);
		final WB_Vector cXa = c.cross(a);
		final WB_Vector aXb = a.cross(b);
		final double bXcLength = bXc.getLength();
		final double cXaLength = cXa.getLength();
		final double aXbLength = aXb.getLength();
		final double dLength = bXc.addSelf(cXa).addSelf(aXb).getLength();
		final WB_Vector O = a.mulSelf(bXcLength);
		O.addSelf(b.mulSelf(cXaLength));
		O.addSelf(c.mulSelf(aXbLength));
		O.divSelf(bXcLength + cXaLength + aXbLength + dLength);
		return p4.add(O);
	}

	public WB_Sphere getInsphere() {
		final WB_Vector a = geometryfactory.createVectorFromTo(p1, p4);
		final WB_Vector b = geometryfactory.createVectorFromTo(p2, p4);
		final WB_Vector c = geometryfactory.createVectorFromTo(p3, p4);
		final WB_Vector bXc = b.cross(c);
		final WB_Vector cXa = c.cross(a);
		final WB_Vector aXb = a.cross(b);
		final double bXcLength = bXc.getLength();
		final double cXaLength = cXa.getLength();
		final double aXbLength = aXb.getLength();
		final double dLength = bXc.addSelf(cXa).addSelf(aXb).getLength();
		final WB_Vector O = a.mulSelf(bXcLength);
		O.addSelf(b.mulSelf(cXaLength));
		O.addSelf(c.mulSelf(aXbLength));
		O.divSelf(bXcLength + cXaLength + aXbLength + dLength);
		return geometryfactory.createSphereWithRadius(p4.add(O), O.getLength());
	}

	public boolean isAcute() {
		return WB_GeometryOp3D.getCosDihedralAngle(p1, p2, p3, p4) > 0.0
				&& WB_GeometryOp3D.getCosDihedralAngle(p1, p2, p4, p3) > 0.0
				&& WB_GeometryOp3D.getCosDihedralAngle(p1, p3, p4, p2) > 0.0
				&& WB_GeometryOp3D.getCosDihedralAngle(p3, p1, p2, p4) > 0.0
				&& WB_GeometryOp3D.getCosDihedralAngle(p2, p1, p3, p4) > 0.0
				&& WB_GeometryOp3D.getCosDihedralAngle(p2, p1, p4, p3) > 0.0;
	}

	@Override
	public WB_Tetrahedron apply2D(final WB_Transform2D T) {
		return new WB_Tetrahedron(p1.apply2D(T), p2.apply2D(T), p3.apply2D(T), p4.apply2D(T));
	}

	@Override
	public WB_Tetrahedron apply2DSelf(final WB_Transform2D T) {
		p1.apply2DSelf(T);
		p2.apply2DSelf(T);
		p3.apply2DSelf(T);
		p4.apply2DSelf(T);
		return this;
	}

	@Override
	public WB_Tetrahedron apply(final WB_Transform3D T) {
		return new WB_Tetrahedron(p1.apply(T), p2.apply(T), p3.apply(T), p4.apply(T));
	}

	@Override
	public WB_Tetrahedron applySelf(final WB_Transform3D T) {
		p1.applySelf(T);
		p2.applySelf(T);
		p3.applySelf(T);
		p4.applySelf(T);
		return this;
	}

	public WB_AABB getAABB() {
		final WB_AABB aabb = new WB_AABB();
		aabb.expandToInclude(p1());
		aabb.expandToInclude(p2());
		aabb.expandToInclude(p3());
		aabb.expandToInclude(p4());
		return aabb;
	}
}