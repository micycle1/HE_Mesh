package wblut.geom;

import lombok.AccessLevel;
import lombok.Data;
import lombok.Getter;
import lombok.Setter;
import lombok.ToString;

@Data
@ToString(includeFieldNames = true)
public class WB_Pentagon implements WB_Transformable3D {
	@Getter(AccessLevel.NONE)
	@Setter(AccessLevel.NONE)
	private final WB_GeometryFactory3D geometryfactory = new WB_GeometryFactory3D();
	@Setter(AccessLevel.NONE)
	private WB_Point p1;
	@Setter(AccessLevel.NONE)
	private WB_Point p2;
	@Setter(AccessLevel.NONE)
	private WB_Point p3;
	@Setter(AccessLevel.NONE)
	private WB_Point p4;
	@Setter(AccessLevel.NONE)
	private WB_Point p5;

	public WB_Pentagon(final WB_Coord p1, final WB_Coord p2, final WB_Coord p3, final WB_Coord p4, final WB_Coord p5) {
		this.p1 = geometryfactory.createPoint(p1);
		this.p2 = geometryfactory.createPoint(p2);
		this.p3 = geometryfactory.createPoint(p3);
		this.p4 = geometryfactory.createPoint(p4);
		this.p5 = geometryfactory.createPoint(p5);
	}

	public void cycle() {
		final WB_Point tmp = p1;
		p1 = p2;
		p2 = p3;
		p3 = p4;
		p4 = p5;
		p5 = tmp;
	}

	public void cycle(int n) {
		while (n >= 5) {
			n -= 5;
		}
		while (n < 0) {
			n += 5;
		}
		for (int i = 0; i < n; i++) {
			cycle();
		}
	}

	@Override
	public WB_Pentagon apply2D(final WB_Transform2D T) {
		return new WB_Pentagon(p1.apply2D(T), p2.apply2D(T), p3.apply2D(T), p4.apply2D(T), p5.apply2D(T));
	}

	@Override
	public WB_Pentagon apply2DSelf(final WB_Transform2D T) {
		p1.apply2DSelf(T);
		p2.apply2DSelf(T);
		p3.apply2DSelf(T);
		p4.apply2DSelf(T);
		p5.apply2DSelf(T);
		return this;
	}

	@Override
	public WB_Pentagon apply(final WB_Transform3D T) {
		return new WB_Pentagon(p1.apply(T), p2.apply(T), p3.apply(T), p4.apply(T), p5.apply(T));
	}

	@Override
	public WB_Pentagon applySelf(final WB_Transform3D T) {
		p1.applySelf(T);
		p2.applySelf(T);
		p3.applySelf(T);
		p4.applySelf(T);
		p5.applySelf(T);
		return this;
	}
}