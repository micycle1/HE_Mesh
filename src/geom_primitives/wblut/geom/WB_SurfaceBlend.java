package wblut.geom;

import wblut.math.WB_Math;

/**
 *
 */
public class WB_SurfaceBlend implements WB_Surface {
	/**  */
	private final WB_Surface surfA;
	/**  */
	private final WB_Surface surfB;

	/**
	 *
	 *
	 * @param surfA
	 * @param surfB
	 */
	public WB_SurfaceBlend(final WB_Surface surfA, final WB_Surface surfB) {
		this.surfA = surfA;
		this.surfB = surfB;
	}

	/**
	 *
	 *
	 * @return
	 */
	@Override
	public double getLowerU() {
		return WB_Math.max(surfA.getLowerU(), surfB.getLowerU());
	}

	/**
	 *
	 *
	 * @return
	 */
	@Override
	public double getLowerV() {
		return WB_Math.max(surfA.getLowerV(), surfB.getLowerV());
	}

	/**
	 *
	 *
	 * @param u
	 * @param v
	 * @return
	 */
	@Override
	public WB_Point surfacePoint(final double u, final double v) {
		return surfA.surfacePoint(u, v).addSelf(surfB.surfacePoint(u, v)).mulSelf(0.5);
	}

	/**
	 *
	 *
	 * @param u
	 * @param v
	 * @param t
	 * @return
	 */
	public WB_Point surfacePoint(final double u, final double v, final double t) {
		if (t == 0) {
			return surfA.surfacePoint(u, v);
		}
		if (t == 1) {
			return surfB.surfacePoint(u, v);
		}
		final WB_Point A = surfA.surfacePoint(u, v);
		return A.addMulSelf(t, surfB.surfacePoint(u, v).subSelf(A));
	}

	/**
	 *
	 *
	 * @return
	 */
	@Override
	public double getUpperU() {
		return WB_Math.min(surfA.getUpperU(), surfB.getUpperU());
	}

	/**
	 *
	 *
	 * @return
	 */
	@Override
	public double getUpperV() {
		return WB_Math.min(surfA.getUpperV(), surfB.getUpperV());
	}
}
