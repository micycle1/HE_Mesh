// 
// Decompiled by Procyon v0.5.36
// 

package wblut.external.poly2Tri;

public class Poly2TriUtils {
	public static final double PI = 3.141592653589793;
	public static final int UNKNOWN = 1;
	public static final int INPUT = 2;
	public static final int INSERT = 3;
	public static final int START = 4;
	public static final int END = 5;
	public static final int MERGE = 6;
	public static final int SPLIT = 7;
	public static final int REGULAR_UP = 8;
	public static final int REGULAR_DOWN = 9;
	public static int l_id;
	public static int p_id;

	static {
		Poly2TriUtils.l_id = 0;
		Poly2TriUtils.p_id = 0;
	}

	public static String typeToString(final int type) {
		switch (type) {
			case 1 : {
				return "UNKNOWN";
			}
			case 2 : {
				return "INPUT";
			}
			case 3 : {
				return "INERT";
			}
			case 4 : {
				return "START";
			}
			case 5 : {
				return "END";
			}
			case 6 : {
				return "MERGE";
			}
			case 7 : {
				return "SPLIT";
			}
			case 8 : {
				return "REGULAR_UP";
			}
			case 9 : {
				return "REGULAR_DOWN";
			}
			default : {
				return "??? (" + type + ")";
			}
		}
	}

	public static double orient2d(final double[] pa, final double[] pb, final double[] pc) {
		final double detleft = (pa[0] - pc[0]) * (pb[1] - pc[1]);
		final double detright = (pa[1] - pc[1]) * (pb[0] - pc[0]);
		return detleft - detright;
	}

	public static void initPoly2TriUtils() {
		Poly2TriUtils.l_id = 0;
		Poly2TriUtils.p_id = 0;
	}
}
