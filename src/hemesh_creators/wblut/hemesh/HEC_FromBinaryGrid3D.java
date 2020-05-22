package wblut.hemesh;

import java.util.List;

import wblut.geom.WB_BinaryGrid3D;
import wblut.geom.WB_List;
import wblut.geom.WB_Point;
import wblut.geom.WB_Quad;

/**
 *
 */
public class HEC_FromBinaryGrid3D extends HEC_Creator {
	/**  */
	private WB_BinaryGrid3D grid;

	/**
	 *
	 *
	 * @param grid
	 */
	public HEC_FromBinaryGrid3D(final WB_BinaryGrid3D grid) {
		this.grid = grid;
	}

	/**
	 *
	 *
	 * @param grid
	 * @return
	 */
	public HEC_FromBinaryGrid3D setGrid(final WB_BinaryGrid3D grid) {
		this.grid = grid;
		return this;
	}

	/**
	 *
	 *
	 * @return
	 */
	@Override
	public HE_Mesh createBase() {
		final List<WB_Quad> quads = new WB_List<>();
		this.setCenter(grid.getCenter());
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
						quads.add(new WB_Quad(new WB_Point(x, y, z), new WB_Point(x, y + grid.getDY(), z),
								new WB_Point(x, y + grid.getDY(), z + grid.getDZ()),
								new WB_Point(x, y, z + grid.getDZ())));
					}
				}
			}
		}
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
						quads.add(new WB_Quad(new WB_Point(x, y, z), new WB_Point(x + grid.getDX(), y, z),
								new WB_Point(x + grid.getDX(), y, z + grid.getDZ()),
								new WB_Point(x, y, z + grid.getDZ())));
					}
				}
			}
		}
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
						quads.add(new WB_Quad(new WB_Point(x, y, z), new WB_Point(x + grid.getDX(), y, z),
								new WB_Point(x + grid.getDX(), y + grid.getDY(), z),
								new WB_Point(x, y + grid.getDY(), z)));
					}
				}
			}
		}
		return new HEC_FromQuads().setQuads(quads).createBase();
	}
}
