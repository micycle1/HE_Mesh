// 
// Decompiled by Procyon v0.5.36
// 

package wblut.external.poly2Tri.testPoly2Tri;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Polygon;
import java.util.ArrayList;

import javax.swing.JPanel;

public class Poly2TriPainting extends JPanel {
	public ArrayList<double[]> polygons;
	public double maxX;
	public double maxY;
	public double minX;
	public double minY;

	public Poly2TriPainting() {
		this.polygons = new ArrayList<>();
		this.maxX = 10.0;
		this.maxY = 10.0;
		this.minX = -1.0;
		this.minY = -1.0;
	}

	public void addPolygon(final double[] xy) {
		this.polygons.add(xy);
	}

	protected int[] getPoint(final double x, final double y) {
		final double xSize = this.getWidth() / (this.maxX - this.minX);
		final double ySize = this.getHeight() / (this.maxY - this.minY);
		final double[] center = { this.minX + (this.maxX - this.minX) / 2.0, this.minY + (this.maxY - this.minY) / 2.0 };
		final int[] realCenter = { this.getWidth() / 2, this.getHeight() / 2 };
		final int[] point = { (int) Math.round(realCenter[0] - center[0] * xSize + x * xSize),
				(int) Math.round(realCenter[1] + center[1] * ySize - y * ySize) };
		return point;
	}

	@Override
	public void paint(final Graphics g) {
		for (final double[] xy : this.polygons) {
			final Polygon p = new Polygon();
			final int[] point1 = this.getPoint(xy[0], xy[1]);
			final int[] point2 = this.getPoint(xy[2], xy[3]);
			final int[] point3 = this.getPoint(xy[4], xy[5]);
			p.addPoint(point1[0], point1[1]);
			p.addPoint(point2[0], point2[1]);
			p.addPoint(point3[0], point3[1]);
			g.setColor(Color.BLUE);
			g.fillPolygon(p);
			g.setColor(Color.BLACK);
			g.drawPolygon(p);
		}
	}
}
