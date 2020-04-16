package wblut.hemesh;

import java.util.Iterator;
import java.util.List;

import wblut.core.WB_ProgressReporter.WB_ProgressCounter;
import wblut.geom.WB_AABB;
import wblut.geom.WB_Coord;
import wblut.geom.WB_Point;

public class HEM_TaubinSmooth extends HEM_Modifier {
	private boolean autoRescale;
	private boolean keepBoundary;
	private double lambda;
	private double mu;
	private int iter;

	public HEM_TaubinSmooth() {
		lambda = 0.5;
		mu = -0.52;
		iter = 1;
		keepBoundary = false;
	}

	public HEM_TaubinSmooth setAutoRescale(final boolean b) {
		autoRescale = b;
		return this;
	}

	public HEM_TaubinSmooth setIterations(final int r) {
		iter = r;
		return this;
	}

	public HEM_TaubinSmooth setKeepBoundary(final boolean b) {
		keepBoundary = b;
		return this;
	}

	public HEM_TaubinSmooth setLambda(final double lambda) {
		this.lambda = lambda;
		return this;
	}

	public HEM_TaubinSmooth setMu(final double mu) {
		this.mu = mu;
		return this;
	}

	@Override
	protected HE_Mesh applySelf(final HE_Mesh mesh) {
		tracker.setStartStatus(this, "Starting HEM_TaubinSmooth.");
		WB_AABB box = new WB_AABB();
		if (autoRescale) {
			box = HE_MeshOp.getAABB(mesh);
		}
		final WB_Coord[] newPositions = new WB_Coord[mesh.getNumberOfVertices()];
		if (iter < 1) {
			iter = 1;
		}
		final WB_ProgressCounter counter = new WB_ProgressCounter(iter * mesh.getNumberOfVertices(), 10);
		tracker.setCounterStatus(this, "Smoothing vertices.", counter);
		for (int r = 0; r < iter; r++) {
			final double f = r % 2 == 0 ? lambda : mu;
			Iterator<HE_Vertex> vItr = mesh.vItr();
			HE_Vertex v;
			List<HE_Vertex> neighbors;
			int id = 0;
			WB_Point p;
			while (vItr.hasNext()) {
				v = vItr.next();
				if (v.isBoundary() && keepBoundary) {
					newPositions[id] = v;
				} else {
					p = new WB_Point(v).mulSelf(1.0 - f);
					neighbors = v.getNeighborVertices();
					for (int i = 0; i < neighbors.size(); i++) {
						p.addMulSelf(f / neighbors.size(), neighbors.get(i));
					}
					newPositions[id] = p;
				}
				id++;
			}
			vItr = mesh.vItr();
			id = 0;
			while (vItr.hasNext()) {
				vItr.next().set(newPositions[id]);
				id++;
				counter.increment();
			}
		}
		if (autoRescale) {
			mesh.fitInAABB(box);
		}
		tracker.setStopStatus(this, "Exiting HEM_TaubinSmooth.");
		return mesh;
	}

	@Override
	protected HE_Mesh applySelf(final HE_Selection selection) {
		tracker.setStartStatus(this, "Starting HEM_TaubinSmooth.");
		selection.collectVertices();
		WB_AABB box = new WB_AABB();
		if (autoRescale) {
			box = HE_MeshOp.getAABB(selection.getParent());
		}
		final WB_Coord[] newPositions = new WB_Coord[selection.getNumberOfVertices()];
		if (iter < 1) {
			iter = 1;
		}
		final WB_ProgressCounter counter = new WB_ProgressCounter(iter * selection.getNumberOfVertices(), 10);
		tracker.setCounterStatus(this, "Smoothing vertices.", counter);
		for (int r = 0; r < iter; r++) {
			final double f = r % 2 == 0 ? lambda : mu;
			Iterator<HE_Vertex> vItr = selection.vItr();
			HE_Vertex v;
			HE_Vertex n;
			List<HE_Vertex> neighbors;
			int id = 0;
			while (vItr.hasNext()) {
				v = vItr.next();
				final WB_Point p = new WB_Point();
				if (v.isBoundary() && keepBoundary) {
					newPositions[id] = v;
				} else {
					neighbors = v.getNeighborVertices();
					final Iterator<HE_Vertex> nItr = neighbors.iterator();
					boolean selectionBoundary = false;
					while (nItr.hasNext()) {
						n = nItr.next();
						if (!selection.contains(n)) {
							if (keepBoundary) {
								selectionBoundary = true;
								break;
							} else {
								nItr.remove();
							}
						}
					}
					if (selectionBoundary) {
						newPositions[id] = v;
					} else {
						for (int i = 0; i < neighbors.size(); i++) {
							p.addMulSelf(f / neighbors.size(), neighbors.get(i));
						}
						newPositions[id] = p.addMulSelf(1.0 - f, v);
					}
				}
				id++;
			}
			vItr = selection.vItr();
			id = 0;
			while (vItr.hasNext()) {
				vItr.next().set(newPositions[id]);
				id++;
				counter.increment();
			}
		}
		if (autoRescale) {
			selection.getParent().fitInAABB(box);
		}
		tracker.setStopStatus(this, "Exiting HEM_TaubinSmooth.");
		return selection.getParent();
	}
}
