package wblut.hemesh;

import java.util.Iterator;
import java.util.List;

import wblut.core.WB_ProgressReporter.WB_ProgressCounter;
import wblut.geom.WB_AABB;
import wblut.geom.WB_Coord;
import wblut.geom.WB_Point;

public class HEM_AntiSmooth extends HEM_Modifier {
	public HEM_AntiSmooth() {
		parameters.set("autorescale", false);
		parameters.set("lambda", 0.5);
		parameters.set("preserveboundary", false);
		parameters.set("iter", 1);
	}

	public HEM_AntiSmooth setAutoRescale(final boolean b) {
		parameters.set("autorescale", b);
		return this;
	}

	public HEM_AntiSmooth setIterations(final int r) {
		parameters.set("iter", r);
		return this;
	}

	public HEM_AntiSmooth setPreserveBoundary(final boolean b) {
		parameters.set("preserveboundary", b);
		return this;
	}

	public HEM_AntiSmooth setLambda(final double lambda) {
		parameters.set("lambda", lambda);
		return this;
	}

	protected boolean getAutoRescale() {
		return parameters.get("autorescale", false);
	}

	protected int getIterations() {
		return parameters.get("iter", 1);
	}

	protected boolean getPreserveBoundary() {
		return parameters.get("preserveboundary", false);
	}

	protected double getLambda() {
		return parameters.get("lambda", 0.5);
	}

	@Override
	protected HE_Mesh applySelf(final HE_Mesh mesh) {
		tracker.setStartStatus(this, "Starting HEM_AntiSmooth.");
		WB_AABB box = new WB_AABB();
		int iter = getIterations();
		final boolean autoRescale = getAutoRescale();
		final boolean preserveBoundary = getPreserveBoundary();
		final double lambda = getLambda();
		if (autoRescale) {
			box = HE_MeshOp.getAABB(mesh);
		}
		final WB_Coord[] newPositions = new WB_Coord[mesh.getNumberOfVertices()];
		if (iter < 1) {
			iter = 1;
		}
		final WB_ProgressCounter counter = new WB_ProgressCounter(iter * mesh.getNumberOfVertices(), 10);
		tracker.setCounterStatus(this, "Anti smoothing vertices.", counter);
		for (int r = 0; r < iter; r++) {
			Iterator<HE_Vertex> vItr = mesh.vItr();
			HE_Vertex v;
			List<HE_Vertex> neighbors;
			int id = 0;
			WB_Point p;
			while (vItr.hasNext()) {
				v = vItr.next();
				if (v.isBoundary() && preserveBoundary) {
					newPositions[id] = v;
				} else {
					p = new WB_Point(v).mulSelf(1.0 - lambda);
					neighbors = v.getNeighborVertices();
					for (int i = 0; i < neighbors.size(); i++) {
						p.addMulSelf(lambda / neighbors.size(), neighbors.get(i));
					}
					newPositions[id] = new WB_Point(v).mulSelf(2.0).subSelf(p);
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
		tracker.setStopStatus(this, "Exiting HEM_AntiSmooth.");
		return mesh;
	}

	@Override
	protected HE_Mesh applySelf(final HE_Selection selection) {
		tracker.setStartStatus(this, "Starting HEM_AntiSmooth.");
		selection.collectVertices();
		WB_AABB box = new WB_AABB();
		int iter = getIterations();
		final boolean autoRescale = getAutoRescale();
		final boolean preserveBoundary = getPreserveBoundary();
		final double lambda = getLambda();
		if (autoRescale) {
			box = HE_MeshOp.getAABB(selection.getParent());
		}
		final WB_Coord[] newPositions = new WB_Coord[selection.getNumberOfVertices()];
		if (iter < 1) {
			iter = 1;
		}
		final WB_ProgressCounter counter = new WB_ProgressCounter(iter * selection.getNumberOfVertices(), 10);
		tracker.setCounterStatus(this, "Anti smoothing vertices.", counter);
		for (int r = 0; r < iter; r++) {
			Iterator<HE_Vertex> vItr = selection.vItr();
			HE_Vertex v;
			HE_Vertex n;
			List<HE_Vertex> neighbors;
			int id = 0;
			while (vItr.hasNext()) {
				v = vItr.next();
				final WB_Point p = new WB_Point();
				if (v.isBoundary() && preserveBoundary) {
					newPositions[id] = v;
				} else {
					neighbors = v.getNeighborVertices();
					final Iterator<HE_Vertex> nItr = neighbors.iterator();
					while (nItr.hasNext()) {
						n = nItr.next();
						if (!selection.contains(n)) {
							nItr.remove();
						}
					}
					for (int i = 0; i < neighbors.size(); i++) {
						p.addMulSelf(lambda / neighbors.size(), neighbors.get(i));
					}
					p.addMulSelf(1.0 - lambda, v);
					newPositions[id] = new WB_Point(v).mulSelf(2.0).subSelf(p);
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
		tracker.setStopStatus(this, "Exiting HEM_AntiSmooth.");
		return selection.getParent();
	}
}
