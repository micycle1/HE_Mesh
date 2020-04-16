package wblut.hemesh;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.eclipse.collections.impl.list.mutable.FastList;

import wblut.core.WB_ProgressReporter.WB_ProgressCounter;
import wblut.geom.WB_AABBTree3D;
import wblut.geom.WB_Classification;
import wblut.geom.WB_GeometryOp3D;
import wblut.geom.WB_Plane;
import wblut.math.WB_Epsilon;

public class HEM_SliceSurface extends HEM_Modifier {
	static final int ON = 1, BACK = 2, FRONT = 3;
	private WB_Plane P;
	private HE_Selection cutFaces;
	private HE_Selection newEdges;
	private List<HE_Path> paths;

	public HEM_SliceSurface() {
		super();
	}

	public HEM_SliceSurface setPlane(final WB_Plane P) {
		this.P = P;
		return this;
	}

	public HEM_SliceSurface setPlane(final double ox, final double oy, final double oz, final double nx,
			final double ny, final double nz) {
		P = new WB_Plane(ox, oy, oz, nx, ny, nz);
		return this;
	}

	private double offset;

	public HEM_SliceSurface setOffset(final double d) {
		offset = d;
		return this;
	}

	@Override
	protected HE_Mesh applySelf(final HE_Mesh mesh) {
		tracker.setStartStatus(this, "Starting HEM_SliceSurface.");
		cutFaces = HE_Selection.getSelection(mesh);
		newEdges = HE_Selection.getSelection(mesh);
		mesh.resetEdgeInternalLabels();
		mesh.resetVertexInternalLabels();
		paths = new FastList<>();
		// no plane defined
		if (P == null) {
			tracker.setStopStatus(this, "No cutplane defined. Exiting HEM_SliceSurface.");
			return mesh;
		}
		// empty mesh
		if (mesh.getNumberOfVertices() == 0) {
			tracker.setStopStatus(this, "Empty mesh. Exiting HEM_SliceSurface.");
			return mesh;
		}
		// check if plane intersects mesh
		final WB_Plane lP = new WB_Plane(P.getNormal(), -P.d() + offset);
		if (!WB_GeometryOp3D.checkIntersection3D(HE_MeshOp.getAABB(mesh), lP)) {
			tracker.setStopStatus(this, "Plane doesn't intersect bounding box. Exiting HEM_SliceSurface.");
			return mesh;
		}
		tracker.setDuringStatus(this, "Creating bounding box tree.");
		final WB_AABBTree3D tree = new WB_AABBTree3D(mesh, Math.max(64, (int) Math.sqrt(mesh.getNumberOfFaces())));
		final HE_Selection faces = HE_Selection.getSelection(mesh);
		tracker.setDuringStatus(this, "Retrieving intersection candidates.");
		faces.addFaces(HE_MeshOp.getPotentialIntersectedFaces(tree, lP));
		faces.collectVertices();
		faces.collectEdgesByFace();
		WB_Classification tmp;
		final HashMap<Long, WB_Classification> vertexClass = new HashMap<>();
		WB_ProgressCounter counter = new WB_ProgressCounter(faces.getNumberOfVertices(), 10);
		tracker.setCounterStatus(this, "Classifying vertices.", counter);
		HE_Vertex v;
		final Iterator<HE_Vertex> vItr = faces.vItr();
		while (vItr.hasNext()) {
			v = vItr.next();
			tmp = WB_GeometryOp3D.classifyPointToPlane3D(v, lP);
			if (tmp == WB_Classification.ON) {
				v.setInternalLabel(ON);
			} else if (tmp == WB_Classification.BACK) {
				v.setInternalLabel(BACK);
			} else if (tmp == WB_Classification.FRONT) {
				v.setInternalLabel(FRONT);
			}
			vertexClass.put(v.getKey(), tmp);
			counter.increment();
		}
		counter = new WB_ProgressCounter(faces.getNumberOfEdges(), 10);
		tracker.setCounterStatus(this, "Classifying edges.", counter);
		final HE_Selection split = HE_Selection.getSelection(mesh);
		final Map<Long, Double> edgeInt = new HashMap<>();
		final Iterator<HE_Halfedge> eItr = faces.eItr();
		HE_Halfedge e;
		while (eItr.hasNext()) {
			e = eItr.next();
			if (vertexClass.get(e.getStartVertex().getKey()) == WB_Classification.ON) {
				if (vertexClass.get(e.getEndVertex().getKey()) == WB_Classification.ON) {
					newEdges.addEdge(e);
					e.setInternalLabel(1);
					e.getPair().setInternalLabel(1);
				} else {
					edgeInt.put(e.getKey(), 0.0);
				}
			} else if (vertexClass.get(e.getStartVertex().getKey()) == WB_Classification.BACK) {
				if (vertexClass.get(e.getEndVertex().getKey()) == WB_Classification.ON) {
					edgeInt.put(e.getKey(), 1.0);
				} else if (vertexClass.get(e.getEndVertex().getKey()) == WB_Classification.FRONT) {
					edgeInt.put(e.getKey(), HE_MeshOp.getIntersection(e, lP));
				}
			} else {
				if (vertexClass.get(e.getEndVertex().getKey()) == WB_Classification.ON) {
					edgeInt.put(e.getKey(), 1.0);
				} else if (vertexClass.get(e.getEndVertex().getKey()) == WB_Classification.BACK) {
					edgeInt.put(e.getKey(), HE_MeshOp.getIntersection(e, lP));
				}
			}
			counter.increment();
		}
		counter = new WB_ProgressCounter(edgeInt.size(), 10);
		tracker.setCounterStatus(this, "Indexing edge intersection.", counter);
		for (final Map.Entry<Long, Double> en : edgeInt.entrySet()) {
			final HE_Halfedge ce = mesh.getHalfedgeWithKey(en.getKey());
			final double u = en.getValue();
			if (ce.getFace() != null) {
				split.add(ce.getFace());
			}
			if (ce.getPair().getFace() != null) {
				split.add(ce.getPair().getFace());
			}
			if (u < WB_Epsilon.EPSILON) {
				split.add(ce.getStartVertex());
			} else if (u > 1.0 - WB_Epsilon.EPSILON) {
				split.add(ce.getEndVertex());
			} else {
				final HE_Vertex vi = HE_MeshOp.splitEdge(mesh, ce, u).vItr().next();
				vi.setInternalLabel(ON);
				split.add(vi);
			}
			counter.increment();
		}
		counter = new WB_ProgressCounter(split.getNumberOfFaces(), 10);
		tracker.setCounterStatus(this, "Splitting faces.", counter);
		HE_Face f;
		HE_FaceIterator fItr = split.fItr();
		while (fItr.hasNext()) {
			f = fItr.next();
			checkConsistency(f, lP);
			splitFace(f, mesh, lP);
			counter.increment();
		}
		tracker.setStopStatus(this, "Exiting HEM_SliceSurface.");
		// mesh.fixDegenerateTriangles();
		newEdges.cleanSelection();
		buildPaths(newEdges);
		fItr = mesh.fItr();
		while (fItr.hasNext()) {
			f = fItr.next();
			if (f.isDegenerate()) {
				// System.out.println(f.getFaceOrder());
			}
		}
		mesh.addSelection("cuts", this, cutFaces);
		mesh.addSelection("edges", this, newEdges);
		return mesh;
	}

	@Override
	protected HE_Mesh applySelf(final HE_Selection selection) {
		tracker.setStartStatus(this, "Starting HEM_SliceSurface.");
		selection.getParent().resetEdgeInternalLabels();
		selection.getParent().resetVertexInternalLabels();
		cutFaces = HE_Selection.getSelection(selection.getParent());
		newEdges = HE_Selection.getSelection(selection.getParent());
		paths = new FastList<>();
		// no plane defined
		if (P == null) {
			tracker.setStopStatus(this, "No cutplane defined. Exiting HEM_SliceSurface.");
			return selection.getParent();
		}
		// empty mesh
		if (selection.getParent().getNumberOfVertices() == 0) {
			tracker.setStopStatus(this, "Empty vertex selection. Exiting HEM_SliceSurface.");
			return selection.getParent();
		}
		final WB_Plane lP = new WB_Plane(P.getNormal(), -P.d() + offset);
		tracker.setDuringStatus(this, "Creating bounding box tree.");
		final WB_AABBTree3D tree = new WB_AABBTree3D(selection.getParent(), 64);
		final HE_Selection faces = HE_Selection.getSelection(selection.getParent());
		tracker.setDuringStatus(this, "Retrieving intersection candidates.");
		faces.addFaces(HE_MeshOp.getPotentialIntersectedFaces(tree, lP));
		final HE_Selection lsel = selection.get();
		lsel.intersect(faces);
		lsel.collectEdgesByFace();
		lsel.collectVertices();
		// empty mesh
		if (lsel.getNumberOfVertices() == 0) {
			tracker.setStopStatus(this, "Plane doesn't intersect bounding box tree. Exiting HEM_SliceSurface.");
			return lsel.getParent();
		}
		// check if plane intersects mesh
		boolean positiveVertexExists = false;
		boolean negativeVertexExists = false;
		WB_Classification tmp;
		final HE_ObjectMap<WB_Classification> vertexClass = new HE_ObjectMap<>();
		HE_Vertex v;
		WB_ProgressCounter counter = new WB_ProgressCounter(lsel.getNumberOfVertices(), 10);
		tracker.setCounterStatus(this, "Classifying vertices.", counter);
		final Iterator<HE_Vertex> vItr = lsel.vItr();
		while (vItr.hasNext()) {
			v = vItr.next();
			tmp = WB_GeometryOp3D.classifyPointToPlane3D(v, lP);
			vertexClass.put(v.getKey(), tmp);
			if (tmp == WB_Classification.FRONT) {
				positiveVertexExists = true;
			}
			if (tmp == WB_Classification.BACK) {
				negativeVertexExists = true;
			}
			counter.increment();
		}
		if (positiveVertexExists && negativeVertexExists) {
			new ArrayList<HE_Vertex>();
			final HE_Selection split = HE_Selection.getSelection(lsel.getParent());
			final HashMap<Long, Double> edgeInt = new HashMap<>();
			final Iterator<HE_Halfedge> eItr = lsel.eItr();
			HE_Halfedge e;
			counter = new WB_ProgressCounter(lsel.getNumberOfEdges(), 10);
			tracker.setCounterStatus(this, "Classifying edges.", counter);
			while (eItr.hasNext()) {
				e = eItr.next();
				if (vertexClass.get(e.getStartVertex().getKey()) == WB_Classification.ON) {
					if (vertexClass.get(e.getEndVertex().getKey()) == WB_Classification.ON) {
						newEdges.addEdge(e);
						e.setInternalLabel(1);
					} else {
						edgeInt.put(e.getKey(), 0.0);
					}
				} else if (vertexClass.get(e.getStartVertex().getKey()) == WB_Classification.BACK) {
					if (vertexClass.get(e.getEndVertex().getKey()) == WB_Classification.ON) {
						edgeInt.put(e.getKey(), 1.0);
					} else if (vertexClass.get(e.getEndVertex().getKey()) == WB_Classification.FRONT) {
						edgeInt.put(e.getKey(), HE_MeshOp.getIntersection(e, lP));
					}
				} else {
					if (vertexClass.get(e.getEndVertex().getKey()) == WB_Classification.ON) {
						edgeInt.put(e.getKey(), 1.0);
					} else if (vertexClass.get(e.getEndVertex().getKey()) == WB_Classification.BACK) {
						edgeInt.put(e.getKey(), HE_MeshOp.getIntersection(e, lP));
					}
				}
				counter.increment();
			}
			counter = new WB_ProgressCounter(edgeInt.size(), 10);
			tracker.setCounterStatus(this, "Indexing edge intersection.", counter);
			for (final Map.Entry<Long, Double> en : edgeInt.entrySet()) {
				final HE_Halfedge ce = lsel.getParent().getHalfedgeWithKey(en.getKey());
				final double u = en.getValue();
				if (lsel.contains(ce.getFace())) {
					split.add(ce.getFace());
				}
				if (lsel.contains(ce.getPair().getFace())) {
					split.add(ce.getPair().getFace());
				}
				if (u < WB_Epsilon.EPSILON) {
					split.add(ce.getStartVertex());
				} else if (u > 1.0 - WB_Epsilon.EPSILON) {
					split.add(ce.getEndVertex());
				} else {
					final HE_Vertex vi = HE_MeshOp.splitEdge(lsel.getParent(), ce, u).vItr().next();
					vi.setInternalLabel(ON);
					split.add(vi);
				}
				counter.increment();
			}
			HE_Face f;
			counter = new WB_ProgressCounter(split.getNumberOfFaces(), 10);
			tracker.setCounterStatus(this, "Splitting faces.", counter);
			final HE_FaceIterator fItr = split.fItr();
			while (fItr.hasNext()) {
				f = fItr.next();
				checkConsistency(f, lP);
				splitFace(f, lsel.getParent(), lP);
				counter.increment();
			}
			paths = new FastList<>();
		}
		if (newEdges.getNumberOfEdges() > 1) {
			buildPaths(newEdges);
		}
		lsel.getParent().addSelection("cuts", this, cutFaces);
		selection.addFaces(cutFaces);
		lsel.getParent().addSelection("edges", this, newEdges);
		tracker.setStopStatus(this, "Exiting HEM_SliceSurface.");
		return lsel.getParent();
	}

	private void buildPaths(final HE_Selection cutEdges) {
		tracker.setDuringStatus(this, "Building slice paths.");
		if (cutEdges.getNumberOfEdges() == 0) {
			return;
		}
		final List<HE_Halfedge> edges = new HE_HalfedgeList();
		for (final HE_Halfedge he : cutEdges.getEdges()) {
			final HE_Face f = he.getFace();
			if (f != null) {
				if (WB_GeometryOp3D.classifyPointToPlane3D(HE_MeshOp.getFaceCenter(f), P) == WB_Classification.FRONT) {
					edges.add(he.getPair());
				} else {
					edges.add(he);
				}
			}
		}
		final WB_ProgressCounter counter = new WB_ProgressCounter(edges.size(), 10);
		tracker.setCounterStatus(this, "Processing slice edges.", counter);
		while (edges.size() > 0) {
			final List<HE_Halfedge> pathedges = new HE_HalfedgeList();
			HE_Halfedge current = edges.get(0);
			pathedges.add(current);
			boolean loop = false;
			for (int i = 0; i < edges.size(); i++) {
				if (edges.get(i).getVertex() == current.getEndVertex()) {
					if (i > 0) {
						current = edges.get(i);
						pathedges.add(current);
						i = -1;
					} else {
						loop = true;
						break;
					}
				}
			}
			if (!loop) {
				final List<HE_Halfedge> reversepathedges = new HE_HalfedgeList();
				current = edges.get(0);
				for (int i = 0; i < edges.size(); i++) {
					if (edges.get(i).getEndVertex() == current.getVertex()) {
						if (i > 0) {
							current = edges.get(i);
							reversepathedges.add(current);
							i = 0;
						}
					}
				}
				final List<HE_Halfedge> finalpathedges = new HE_HalfedgeList();
				for (int i = reversepathedges.size() - 1; i > -1; i--) {
					finalpathedges.add(reversepathedges.get(i));
				}
				finalpathedges.addAll(pathedges);
				paths.add(new HE_Path(finalpathedges, loop));
				edges.removeAll(finalpathedges);
			} else {
				paths.add(new HE_Path(pathedges, loop));
				edges.removeAll(pathedges);
			}
			counter.increment(pathedges.size());
		}
	}

	public List<HE_Path> getPaths() {
		return paths;
	}

	void checkConsistency(final HE_Face f, final WB_Plane P) {
		final HE_FaceEdgeCirculator feCrc = f.feCrc();
		HE_Halfedge e;
		while (feCrc.hasNext()) {
			e = feCrc.next();
			if (e.getVertex().getInternalLabel() == FRONT && e.getEndVertex().getInternalLabel() == BACK
					|| e.getVertex().getInternalLabel() == BACK && e.getEndVertex().getInternalLabel() == FRONT) {
				final double e0 = WB_GeometryOp3D.getDistanceToPlane3D(e.getVertex(), P);
				final double e1 = WB_GeometryOp3D.getDistanceToPlane3D(e.getEndVertex(), P);
				if (e0 < e1) {
					e.getVertex().setInternalLabel(ON);
					e.getVertex().set(WB_GeometryOp3D.projectOnPlane(e.getVertex(), P));
				} else {
					e.getEndVertex().setInternalLabel(ON);
					e.getEndVertex().set(WB_GeometryOp3D.projectOnPlane(e.getEndVertex(), P));
				}
			}
		}
	}

	void splitFace(final HE_Face f, final HE_Mesh mesh, final WB_Plane P) {
		int intersectionCount = 0;
		final HE_FaceVertexCirculator fvCrc = f.fvCrc();
		while (fvCrc.hasNext()) {
			if (fvCrc.next().getInternalLabel() == ON) {
				intersectionCount++;
			}
		}
		if (intersectionCount < 2) {
			return;
		} else {
			final List<HE_Vertex[]> subPolygons = new HET_FaceSplitter().splitFace(f, P);
			if (subPolygons.size() > 1) {
				final HE_HalfedgeList allhalfedges = new HE_HalfedgeList();
				final HE_ObjectMap<HE_TextureCoordinate> UVWs = new HE_ObjectMap<>();
				final List<HE_Halfedge> hes = f.getFaceHalfedges();
				for (final HE_Halfedge he : hes) {
					UVWs.put(he.getVertex().getKey(), he.getUVW());
				}
				for (final HE_Vertex[] subPoly : subPolygons) {
					final HE_HalfedgeList halfedges = new HE_HalfedgeList();
					HE_Halfedge he;
					final HE_Face subFace = new HE_Face();
					subFace.copyProperties(f);
					for (int j = 0; j < subPoly.length; j++) {
						he = new HE_Halfedge();
						if (subPoly[j] != subPoly[(j + 1) % subPoly.length]) {
							he.setUVW(UVWs.get(subPoly[j].getKey()));
							mesh.setVertex(he, subPoly[j]);
							mesh.setHalfedge(subPoly[j], he);
							mesh.setFace(he, subFace);
							halfedges.add(he);
							allhalfedges.add(he);
						}
					}
					if (halfedges.size() > 2) {
						mesh.setHalfedge(subFace, halfedges.get(0));
						for (int j = 0, k = halfedges.size() - 1; j < halfedges.size(); k = j, j++) {
							mesh.setNext(halfedges.get(k), halfedges.get(j));
						}
						mesh.addDerivedElement(subFace, f);
						cutFaces.add(subFace);
						mesh.addHalfedges(halfedges);
					}
				}
				mesh.cutFace(f);
				HE_MeshOp.pairHalfedges(mesh);
				for (final HE_Halfedge he : allhalfedges) {
					if (he.isEdge()) {
						if (he.getVertex().getInternalLabel() == ON && he.getEndVertex().getInternalLabel() == ON) {
							newEdges.addEdge(he);
							he.setInternalLabel(1);
						}
					}
				}
			}
		}
	}

	public static void main(final String[] args) {
		final HEC_Torus creator = new HEC_Torus(80, 200, 6, 16);
		final HE_Mesh mesh = new HE_Mesh(creator);
		final HEM_SliceSurface modifier = new HEM_SliceSurface();
		final WB_Plane P = new WB_Plane(0, 0, 0, 0, 0, 1);
		modifier.setPlane(P);
		modifier.setOffset(0);
		mesh.modify(modifier);
	}
}
