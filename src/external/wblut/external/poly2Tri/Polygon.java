// 
// Decompiled by Procyon v0.5.36
// 

package wblut.external.poly2Tri;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.Stack;

import wblut.external.poly2Tri.splayTree.BTreeNode;
import wblut.external.poly2Tri.splayTree.SplayTree;

public class Polygon {
	protected int _ncontours;
	protected int[] _nVertices;
	protected HashMap<Integer, Pointbase> _points;
	protected int[] _pointsKeys;
	protected HashMap<Integer, Linebase> _edges;
	protected int[] _edgesKeys;
	private PriorityQueue<Pointbase> _qpoints;
	private SplayTree _edgebst;
	private ArrayList<ArrayList<Integer>> _mpolys;
	private ArrayList<ArrayList<Integer>> _triangles;
	private HashMap<Integer, Set<Integer>> _startAdjEdgeMap;
	private HashMap _diagonals;
	private boolean _debug;
	private FileWriter _logfile;
	private UpdateKey updateKey;
	private String _debugFileName;

	public Map<Integer, Pointbase> points() {
		return this._points;
	}

	public Map<Integer, Linebase> edges() {
		return this._edges;
	}

	private void initPolygon(final int numContours, final int[] numVerticesInContours, final double[][] vertices) {
		int nextNumber = 1;
		this._ncontours = numContours;
		this._nVertices = new int[this._ncontours];
		for (int i = 0; i < numContours; ++i) {
			for (int j = 0; j < numVerticesInContours[i]; ++j) {
				this._points.put(nextNumber, new Pointbase(nextNumber, vertices[nextNumber - 1][0], vertices[nextNumber - 1][1], 2));
				++nextNumber;
			}
		}
		this._nVertices[0] = numVerticesInContours[0];
		for (int i = 1; i < this._ncontours; ++i) {
			this._nVertices[i] = this._nVertices[i - 1] + numVerticesInContours[i];
		}
		int i = 0;
		int j = 1;
		int first = 1;
		while (i < this._ncontours) {
			while (j + 1 <= this._nVertices[i]) {
				final Linebase edge = new Linebase(this._points.get(j), this._points.get(j + 1), 2);
				this._edges.put(Poly2TriUtils.l_id, edge);
				++j;
			}
			final Linebase edge = new Linebase(this._points.get(j), this._points.get(first), 2);
			this._edges.put(Poly2TriUtils.l_id, edge);
			j = this._nVertices[i] + 1;
			first = this._nVertices[i] + 1;
			++i;
		}
		Poly2TriUtils.p_id = this._nVertices[this._ncontours - 1];
	}

	Polygon(final int numContures, final int[] numVerticesInContures, final double[][] vertices) {
		this._ncontours = 0;
		this._nVertices = null;
		this._points = new HashMap<>();
		this._pointsKeys = null;
		this._edges = new HashMap<>();
		this._edgesKeys = null;
		this._qpoints = new PriorityQueue<>(30, new PointbaseComparatorCoordinatesReverse());
		this._edgebst = new SplayTree();
		this._mpolys = new ArrayList<>();
		this._triangles = new ArrayList<>();
		this._startAdjEdgeMap = new HashMap<>();
		this._diagonals = new HashMap<>();
		this._debug = false;
		this._logfile = null;
		this.updateKey = new UpdateKey();
		this._debugFileName = "polygon_triangulation_log.txt";
		Poly2TriUtils.initPoly2TriUtils();
		this.initPolygon(numContures, numVerticesInContures, vertices);
		this.initializate();
		this._debug = false;
	}

	public void writeToLog(final String s) {
		if (!this._debug) {
			return;
		}
		try {
			this._logfile.write(s);
		} catch (IOException e) {
			this._debug = false;
			System.out.println("Writing to LogFile (debugging) failed.");
			e.printStackTrace();
			System.out.println("Setting _debug = false, continuing the work.");
		}
	}

	public Pointbase getPoint(final int index) {
		return this._points.get(index);
	}

	public Linebase getEdge(final int index) {
		return this._edges.get(index);
	}

	public Pointbase qpointsTop() {
		return this._qpoints.peek();
	}

	public Pointbase qpointsPop() {
		return this._qpoints.poll();
	}

	public void destroy() {
	}

	public boolean is_exist(final double x, final double y) {
		final Iterator<Integer> iter = this._points.keySet().iterator();
		while (iter.hasNext()) {
			final Pointbase pb = this.getPoint(iter.next());
			if (pb.x == x && pb.y == y) {
				return true;
			}
		}
		return false;
	}

	private int prev(final int i) {
		int j = 0;
		int prevLoop = 0;
		int currentLoop;
		for (currentLoop = 0; i > this._nVertices[currentLoop]; ++currentLoop) {
			prevLoop = currentLoop;
		}
		if (i == 1 || i == this._nVertices[prevLoop] + 1) {
			j = this._nVertices[currentLoop];
		} else if (i <= this._nVertices[currentLoop]) {
			j = i - 1;
		}
		return j;
	}

	private int next(final int i) {
		int j = 0;
		int prevLoop = 0;
		int currentLoop;
		for (currentLoop = 0; i > this._nVertices[currentLoop]; ++currentLoop) {
			prevLoop = currentLoop;
		}
		if (i < this._nVertices[currentLoop]) {
			j = i + 1;
		} else if (i == this._nVertices[currentLoop]) {
			if (currentLoop == 0) {
				j = 1;
			} else {
				j = this._nVertices[prevLoop] + 1;
			}
		}
		return j;
	}

	private int[] getSorted(final Set s) {
		final Object[] temp = s.toArray();
		final int[] result = new int[temp.length];
		for (int i = 0; i < temp.length; ++i) {
			result[i] = (int) temp[i];
		}
		Arrays.sort(result);
		return result;
	}

	private void initializePointsKeys() {
		this._pointsKeys = this.getSorted(this._points.keySet());
	}

	private void initializeEdgesKeys() {
		this._edgesKeys = this.getSorted(this._edges.keySet());
	}

	private Set<Integer> getSetFromStartAdjEdgeMap(final int index) {
		Set<Integer> s = this._startAdjEdgeMap.get(index);
		if (s != null) {
			return s;
		}
		s = new HashSet<>();
		this._startAdjEdgeMap.put(index, s);
		return s;
	}

	private void initializate() {
		this.initializePointsKeys();
		for (final int id : this._pointsKeys) {
			final int idp = this.prev(id);
			final int idn = this.next(id);
			final Pointbase p = this.getPoint(id);
			final Pointbase pnext = this.getPoint(idn);
			final Pointbase pprev = this.getPoint(idp);
			if (p.compareTo(pnext) > 0 && pprev.compareTo(p) > 0) {
				p.type = 9;
			} else if (p.compareTo(pprev) > 0 && pnext.compareTo(p) > 0) {
				p.type = 8;
			} else {
				final double area = Poly2TriUtils.orient2d(new double[] { pprev.x, pprev.y }, new double[] { p.x, p.y },
						new double[] { pnext.x, pnext.y });
				if (pprev.compareTo(p) > 0 && pnext.compareTo(p) > 0) {
					p.type = ((area > 0.0) ? 5 : 6);
				}
				if (pprev.compareTo(p) < 0 && pnext.compareTo(p) < 0) {
					p.type = ((area > 0.0) ? 4 : 7);
				}
			}
			this._qpoints.add(new Pointbase(p));
			this.getSetFromStartAdjEdgeMap(id).add(id);
		}
	}

	private void addDiagonal(final int i, final int j) {
		final int type = 3;
		final Linebase diag = new Linebase(this.getPoint(i), this.getPoint(j), type);
		this._edges.put(diag.id(), diag);
		this.getSetFromStartAdjEdgeMap(i).add(diag.id());
		this.getSetFromStartAdjEdgeMap(j).add(diag.id());
		this._diagonals.put(diag.id(), diag);
		this.writeToLog("Add Diagonal from " + i + " to " + j + "\n");
	}

	private void handleStartVertex(final int i) {
		final double y = this._points.get(i).y;
		this._edgebst.inOrder(this.updateKey, y);
		final Linebase edge = this.getEdge(i);
		edge.setHelper(i);
		edge.setKeyValue(y);
		this._edgebst.insert(edge);
		if (this._debug) {
			this.writeToLog("set e" + i + " helper to " + i + "\n");
			this.writeToLog("Insert e" + i + " to splay tree\n");
			this.writeToLog("key:" + edge.keyValue() + "\n");
		}
	}

	private void handleEndVertex(final int i) {
		final double y = this.getPoint(i).y;
		this._edgebst.inOrder(this.updateKey, y);
		final int previ = this.prev(i);
		final Linebase edge = this.getEdge(previ);
		final int helper = edge.helper();
		if (this.getPoint(helper).type == 6) {
			this.addDiagonal(i, helper);
		}
		this._edgebst.delete(edge.keyValue());
		if (this._debug) {
			this.writeToLog("Remove e" + previ + " from splay tree\n");
			this.writeToLog("key:" + edge.keyValue() + "\n");
		}
	}

	private void handleSplitVertex(final int i) {
		final Pointbase point = this.getPoint(i);
		final double x = point.x;
		final double y = point.y;
		this._edgebst.inOrder(this.updateKey, y);
		final BTreeNode leftnode = this._edgebst.findMaxSmallerThan(x);
		final Linebase leftedge = (Linebase) leftnode.data();
		final int helper = leftedge.helper();
		this.addDiagonal(i, helper);
		if (this._debug) {
			this.writeToLog("Search key:" + x + " edge key:" + leftedge.keyValue() + "\n");
			this.writeToLog("e" + leftedge.id() + " is directly left to v" + i + "\n");
			this.writeToLog("Set e" + leftedge.id() + " helper to " + i + "\n");
			this.writeToLog("set e" + i + " helper to " + i + "\n");
			this.writeToLog("Insert e" + i + " to splay tree\n");
			this.writeToLog("Insert key:" + this.getEdge(i).keyValue() + "\n");
		}
		leftedge.setHelper(i);
		final Linebase edge = this.getEdge(i);
		edge.setHelper(i);
		edge.setKeyValue(y);
		this._edgebst.insert(edge);
	}

	private void handleMergeVertex(final int i) {
		final Pointbase point = this.getPoint(i);
		final double x = point.x;
		final double y = point.y;
		this._edgebst.inOrder(this.updateKey, y);
		final int previ = this.prev(i);
		final Linebase previEdge = this.getEdge(previ);
		int helper = previEdge.helper();
		Pointbase helperPoint = this.getPoint(helper);
		if (helperPoint.type == 6) {
			this.addDiagonal(i, helper);
		}
		this._edgebst.delete(previEdge.keyValue());
		if (this._debug) {
			this.writeToLog("e" + previ + " helper is " + helper + "\n");
			this.writeToLog("Remove e" + previ + " from splay tree.\n");
		}
		final BTreeNode leftnode = this._edgebst.findMaxSmallerThan(x);
		final Linebase leftedge = (Linebase) leftnode.data();
		helper = leftedge.helper();
		helperPoint = this.getPoint(helper);
		if (helperPoint.type == 6) {
			this.addDiagonal(i, helper);
		}
		leftedge.setHelper(i);
		if (this._debug) {
			this.writeToLog("Search key:" + x + " found:" + leftedge.keyValue() + "\n");
			this.writeToLog("e" + leftedge.id() + " is directly left to v" + i + "\n");
			this.writeToLog("Set e" + leftedge.id() + " helper to " + i + "\n");
		}
	}

	private void handleRegularVertexDown(final int i) {
		final Pointbase point = this.getPoint(i);
		final double y = point.y;
		this._edgebst.inOrder(this.updateKey, y);
		final int previ = this.prev(i);
		final Linebase previEdge = this.getEdge(previ);
		final int helper = previEdge.helper();
		final Pointbase helperPoint = this.getPoint(helper);
		if (helperPoint.type == 6) {
			this.addDiagonal(i, helper);
		}
		this._edgebst.delete(previEdge.keyValue());
		final Linebase edge = this.getEdge(i);
		edge.setHelper(i);
		edge.setKeyValue(y);
		this._edgebst.insert(edge);
		if (this._debug) {
			this.writeToLog("e" + previ + " helper is " + helper + "\n");
			this.writeToLog("Remove e" + previ + " from splay tree.\n");
			this.writeToLog("Set e" + i + " helper to " + i + "\n");
			this.writeToLog("Insert e" + i + " to splay tree\n");
			this.writeToLog("Insert key:" + edge.keyValue() + "\n");
		}
	}

	private void handleRegularVertexUp(final int i) {
		final Pointbase point = this.getPoint(i);
		final double x = point.x;
		final double y = point.y;
		this._edgebst.inOrder(this.updateKey, y);
		final BTreeNode leftnode = this._edgebst.findMaxSmallerThan(x);
		final Linebase leftedge = (Linebase) leftnode.data();
		final int helper = leftedge.helper();
		final Pointbase helperPoint = this.getPoint(helper);
		if (helperPoint.type == 6) {
			this.addDiagonal(i, helper);
		}
		leftedge.setHelper(i);
		if (this._debug) {
			this.writeToLog("Search key:" + x + " found:" + leftedge.keyValue() + "\n");
			this.writeToLog("e" + leftedge.id() + " is directly left to v" + i + " and its helper is:" + helper + "\n");
			this.writeToLog("Set e" + leftedge.id() + " helper to " + i + "\n");
		}
	}

	public boolean partition2Monotone() {
		if (this.qpointsTop().type != 4) {
			System.out.println("Please check your input polygon:\n1)orientations?\n2)duplicated points?\n");
			System.out.println("poly2tri stopped.\n");
			return false;
		}
		while (this._qpoints.size() > 0) {
			final Pointbase vertex = this.qpointsPop();
			final int id = vertex.id;
			if (this._debug) {
				String stype = null;
				switch (vertex.type) {
					case 4 : {
						stype = "START";
						break;
					}
					case 5 : {
						stype = "END";
						break;
					}
					case 6 : {
						stype = "MERGE";
						break;
					}
					case 7 : {
						stype = "SPLIT";
						break;
					}
					case 8 : {
						stype = "REGULAR_UP";
						break;
					}
					case 9 : {
						stype = "REGULAR_DOWN";
						break;
					}
					default : {
						System.out.println("No duplicated points please! poly2tri stopped\n");
						return false;
					}
				}
				this.writeToLog("\n\nHandle vertex:" + vertex.id + " type:" + stype + "\n");
			}
			switch (vertex.type) {
				case 4 : {
					this.handleStartVertex(id);
					continue;
				}
				case 5 : {
					this.handleEndVertex(id);
					continue;
				}
				case 6 : {
					this.handleMergeVertex(id);
					continue;
				}
				case 7 : {
					this.handleSplitVertex(id);
					continue;
				}
				case 8 : {
					this.handleRegularVertexUp(id);
					continue;
				}
				case 9 : {
					this.handleRegularVertexDown(id);
					continue;
				}
				default : {
					System.out.println("No duplicated points please! poly2tri stopped\n");
					return false;
				}
			}
		}
		return true;
	}

	private double angleCosb(final double[] pa, final double[] pb, final double[] pc) {
		final double dxab = pa[0] - pb[0];
		final double dyab = pa[1] - pb[1];
		final double dxcb = pc[0] - pb[0];
		final double dycb = pc[1] - pb[1];
		final double dxab2 = dxab * dxab;
		final double dyab2 = dyab * dyab;
		final double dxcb2 = dxcb * dxcb;
		final double dycb2 = dycb * dycb;
		final double ab = dxab2 + dyab2;
		final double cb = dxcb2 + dycb2;
		double cosb = dxab * dxcb + dyab * dycb;
		final double denom = Math.sqrt(ab * cb);
		cosb /= denom;
		return cosb;
	}

	private int selectNextEdge(final Linebase edge) {
		final int eid = edge.endPoint(1).id;
		final Set<Integer> edges = this.getSetFromStartAdjEdgeMap(eid);
		assert edges.size() != 0;
		int nexte = 0;
		if (edges.size() == 1) {
			nexte = edges.iterator().next();
		} else {
			final int[] edgesKeys = this.getSorted(edges);
			int nexte_ccw = 0;
			int nexte_cw = 0;
			double max = -2.0;
			double min = 2.0;
			for (final int it : edges) {
				if (it == edge.id()) {
					continue;
				}
				final Linebase iEdge = this.getEdge(it);
				final double[] A = { 0.0, 0.0 };
				final double[] B = { 0.0, 0.0 };
				final double[] C = { 0.0, 0.0 };
				A[0] = edge.endPoint(0).x;
				A[1] = edge.endPoint(0).y;
				B[0] = edge.endPoint(1).x;
				B[1] = edge.endPoint(1).y;
				if (!edge.endPoint(1).equals(iEdge.endPoint(0))) {
					iEdge.reverse();
				}
				C[0] = iEdge.endPoint(1).x;
				C[1] = iEdge.endPoint(1).y;
				final double area = Poly2TriUtils.orient2d(A, B, C);
				final double cosb = this.angleCosb(A, B, C);
				if (area > 0.0 && max < cosb) {
					nexte_ccw = it;
					max = cosb;
				} else {
					if (min <= cosb) {
						continue;
					}
					nexte_cw = it;
					min = cosb;
				}
			}
			nexte = ((nexte_ccw != 0) ? nexte_ccw : nexte_cw);
		}
		return nexte;
	}

	public boolean searchMonotones() {
		int loop = 0;
		final HashMap<Integer, Linebase> edges = (HashMap<Integer, Linebase>) this._edges.clone();
		while (edges.size() > this._diagonals.size()) {
			++loop;
			final ArrayList<Integer> poly = new ArrayList<>();
			final int[] edgesKeys = this.getSorted(edges.keySet());
			final int it = edgesKeys[0];
			final Linebase itEdge = edges.get(it);
			final Pointbase startp = itEdge.endPoint(0);
			Pointbase endp = null;
			Linebase next = itEdge;
			poly.add(startp.id);
			if (this._debug) {
				this.writeToLog("Searching for loops:" + loop + "\n");
				this.writeToLog("vertex index:" + startp.id + " ");
			}
			while (true) {
				endp = next.endPoint(1);
				if (next.type() != 3) {
					edges.remove(next.id());
					this.getSetFromStartAdjEdgeMap(next.endPoint(0).id).remove(next.id());
				}
				if (endp == startp) {
					this.writeToLog("\nloop closed!\n\n");
					this._mpolys.add(poly);
					break;
				}
				poly.add(endp.id);
				this.writeToLog(String.valueOf(endp.id) + " ");
				final int nexte = this.selectNextEdge(next);
				if (nexte == 0) {
					System.out.println("Please check your input polygon:\n");
					System.out.println("1)orientations?\n2)with duplicated points?\n3)is a simple one?\n");
					System.out.println("poly2tri stopped.\n");
					return false;
				}
				next = edges.get(nexte);
				if (next.endPoint(0).equals(endp)) {
					continue;
				}
				next.reverse();
			}
		}
		return true;
	}

	private void triangulateMonotone(final ArrayList<Integer> mpoly) {
		final PriorityQueue<Pointbase> qvertex = new PriorityQueue<>(30, new PointbaseComparatorCoordinatesReverse());
		for (int it = 0; it < mpoly.size(); ++it) {
			int itnext = it + 1;
			if (itnext == mpoly.size()) {
				itnext = 0;
			}
			final Pointbase point = new Pointbase(this.getPoint(mpoly.get(it)));
			final Pointbase pointnext = new Pointbase(this.getPoint(mpoly.get(itnext)));
			point.left = (point.compareTo(pointnext) > 0);
			qvertex.add(point);
		}
		final Stack<Pointbase> spoint = new Stack<>();
		for (int i = 0; i < 2; ++i) {
			spoint.push(qvertex.poll());
		}
		final double[] pa = { 0.0, 0.0 };
		final double[] pb = { 0.0, 0.0 };
		final double[] pc = { 0.0, 0.0 };
		while (qvertex.size() > 1) {
			final Pointbase topQueuePoint = qvertex.peek();
			final Pointbase topStackPoint = spoint.peek();
			if (topQueuePoint.left != topStackPoint.left) {
				while (spoint.size() > 1) {
					final Pointbase p1 = spoint.peek();
					spoint.pop();
					final Pointbase p2 = spoint.peek();
					final ArrayList<Integer> v = new ArrayList<>(3);
					v.add(topQueuePoint.id - 1);
					v.add(p1.id - 1);
					v.add(p2.id - 1);
					this._triangles.add(v);
					this.writeToLog("Add triangle:" + (v.get(0) + 1) + " " + (v.get(1) + 1) + " " + (v.get(2) + 1) + "\n");
				}
				spoint.pop();
				spoint.push(topStackPoint);
				spoint.push(topQueuePoint);
			} else {
				while (spoint.size() > 1) {
					final Pointbase stack1Point = spoint.peek();
					spoint.pop();
					final Pointbase stack2Point = spoint.peek();
					spoint.push(stack1Point);
					pa[0] = topQueuePoint.x;
					pa[1] = topQueuePoint.y;
					pb[0] = stack2Point.x;
					pb[1] = stack2Point.y;
					pc[0] = stack1Point.x;
					pc[1] = stack1Point.y;
					if (this._debug) {
						this.writeToLog("current top queue vertex index=" + topQueuePoint.id + "\n");
						this.writeToLog("Current top stack vertex index=" + stack1Point.id + "\n");
						this.writeToLog("Second stack vertex index=" + stack2Point.id + "\n");
					}
					final double area = Poly2TriUtils.orient2d(pa, pb, pc);
					final boolean left = stack1Point.left;
					if ((area <= 0.0 || !left) && (area >= 0.0 || left)) {
						break;
					}
					final ArrayList<Integer> v = new ArrayList<>(3);
					v.add(topQueuePoint.id - 1);
					v.add(stack2Point.id - 1);
					v.add(stack1Point.id - 1);
					this._triangles.add(v);
					this.writeToLog("Add triangle:" + (v.get(0) + 1) + " " + (v.get(1) + 1) + " " + (v.get(2) + 1) + "\n");
					spoint.pop();
				}
				spoint.push(topQueuePoint);
			}
			qvertex.poll();
		}
		final Pointbase lastQueuePoint = qvertex.peek();
		while (spoint.size() != 1) {
			final Pointbase topPoint = spoint.peek();
			spoint.pop();
			final Pointbase top2Point = spoint.peek();
			final ArrayList<Integer> v = new ArrayList<>(3);
			v.add(lastQueuePoint.id - 1);
			v.add(topPoint.id - 1);
			v.add(top2Point.id - 1);
			this._triangles.add(v);
			this.writeToLog("Add triangle:" + (v.get(0) + 1) + " " + (v.get(1) + 1) + " " + (v.get(2) + 1) + "\n");
		}
	}

	public boolean triangulation() {
		if (!this.partition2Monotone() || !this.searchMonotones()) {
			return false;
		}
		for (ArrayList<Integer> element : this._mpolys) {
			this.triangulateMonotone(element);
		}
		this.setDebugOption(false);
		return true;
	}

	public ArrayList<ArrayList<Integer>> triangles() {
		return this._triangles;
	}

	public void setDebugOption(final boolean debug) {
		if (debug == this._debug) {
			return;
		}
		if (this._debug) {
			try {
				this._logfile.close();
			} catch (IOException e) {
				System.out.println("Problem closing logfile.");
				e.printStackTrace();
				System.out.println("Continueing the work");
			}
		} else {
			try {
				this._logfile = new FileWriter(this._debugFileName);
			} catch (IOException e) {
				System.out.println("Error creating file polygon_triangulation_log.txt, switchin debug off, continuing.");
				e.printStackTrace();
				this._debug = false;
			}
		}
		this._debug = debug;
	}

	public void setDebugFile(final String debugFileName) {
		this._debugFileName = debugFileName;
	}
}
