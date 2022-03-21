// 
// Decompiled by Procyon v0.5.36
// 

package wblut.external.poly2Tri.splayTree;

public class BTreeNode {
	protected SplayTreeItem _data;
	protected BTreeNode _left;
	protected BTreeNode _right;
	protected boolean _visited;

	public BTreeNode() {
		this._data = null;
		this._left = null;
		this._right = null;
		this._visited = false;
	}

	public BTreeNode(final SplayTreeItem data, final BTreeNode left, final BTreeNode right) {
		this._data = null;
		this._left = null;
		this._right = null;
		this._visited = false;
		this._data = data;
		this._left = left;
		this._right = right;
	}

	public SplayTreeItem data() {
		return this._data;
	}

	public BTreeNode left() {
		return this._left;
	}

	public BTreeNode right() {
		return this._right;
	}

	void setVisited(final boolean visited) {
		this._visited = visited;
	}

	boolean getVisited() {
		return this._visited;
	}

	Comparable keyValue() {
		return this._data.keyValue();
	}
}
