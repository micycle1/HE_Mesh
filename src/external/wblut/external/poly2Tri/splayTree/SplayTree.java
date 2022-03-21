// 
// Decompiled by Procyon v0.5.36
// 

package wblut.external.poly2Tri.splayTree;

public class SplayTree
{
    private BTreeNode root;
    private long size;
    private static BTreeNode header;
    
    static {
        SplayTree.header = new BTreeNode();
    }
    
    public SplayTree() {
        this.root = null;
        this.size = 0L;
    }
    
    public SplayTree(final SplayTree rhs) {
        final SplayTree st = clone(rhs);
        this.root = st.root;
        this.size = st.size;
    }
    
    void makeEmpty() {
        this.root = null;
        this.size = 0L;
    }
    
    boolean isEmpty() {
        return this.size == 0L;
    }
    
    long size() {
        return this.size;
    }
    
    public BTreeNode root() {
        return this.root;
    }
    
    public void insert(final SplayTreeItem x) {
        final BTreeNode newNode = new BTreeNode();
        newNode._data = x;
        if (this.root == null) {
            this.root = newNode;
            ++this.size;
            return;
        }
        Comparable keys = x.keyValue();
        while (true) {
            this.root = this.splay(keys, this.root);
            final Comparable rootk = this.root.keyValue();
            if (keys.compareTo(rootk) < 0) {
                newNode._left = this.root._left;
                newNode._right = this.root;
                this.root._left = null;
                this.root = newNode;
                ++this.size;
                return;
            }
            if (keys.compareTo(rootk) > 0) {
                newNode._right = this.root._right;
                newNode._left = this.root;
                this.root._right = null;
                this.root = newNode;
                ++this.size;
                return;
            }
            x.increaseKeyValue(1.0E-10);
            keys = x.keyValue();
        }
    }
    
    public BTreeNode delete(final Comparable keys) {
        this.root = this.splay(keys, this.root);
        if (!this.root.keyValue().equals(keys)) {
            return null;
        }
        final BTreeNode result = this.root;
        BTreeNode newTree;
        if (this.root._left == null) {
            newTree = this.root._right;
        }
        else {
            newTree = this.root._left;
            newTree = this.splay(keys, newTree);
            newTree._right = this.root._right;
        }
        this.root = newTree;
        --this.size;
        return result;
    }
    
    public BTreeNode deleteMax() {
        if (this.isEmpty()) {
            return null;
        }
        final double keys = Double.MAX_VALUE;
        this.root = this.splay(keys, this.root);
        final BTreeNode maxResult = this.root;
        BTreeNode newTree;
        if (this.root._left == null) {
            newTree = this.root._right;
        }
        else {
            newTree = this.root._left;
            newTree = this.splay(keys, newTree);
            newTree._right = this.root._right;
        }
        --this.size;
        this.root = newTree;
        return maxResult;
    }
    
    public static SplayTree clone(final SplayTree rhs) {
        final SplayTree st = new SplayTree();
        st.root = rhs.clone(rhs.root);
        st.size = rhs.size;
        return st;
    }
    
    public BTreeNode find(final Comparable keys) {
        if (this.isEmpty()) {
            return null;
        }
        this.root = this.splay(keys, this.root);
        if (!this.root.keyValue().equals(keys)) {
            return null;
        }
        return this.root;
    }
    
    public BTreeNode findMaxSmallerThan(final Comparable keys) {
        if (this.isEmpty()) {
            return null;
        }
        this.root = this.splay(keys, this.root);
        if (this.root.data().keyValue().compareTo(keys) < 0) {
            return this.root;
        }
        if (this.root._left != null) {
            BTreeNode result;
            for (result = this.root._left; result._right != null; result = result._right) {}
            return result;
        }
        assert false;
        return null;
    }
    
    public void inOrder(final SplayTreeAction action, final double y) {
        this.inOrder(action, this.root, y);
    }
    
    public int height() {
        return this.height(this.root);
    }
    
    public int height(final BTreeNode t) {
        if (t == null) {
            return 0;
        }
        int lh = this.height(t._left);
        int rh = this.height(t._right);
        return (lh > rh) ? (++lh) : (++rh);
    }
    
    public BTreeNode left(final BTreeNode node) {
        return node.left();
    }
    
    public BTreeNode right(final BTreeNode node) {
        return node.right();
    }
    
    private BTreeNode clone(final BTreeNode t) {
        if (t == null) {
            return null;
        }
        if (t == t._left) {
            return null;
        }
        return new BTreeNode(t._data, this.clone(t._left), this.clone(t._right));
    }
    
    private void inOrder(final SplayTreeAction action, final BTreeNode t, final double y) {
        if (t != null) {
            this.inOrder(action, t._left, y);
            action.action(t, y);
            this.inOrder(action, t._right, y);
        }
    }
    
    private BTreeNode rotateWithLeftChild(final BTreeNode k2) {
        final BTreeNode k3 = k2._left;
        k2._left = k3._right;
        k3._right = k2;
        return k3;
    }
    
    private BTreeNode rotateWithRightChild(final BTreeNode k1) {
        final BTreeNode k2 = k1._right;
        k1._right = k2._left;
        k2._left = k1;
        return k2;
    }
    
    private BTreeNode splay(final Comparable keys, BTreeNode t) {
        final BTreeNode header = SplayTree.header;
        final BTreeNode header2 = SplayTree.header;
        final BTreeNode bTreeNode = null;
        header2._right = bTreeNode;
        header._left = bTreeNode;
        BTreeNode _leftTreeMax;
        BTreeNode _rightTreeMin = _leftTreeMax = SplayTree.header;
        while (true) {
            final Comparable rKey = t.keyValue();
            if (keys.compareTo(rKey) < 0) {
                if (t._left == null) {
                    break;
                }
                if (keys.compareTo(t._left.keyValue()) < 0) {
                    t = this.rotateWithLeftChild(t);
                }
                if (t._left == null) {
                    break;
                }
                _rightTreeMin._left = t;
                _rightTreeMin = t;
                t = t._left;
            }
            else {
                if (keys.compareTo(rKey) <= 0) {
                    break;
                }
                if (t._right == null) {
                    break;
                }
                if (keys.compareTo(t._right.keyValue()) > 0) {
                    t = this.rotateWithRightChild(t);
                }
                if (t._right == null) {
                    break;
                }
                _leftTreeMax._right = t;
                _leftTreeMax = t;
                t = t._right;
            }
        }
        _leftTreeMax._right = t._left;
        _rightTreeMin._left = t._right;
        t._left = SplayTree.header._right;
        t._right = SplayTree.header._left;
        return t;
    }
}
