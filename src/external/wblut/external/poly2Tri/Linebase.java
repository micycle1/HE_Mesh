// 
// Decompiled by Procyon v0.5.36
// 

package wblut.external.poly2Tri;

import wblut.external.poly2Tri.splayTree.SplayTreeItem;

public class Linebase implements SplayTreeItem
{
    protected int _id;
    protected Pointbase[] _endp;
    protected int _type;
    protected double _key;
    protected int _helper;
    
    public Linebase() {
        this._id = -1;
        this._endp = new Pointbase[2];
        this._type = 1;
        this._key = 0.0;
        this._helper = -1;
        for (int i = 0; i < 2; ++i) {
            this._endp[i] = null;
        }
        this._id = 0;
    }
    
    public Linebase(final Pointbase ep1, final Pointbase ep2, final int iType) {
        this._id = -1;
        this._endp = new Pointbase[2];
        this._type = 1;
        this._key = 0.0;
        this._helper = -1;
        this._endp[0] = ep1;
        this._endp[1] = ep2;
        this._id = ++Poly2TriUtils.l_id;
        this._type = iType;
    }
    
    public Linebase(final Linebase line) {
        this._id = -1;
        this._endp = new Pointbase[2];
        this._type = 1;
        this._key = 0.0;
        this._helper = -1;
        this._id = line._id;
        this._endp[0] = line._endp[0];
        this._endp[1] = line._endp[1];
        this._key = line._key;
        this._helper = line._helper;
    }
    
    public int id() {
        return this._id;
    }
    
    public Pointbase endPoint(final int i) {
        return this._endp[i];
    }
    
    public int type() {
        return this._type;
    }
    
    @Override
    public Comparable keyValue() {
        return this._key;
    }
    
    public void setKeyValue(final double y) {
        if (this._endp[1].y == this._endp[0].y) {
            this._key = ((this._endp[0].x < this._endp[1].x) ? this._endp[0].x : this._endp[1].x);
        }
        else {
            this._key = (y - this._endp[0].y) * (this._endp[1].x - this._endp[0].x) / (this._endp[1].y - this._endp[0].y) + this._endp[0].x;
        }
    }
    
    public void reverse() {
        assert this._type == 3;
        final Pointbase tmp = this._endp[0];
        this._endp[0] = this._endp[1];
        this._endp[1] = tmp;
    }
    
    public void setHelper(final int i) {
        this._helper = i;
    }
    
    public int helper() {
        return this._helper;
    }
    
    @Override
    public String toString() {
        final StringBuffer sb = new StringBuffer();
        sb.append("Linebase(");
        sb.append("ID = " + this._id);
        sb.append(", " + Poly2TriUtils.typeToString(this._type));
        sb.append(", [");
        sb.append(this._endp[0]);
        sb.append(", ");
        sb.append(this._endp[1]);
        sb.append("], type = " + this._type);
        sb.append(", keyValue =" + this.keyValue());
        return sb.toString();
    }
    
    @Override
    public void increaseKeyValue(final double delta) {
        this._key += delta;
    }
}
