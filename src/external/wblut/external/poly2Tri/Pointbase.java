// 
// Decompiled by Procyon v0.5.36
// 

package wblut.external.poly2Tri;

public class Pointbase implements Comparable
{
    public int id;
    public double x;
    public double y;
    public int type;
    public boolean left;
    
    public Pointbase() {
        this.id = -1;
        this.x = 0.0;
        this.y = 0.0;
        this.type = 1;
        this.left = false;
    }
    
    public Pointbase(final Pointbase pb) {
        this.id = -1;
        this.x = 0.0;
        this.y = 0.0;
        this.type = 1;
        this.left = false;
        this.id = pb.id;
        this.x = pb.x;
        this.y = pb.y;
        this.type = pb.type;
        this.left = pb.left;
    }
    
    public Pointbase(final double xx, final double yy) {
        this.id = -1;
        this.x = 0.0;
        this.y = 0.0;
        this.type = 1;
        this.left = false;
        this.x = xx;
        this.y = yy;
    }
    
    public Pointbase(final int idd, final double xx, final double yy) {
        this.id = -1;
        this.x = 0.0;
        this.y = 0.0;
        this.type = 1;
        this.left = false;
        this.id = idd;
        this.x = xx;
        this.y = yy;
    }
    
    public Pointbase(final double xx, final double yy, final int ttype) {
        this.id = -1;
        this.x = 0.0;
        this.y = 0.0;
        this.type = 1;
        this.left = false;
        this.id = 0;
        this.x = xx;
        this.y = yy;
        this.type = ttype;
    }
    
    public Pointbase(final int idd, final double xx, final double yy, final int ttype) {
        this.id = -1;
        this.x = 0.0;
        this.y = 0.0;
        this.type = 1;
        this.left = false;
        this.id = idd;
        this.x = xx;
        this.y = yy;
        this.type = ttype;
    }
    
    @Override
    public boolean equals(final Object o) {
        return o instanceof Pointbase && this.equals((Pointbase)o);
    }
    
    public boolean equals(final Pointbase pb) {
        return this.x == pb.x && this.y == pb.y;
    }
    
    @Override
    public int compareTo(final Object o) {
        if (!(o instanceof Pointbase)) {
            return -1;
        }
        final Pointbase pb = (Pointbase)o;
        if (this.equals(pb)) {
            return 0;
        }
        if (this.y > pb.y) {
            return 1;
        }
        if (this.y < pb.y) {
            return -1;
        }
        if (this.x < pb.x) {
            return 1;
        }
        return -1;
    }
    
    @Override
    public String toString() {
        return "Pointbase([" + this.x + ", " + this.y + "], ID = " + this.id + ", " + Poly2TriUtils.typeToString(this.type) + ")";
    }
}
