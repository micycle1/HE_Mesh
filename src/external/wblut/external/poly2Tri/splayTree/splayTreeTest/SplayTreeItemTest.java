// 
// Decompiled by Procyon v0.5.36
// 

package wblut.external.poly2Tri.splayTree.splayTreeTest;

import wblut.external.poly2Tri.splayTree.SplayTreeItem;

public class SplayTreeItemTest implements SplayTreeItem
{
    public int data;
    
    public SplayTreeItemTest(final int a) {
        this.data = 1;
        this.data = a;
    }
    
    @Override
    public void increaseKeyValue(final double delta) {
        ++this.data;
    }
    
    @Override
    public Comparable keyValue() {
        return this.data;
    }
}
