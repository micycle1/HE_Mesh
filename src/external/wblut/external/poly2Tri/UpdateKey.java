// 
// Decompiled by Procyon v0.5.36
// 

package wblut.external.poly2Tri;

import wblut.external.poly2Tri.splayTree.BTreeNode;
import wblut.external.poly2Tri.splayTree.SplayTreeAction;

public class UpdateKey implements SplayTreeAction
{
    @Override
    public void action(final BTreeNode node, final double y) {
        ((Linebase)node.data()).setKeyValue(y);
    }
}
