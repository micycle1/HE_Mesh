// 
// Decompiled by Procyon v0.5.36
// 

package wblut.external.poly2Tri.splayTree.splayTreeTest;

import wblut.external.poly2Tri.splayTree.BTreeNode;
import wblut.external.poly2Tri.splayTree.SplayTree;

public class SplayTreeTest {
	public static void main(final String[] argv) {
		final SplayTree t = new SplayTree();
		final int NUMS = 40000;
		final int GAP = 307;
		System.out.println("Checking... (no bad output means success)");
		for (int i = 307; i != 0; i = (i + 307) % 40000) {
			t.insert(new SplayTreeItemTest(i));
		}
		System.out.println("Inserts complete");
		for (int i = 1; i < 40000; i += 2) {
			t.delete(i);
		}
		System.out.println("Removes complete");
		for (int j = 2; j < 40000; j += 2) {
			final BTreeNode node = t.find(j);
			if (node == null || ((SplayTreeItemTest) node.data()).data != j) {
				System.out.println("Error: find fails for " + j);
			}
		}
		for (int j = 1; j < 40000; j += 2) {
			if (t.find(j) != null) {
				System.out.println("Error: Found deleted item " + j);
			}
		}
		for (int j = 4; j < 40000; j += 2) {
			final BTreeNode node = t.findMaxSmallerThan(j);
			if (((SplayTreeItemTest) node.data()).data != j - 2) {
				System.out.println("Error: findMaxSmallerThan(i) == " + ((SplayTreeItemTest) node.data()).data + " != " + (j - 2));
			}
		}
		System.out.println("FindMaxSmallerThen complete");
	}
}
