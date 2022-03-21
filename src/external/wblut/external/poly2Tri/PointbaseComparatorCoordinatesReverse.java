// 
// Decompiled by Procyon v0.5.36
// 

package wblut.external.poly2Tri;

import java.util.Comparator;

public class PointbaseComparatorCoordinatesReverse implements Comparator {
	@Override
	public int compare(final Object o1, final Object o2) {
		final Pointbase pb1 = (Pointbase) o1;
		final Pointbase pb2 = (Pointbase) o2;
		return -pb1.compareTo(pb2);
	}
}
