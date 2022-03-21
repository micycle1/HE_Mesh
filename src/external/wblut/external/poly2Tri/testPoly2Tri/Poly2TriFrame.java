// 
// Decompiled by Procyon v0.5.36
// 

package wblut.external.poly2Tri.testPoly2Tri;

import java.awt.LayoutManager;
import java.awt.Container;
import javax.swing.JPanel;
import javax.swing.JFrame;

public class Poly2TriFrame extends JFrame
{
    private static final long serialVersionUID = 1L;
    private JPanel jContentPane;
    private double maxX;
    private double maxY;
    private double minX;
    private double minY;
    
    public Poly2TriFrame() {
        this.jContentPane = null;
        this.maxX = 10.0;
        this.maxY = 10.0;
        this.minX = -10.0;
        this.minY = -10.0;
        this.initialize();
    }
    
    private void initialize() {
        this.setSize(300, 300);
        this.setContentPane(this.getJContentPane());
        this.setTitle("Poly2Tri");
        this.setDefaultCloseOperation(3);
    }
    
    private JPanel getJContentPane() {
        if (this.jContentPane == null) {
            (this.jContentPane = new Poly2TriPainting()).setLayout(null);
        }
        return this.jContentPane;
    }
    
    public void setMaxX(final double newMaxX) {
        if (newMaxX > 0.0) {
            this.maxX = newMaxX + 0.5;
        }
        else {
            this.maxX = newMaxX + 0.5;
        }
        ((Poly2TriPainting)this.getContentPane()).maxX = this.maxX;
    }
    
    public void setMaxY(final double newMaxY) {
        if (newMaxY > 0.0) {
            this.maxY = newMaxY + 0.5;
        }
        else {
            this.maxY = newMaxY + 0.5;
        }
        ((Poly2TriPainting)this.getContentPane()).maxY = this.maxY;
    }
    
    public void setMinX(final double newMinX) {
        if (newMinX > 0.0) {
            this.minX = newMinX - 0.5;
        }
        else {
            this.minX = newMinX - 0.5;
        }
        ((Poly2TriPainting)this.getContentPane()).minX = this.minX;
    }
    
    public void setMinY(final double newMinY) {
        if (newMinY > 0.0) {
            this.minY = newMinY - 0.5;
        }
        else {
            this.minY = newMinY - 0.5;
        }
        ((Poly2TriPainting)this.getContentPane()).minY = this.minY;
    }
    
    public void addTriangle(final double x1, final double y1, final double x2, final double y2, final double x3, final double y3) {
        ((Poly2TriPainting)this.getContentPane()).addPolygon(new double[] { x1, y1, x2, y2, x3, y3 });
    }
}
