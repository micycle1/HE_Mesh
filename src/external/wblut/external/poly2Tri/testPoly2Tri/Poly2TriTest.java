// 
// Decompiled by Procyon v0.5.36
// 

package wblut.external.poly2Tri.testPoly2Tri;

import java.awt.Component;
import javax.swing.JOptionPane;
import javax.swing.JFrame;
import java.io.IOException;
import java.io.Reader;
import java.io.LineNumberReader;
import java.io.FileReader;
import wblut.external.poly2Tri.Triangulation;
import java.util.Date;
import java.util.ArrayList;

public class Poly2TriTest
{
    protected static double[] temp;
    protected static int i;
    public static ArrayList<ArrayList<Integer>> triangles;
    public static int numContours;
    public static int[] contours;
    public static double[][] vertices;
    public static int number;
    public static ArrayList<Component> tstfs;
    public static String name;
    
    static {
        Poly2TriTest.triangles = null;
        Poly2TriTest.number = 0;
        Poly2TriTest.tstfs = new ArrayList<Component>();
        Poly2TriTest.name = "";
    }
    
    public static void printVertices(final double[][] vertices) {
        System.out.println("Vertices:");
        Poly2TriTest.i = 1;
        while (Poly2TriTest.i < vertices.length) {
            Poly2TriTest.temp = vertices[Poly2TriTest.i];
            if (Poly2TriTest.i % 6 == 0) {
                System.out.println();
            }
            System.out.print("[" + Poly2TriTest.temp[0] + ", " + Poly2TriTest.temp[1] + "]   ");
            ++Poly2TriTest.i;
        }
        System.out.println();
        if (Poly2TriTest.i % 5 != 0) {
            System.out.println();
        }
    }
    
    public static void doTriangulation() {
        ++Poly2TriTest.number;
        Poly2TriTest.triangles = null;
        final Date begin = new Date();
        Poly2TriTest.triangles = Triangulation.triangulate(Poly2TriTest.numContours, Poly2TriTest.contours, Poly2TriTest.vertices);
        final Date end = new Date();
        final Poly2TriFrame tstf = new Poly2TriFrame();
        final long ms = end.getTime() - begin.getTime();
        tstf.setTitle("Poly2Tri - " + Poly2TriTest.name + " - Time: " + ms + " miliseconds");
        System.out.println(String.valueOf(Poly2TriTest.name) + ", time: " + ms + " miliseconds");
        double maxX = Double.MIN_VALUE;
        double maxY = Double.MIN_VALUE;
        double minX = Double.MAX_VALUE;
        double minY = Double.MAX_VALUE;
        double[] xy1 = { 0.0, 0.0 };
        double[] xy2 = { 0.0, 0.0 };
        double[] xy3 = { 0.0, 0.0 };
        for (int i = 0; i < Poly2TriTest.triangles.size(); ++i) {
            final ArrayList<Integer> t = Poly2TriTest.triangles.get(i);
            for (int j = 0; j < 3; ++j) {
                xy1 = Poly2TriTest.vertices[t.get(j)];
                if (xy1[0] > maxX) {
                    maxX = xy1[0];
                }
                if (xy1[0] < minX) {
                    minX = xy1[0];
                }
                if (xy1[1] > maxY) {
                    maxY = xy1[1];
                }
                if (xy1[1] < minY) {
                    minY = xy1[1];
                }
            }
        }
        tstf.setMaxX(maxX);
        tstf.setMinX(minX);
        tstf.setMaxY(maxY);
        tstf.setMinY(minY);
        for (int i = 0; i < Poly2TriTest.triangles.size(); ++i) {
            final ArrayList<Integer> t = Poly2TriTest.triangles.get(i);
            xy1 = Poly2TriTest.vertices[t.get(0)];
            xy2 = Poly2TriTest.vertices[t.get(1)];
            xy3 = Poly2TriTest.vertices[t.get(2)];
            tstf.addTriangle(xy1[0], xy1[1], xy2[0], xy2[1], xy3[0], xy3[1]);
        }
        tstf.setLocation(Poly2TriTest.number * 30, Poly2TriTest.number * 30);
        tstf.setVisible(true);
        tstf.toFront();
        Poly2TriTest.tstfs.add(tstf);
    }
    
    public static int skipWhitespaces(final String str, int index) {
        while (index < str.length() && (str.charAt(index) == ' ' || str.charAt(index) == '\t')) {
            ++index;
        }
        return index;
    }
    
    public static Double getNumber(final String str, final int[] index) {
        index[0] = skipWhitespaces(str, index[0]);
        if (index[0] >= str.length()) {
            return Double.NaN;
        }
        final StringBuffer temp = new StringBuffer();
        while (index[0] < str.length() && (str.charAt(index[0]) == 'e' || str.charAt(index[0]) == 'E' || str.charAt(index[0]) == '+' || str.charAt(index[0]) == '-' || str.charAt(index[0]) == '.' || (str.charAt(index[0]) >= '0' && str.charAt(index[0]) <= '9'))) {
            temp.append(str.charAt(index[0]));
            final int n = 0;
            ++index[n];
        }
        if (temp.length() == 0) {
            return Double.NaN;
        }
        return new Double(temp.toString());
    }
    
    public static boolean loadBDMFile(final String fileName) {
        final ArrayList<ArrayList<double[]>> contoursAL = new ArrayList<>();
        int numVertices = 0;
        final int[] index = { 0 };
        try {
            final LineNumberReader fr = new LineNumberReader(new FileReader(fileName));
            while (fr.ready()) {
                String line = fr.readLine();
                if (line == null) {
                    continue;
                }
                if (line.length() == 0) {
                    continue;
                }
                if (line.charAt(0) < '0' || line.charAt(0) > '9') {
                    continue;
                }
                Double num = getNumber(line, new int[1]);
                final int curNumVertices = (int)Math.round(num);
                numVertices += curNumVertices;
                final ArrayList<double[]> verticesAL = new ArrayList<>();
                for (int i = 0; i < curNumVertices; ++i) {
                    line = fr.readLine();
                    if (line == null) {
                        return false;
                    }
                    index[0] = 0;
                    num = getNumber(line, index);
                    if (num == Double.NaN) {
                        return false;
                    }
                    final Double num2 = getNumber(line, index);
                    if (num2 == Double.NaN) {
                        return false;
                    }
                    verticesAL.add(new double[] { num, num2 });
                }
                contoursAL.add(verticesAL);
            }
            fr.close();
        }
        catch (IOException e) {
            return false;
        }
        Poly2TriTest.numContours = contoursAL.size();
        Poly2TriTest.contours = new int[Poly2TriTest.numContours];
        Poly2TriTest.vertices = new double[numVertices][2];
        int k = 0;
        for (int j = 0; j < Poly2TriTest.numContours; ++j) {
            final ArrayList<double[]> verticesAL = contoursAL.get(j);
            Poly2TriTest.contours[j] = verticesAL.size();
            for (int l = 0; l < verticesAL.size(); ++l) {
                Poly2TriTest.vertices[k++] = verticesAL.get(l);
            }
        }
        return true;
    }
    
    public static void main(final String[] args) {
        Poly2TriTest.numContours = 1;
        Poly2TriTest.contours = new int[] { 5 };
        Poly2TriTest.vertices = new double[][] { { 1.0, 1.0 }, { 4.0, 0.0 }, { 3.0, 5.0 }, { 2.0, 8.0 }, { 0.0, 0.0 } };
        Poly2TriTest.name = "1st";
        doTriangulation();
        Poly2TriTest.numContours = 2;
        Poly2TriTest.contours = new int[] { 3, 3 };
        Poly2TriTest.vertices = new double[][] { { 0.0, 0.0 }, { 7.0, 0.0 }, { 3.0, 4.0 }, { 2.0, 1.0 }, { 2.0, 2.0 }, { 3.0, 1.0 } };
        Poly2TriTest.name = "2nd";
        doTriangulation();
        Poly2TriTest.numContours = 3;
        Poly2TriTest.contours = new int[] { 3, 3, 3 };
        Poly2TriTest.vertices = new double[][] { { 0.0, 0.0 }, { 9.0, 0.0 }, { 3.0, 8.0 }, { 2.0, 1.0 }, { 2.0, 2.0 }, { 3.0, 1.0 }, { 6.0, 2.0 }, { 6.0, 1.0 }, { 5.0, 2.0 } };
        Poly2TriTest.name = "3rd";
        doTriangulation();
        final String prefix = "c:/temp/";
        final String[] bdmFiles = { "boxc100.bdm", "circle1.bdm", "circle2.bdm", "circle3.bdm", "crazybox1.bdm", "crazybox2.bdm", "guitar.bdm", "sample1.bdm", "sample2.bdm", "sample3.bdm", "superior.bdm" };
        boolean end = false;
        boolean next = false;
        boolean all = false;
        int option = JOptionPane.showConfirmDialog(Poly2TriTest.tstfs.get(Poly2TriTest.tstfs.size() - 1), "Do you want to do ALL bdm files triangulation?");
        all = (option == 0);
        for (int i = 0; i < bdmFiles.length; ++i) {
            if (!all) {
                next = false;
                option = JOptionPane.showConfirmDialog(Poly2TriTest.tstfs.get(Poly2TriTest.tstfs.size() - 1), "Do triangulation of " + bdmFiles[i] + "?");
                switch (option) {
                    case 2: {
                        end = true;
                        break;
                    }
                    case 1: {
                        next = true;
                        break;
                    }
                }
                if (end) {
                    break;
                }
                if (next) {
                    continue;
                }
            }
            if (loadBDMFile(String.valueOf(prefix) + bdmFiles[i])) {
                Poly2TriTest.name = bdmFiles[i];
                doTriangulation();
            }
            else {
                JOptionPane.showMessageDialog(Poly2TriTest.tstfs.get(Poly2TriTest.tstfs.size() - 1), "Failed to load " + bdmFiles[i]);
            }
        }
    }
}
