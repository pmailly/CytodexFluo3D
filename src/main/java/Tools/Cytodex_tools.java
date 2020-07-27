/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Tools;

import CytoDex.Cytodex_Fluo3D;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.EllipseRoi;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.gui.WaitForUserDialog;
import ij.io.FileSaver;
import ij.io.Opener;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.ZProjector;
import ij.plugin.filter.RankFilters;
import ij.process.AutoThresholder;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import ij.util.ArrayUtil;
import ij.util.Tools;
import java.awt.Polygon;
import java.awt.image.IndexColorModel;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import static java.lang.Math.abs;
import static java.lang.Math.round;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import mcib3d.geom.ObjectCreator3D;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Point3D;
import mcib3d.geom.Voxel3D;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;
import sc.fiji.analyzeSkeleton.AnalyzeSkeleton_;
import sc.fiji.analyzeSkeleton.Edge;
import sc.fiji.analyzeSkeleton.Graph;
import sc.fiji.analyzeSkeleton.Point;
import sc.fiji.analyzeSkeleton.SkeletonResult;
import sc.fiji.localThickness.Clean_Up_Local_Thickness;
import sc.fiji.localThickness.EDT_S1D;
import sc.fiji.localThickness.Local_Thickness_Parallel;

/**
 *
 * @author phm
 */
public class Cytodex_tools {
    

    public static final boolean canceled = false;
    public static int nbNucleus;   // total number of nucleus
    public static String imgOutDir;
    public static String fileNameWithOutExt;
    
    public static final int shollStep = 10;
    public static double cytoDexRad;
    public static int cytoDexCenterZ;
    public static int cytoDexCenterX;
    public static int cytoDexCenterY;
    public static double cytoDexEnlarge = 0;    // thickness to enlarge cytoDexRad
    public static double minThreshold_Nucleus = 0;
    public static double maxThreshold_Nucleus = 0;
    private static double minThreshold_Vessel = 0;
    private static double maxThreshold_Vessel = 0;
    
    
    /**
     * Returns an IndexColorModel similar to MATLAB's jet color map.
     *
     * @param backgroundGray
     *            the gray value (8-bit scale) to be used as the first entry of
     *            the LUT. It is ignored if negative.
     * @return The "Jet" LUT with the specified background entry
     * @see <a href=
     *      "https://list.nih.gov/cgi-bin/wa.exe?A2=IMAGEJ;c8cb4d8d.1306">Jerome
     *      Mutterer alternative</a>
     */
    public static IndexColorModel matlabJetColorMap(final int backgroundGray, final int foregroundGray) {

            // Initialize colors arrays (zero-filled by default)
            final byte[] reds = new byte[256];
            final byte[] greens = new byte[256];
            final byte[] blues = new byte[256];

            // Set greens, index 0-32; 224-255: 0
            for (int i = 0; i < 256 / 4; i++) // index 32-96
                    greens[i + 256 / 8] = (byte) (i * 255 * 4 / 256);
            for (int i = 256 * 3 / 8; i < 256 * 5 / 8; ++i) // index 96-160
                    greens[i] = (byte) 255;
            for (int i = 0; i < 256 / 4; i++) // index 160-224
                    greens[i + 256 * 5 / 8] = (byte) (255 - (i * 255 * 4 / 256));

            // Set blues, index 224-255: 0
            for (int i = 0; i < 256 * 7 / 8; i++) // index 0-224
                    blues[i] = greens[(i + 256 / 4) % 256];

            // Set reds, index 0-32: 0
            for (int i = 256 / 8; i < 256; i++) // index 32-255
                    reds[i] = greens[(i + 256 * 6 / 8) % 256];

            // Set background and foreground colors
            if (backgroundGray >= 0) {
                    reds[0] = greens[0] = blues[0] = (byte) backgroundGray;
            }
            if (foregroundGray >= 0) {
                    reds[255] = greens[255] = blues[255] = (byte) foregroundGray;
            }

            return new IndexColorModel(8, 256, reds, greens, blues);

    }
    
    /**
     * Find nucleus distance to cytodex
     * @param x nucleus positions
     * @param y nucleus positions
     * @param points roi
     * @return 
     */
    public double getClosestDistance(int x, int y, Polygon points) {
        double distance = Double.MAX_VALUE;
        for (int i = 0; i < points.npoints; i++) {
                double dx = points.xpoints[i] - x;
                double dy = points.ypoints[i] - y;
                double distance2 = Math.sqrt(dx*dx+dy*dy);
                if (distance2 < distance) {
                        distance = distance2;
                }
        }
        return distance;
    }
    
    /**
     * Find Z with max intensity in stack
     * @param img
     * @return z
     */
    
    private static int find_max(ImagePlus img) {
        double max = 0;
        int zmax = 0;
        for (int z = 1; z <= img.getNSlices(); z++) {
            ImageProcessor ip = img.getStack().getProcessor(z);
            ImageStatistics statistics = new ImageStatistics().getStatistics(ip, ImageStatistics.MEAN, img.getCalibration());
            double meanInt = statistics.mean;
            if (meanInt > max) {
                max = meanInt;
                zmax = z;
            }
        }
        return(zmax);
    }
    
     /*Median filter 
     * 
     * @param img
     * @param size
     */ 
    public static void median_filter(ImagePlus img, double size) {
        RankFilters median = new RankFilters();
        for (int s = 1; s <= img.getNSlices(); s++) {
            img.setZ(s);
            median.rank(img.getProcessor(), size, RankFilters.MEDIAN);
            img.updateAndDraw();
        }
    }
    
        // Threshold images and fill holes
    public static void threshold(ImagePlus img, AutoThresholder.Method thMed, boolean fill, boolean calcul) {
        //  Threshold and binarize
       String cal = "";
       img.setZ(img.getNSlices()/2);
       img.updateAndDraw();
       IJ.setAutoThreshold(img, thMed.toString()+" dark");
       Prefs.blackBackground = false;
       if (calcul)
           cal = " calculate";
        IJ.run(img, "Convert to Mask", "method="+thMed.toString()+cal+" background=Dark");
        if (fill) {
            IJ.run(img,"Fill Holes", "stack");
        }
    }
    
    /**
    * Detect bead center sphere position, compute radius
     * @param img
     * @param roi
    */
    public static void findSpheroid(ImagePlus img, Roi roi) {
        // Find center position of cytodex measuring max area of detected object
        IJ.showStatus("Finding bead");
        ImagePlus imgDup = img.duplicate();
        Calibration cal = imgDup.getCalibration();
        imgDup.setCalibration(cal);
        IJ.run(imgDup, "Gaussian Blur...", "radius=4 stack");
        if (roi != null) {
            imgDup.setRoi(roi);
            imgDup.updateAndDraw();
            IJ.run("Colors...", "foreground=white background=black selection=yellow");
            IJ.run(imgDup, "Clear Outside","stack");
            imgDup.deleteRoi();
        }
        // Find bead diameter on Z projection
        ImagePlus imgProj = doZProjection(imgDup);
        imgProj.setCalibration(cal);
        IJ.setAutoThreshold(imgProj, "Otsu dark");
        Prefs.blackBackground = false;
        IJ.run(imgProj, "Convert to Mask", "method=Otsu background=Dark");
        IJ.run(imgProj,"Fill Holes", "");
        IJ.run("Set Measurements...","centroid feret's stack redirect=None decimal=4");
        IJ.run(imgProj,"Analyze Particles...","size=0-infinity circularity=0.7-1.00 include exclude clear");
        ResultsTable rt = ResultsTable.getResultsTable();
        double centerX = 0;
        double centerY = 0;
        double centerZ = 0;
        double feretProj = 0; 
        imgProj.close();
        // no object
        if (rt.size() == 0) {
            System.out.println("No bead");
        }
        // one object
        else {
            feretProj = rt.getValue("Feret", 0);
            centerX = rt.getValue("X", 0);
            centerY = rt.getValue("Y", 0);
            rt.reset();
            // find center position
            imgDup.setSlice(find_max(imgDup));
            imgDup.updateAndDraw();
            IJ.setAutoThreshold(imgDup, "Otsu dark");
            IJ.run(imgDup, "Convert to Mask", "method=Otsu background=Dark");
            IJ.run(imgDup,"Fill Holes", "stack");
            IJ.run(imgDup, "Options...", "iterations=6 count=1 do=Close stack");
            IJ.run(imgDup,"Analyze Particles...","size=500-infinity circularity=0.7-1.00 include exclude clear stack");
            // find the nearest feret value from feretProj
            double[] ferets = new double[rt.size()];
            for (int r = 0; r < rt.size(); r++) {
                ferets[r] = rt.getValue("Feret", r);
            }
            double diff = Math.abs(ferets[0] - feretProj);
            int idx = 0;
            for(int c = 1; c < ferets.length; c++){
                double cdiff = Math.abs(ferets[c] - feretProj);
                if(cdiff < diff){
                    idx = c;
                    diff = cdiff;
                }
            }
            centerZ = rt.getValue("Slice",idx);
            cytoDexCenterX = (int)(centerX/cal.pixelWidth);
            cytoDexCenterY = (int)(centerY/cal.pixelWidth);
            cytoDexCenterZ = (int)centerZ;
            cytoDexRad = feretProj/2;
            System.out.println("x = "+cytoDexCenterX+" Y= "+cytoDexCenterY+" Z= "+cytoDexCenterZ+" Rad = "+cytoDexRad);
        }
        imgDup.close();
    } 
    
    
    /*
    Crop Z stack
    */
    public static void cropStack(ImagePlus img, String dir) {
        Calibration cal = img.getCalibration();
        double cropZ = round(abs((cytoDexCenterZ - cytoDexRad/cal.pixelDepth) - 1));
        //System.out.println("removing slices from 1 to "+cropZ+ " stack size = "+img.getNSlices());
        for (int z = 1; z <= cropZ; z++) {
            if ("Top".equals(dir))
                img.getStack().deleteSlice(1);
            else
                img.getStack().deleteLastSlice();
        }
    }
    
    
    // Find vessels with DOG
    public static ImagePlus vesselDOG(ImagePlus img, String seriesName, Roi roi, int first) {
        
        IJ.run(img, "Median...", "radius=2 stack");
        IJ.run(img, "Difference of Gaussians", "  sigma1=15 sigma2=10 stack");
        img.setSlice(cytoDexCenterZ);
        img.updateAndDraw();
        if (first == 1) {
            img.getProcessor().resetThreshold();
            img.show();
            IJ.setAutoThreshold(img, "Default dark stack");
            IJ.run("Threshold...");
            if (!WindowManager.getWindow("Threshold").isVisible());
            WindowManager.getWindow("Threshold").setVisible(true);
            new WaitForUserDialog("Set threshold values using set button and press OK").show();
            minThreshold_Vessel = img.getProcessor().getMinThreshold();
            maxThreshold_Vessel = img.getProcessor().getMaxThreshold();
            WindowManager.getWindow("Threshold").setVisible(false);
            EllipseRoi roi_Circle = new EllipseRoi(cytoDexCenterX - cytoDexRad, cytoDexCenterY - cytoDexRad, 
                    cytoDexCenterX + cytoDexRad, cytoDexCenterY + cytoDexRad, 1);
            IJ.setTool("oval");
            img.setRoi(roi_Circle);
            img.updateAndDraw();
            new WaitForUserDialog("Resize circle to define bead diameter and press OK").show();
            img.updateAndDraw();
            double cytoDexRadNew = img.getRoi().getFeretsDiameter()/2;
            cytoDexEnlarge = cytoDexRadNew - cytoDexRad;
            System.out.println("New cyto radius = "+cytoDexRadNew+ " enlarge thickness = "+cytoDexEnlarge);
            img.deleteRoi();   
            img.hide();  
        }
        System.out.println("MinTh = "+minThreshold_Vessel+" MaxTh = "+maxThreshold_Vessel);
        Prefs.blackBackground = false;
        IJ.setThreshold(img, minThreshold_Vessel, maxThreshold_Vessel);
        IJ.run(img, "Convert to Mask","background=Dark black");
        IJ.run(img, "Options...", "iterations=5 count=1 do=Close stack");
        removeSpheroid(img, 0);
        if (roi != null) {
            img.setRoi(roi);
            IJ.run("Colors...", "foreground=white background=black selection=yellow");
            IJ.run(img, "Clear Outside","stack");
            img.deleteRoi();
        }
        FileSaver imgObjectsSave = new FileSaver(img);
        imgObjectsSave.saveAsTiff(imgOutDir+fileNameWithOutExt+"_VesselMask.tif");
        return(img);
    }
    
    /*
    Histogramm frequency
    */
    public static int[] calcHistogram(float[] data, double min, double max, int binSize) {
        int numBins = (int) Math.round((max - min)/binSize);
        final int[] result = new int[numBins];

        for (double d : data) {
          int bin = (int) ((d - min) / binSize);
          if (bin < 0) { /* this data is smaller than min */ }
          else if (bin >= numBins) { /* this data point is bigger than max */ }
          else {
            result[bin] += 1;
          }
        }
        return result;
    }

    public static double calculateDistance(Point point1, Point point2) {
        return Math.sqrt( Math.pow(point1.x - point2.x, 2) 
                          + Math.pow(point1.y - point2.y, 2)
                          + Math.pow(point1.z - point2.z, 2));
    }
    
    /**
     * Returns the location of pixels clockwise along circumference
     * using Bresenham's Circle Algorithm
     * keep point only if pixel is inside image
     */
    public static ArrayList<Point> BresenhamCircle(int xc,int yc,int r) {    
        ArrayList<Point> ret = new ArrayList<>();
        int x,y,p;
        x=0;
        y=r;
        ret.add(new Point(xc+x,yc-y,0));
        p=3-(2*r);
        for(x=0;x<=y;x++) {
            if (p<0) {
                p=(p+(4*x)+6);
            }
            else {
                y=y-1;
                p=p+((4*(x-y)+10));
            }
            ret.add(new Point(xc+x,yc-y,0));
            ret.add(new Point(xc-x,yc-y,0));
            ret.add(new Point(xc+x,yc+y,0));
            ret.add(new Point(xc-x,yc+y,0));
            ret.add(new Point(xc+y,yc-x,0));
            ret.add(new Point(xc-y,yc-x,0));
            ret.add(new Point(xc+y,yc+x,0));
            ret.add(new Point(xc-y,yc+x,0));
        }
        return ret;
    }
    
    /**
     * Draw the capillary diameter value inside two circles on the skeleton image
     */
    public static void drawInsideCircle(float rad1, float rad2, double value, ImagePlus img) {
        for(int i = (int)rad1; i < rad2; i++) {
            ArrayList<Point> points = BresenhamCircle(cytoDexCenterX, cytoDexCenterY,i);
            for(int n = 0; n < points.size(); n++) {
                double pixelValue = img.getProcessor().getPixel(points.get(n).x, points.get(n).y);
                if (pixelValue == 255) {
                    img.getProcessor().putPixel(points.get(n).x, points.get(n).y, (int)value);
                }
            }
        }
    }
    
    
    /*
    * clear inside cytodex
    */
    public static void removeSpheroid(ImagePlus img, int color) {
        IJ.showStatus("Removing bead ...");
        Calibration cal = img.getCalibration();
        ObjectCreator3D sphereC0 = new ObjectCreator3D(img.getStack());
        double sphereRad = (cytoDexRad  + cytoDexEnlarge) / cal.pixelWidth;
        //System.out.println("sphere rad = " + (cytoDexRad  + cytoDexEnlarge));
        sphereC0.createEllipsoid(cytoDexCenterX, cytoDexCenterY, cytoDexCenterZ, sphereRad, sphereRad, sphereRad/cal.pixelDepth, color, false);
        img.updateAndDraw();
    } 
    
/**
     * compute local thickness
     * @param imgMask
     * @return imgMap
    **/
    public static ImagePlus localThickness3D (ImagePlus imgMask) {
        Calibration cal = imgMask.getCalibration();
        ImagePlus img = imgMask.duplicate();
        img.setCalibration(cal);
        EDT_S1D edt = new EDT_S1D();
        edt.runSilent = true;
        edt.thresh = 1;
        edt.inverse = false;
        edt.showOptions = false;
        edt.setup("", img);
        edt.run(img.getProcessor());
        ImagePlus imgEDT = edt.getResultImage();
        imgEDT.setCalibration(cal);
        Local_Thickness_Parallel locThk = new Local_Thickness_Parallel();
        locThk.runSilent = true;
        locThk.setup("", imgEDT);
        locThk.run(imgEDT.getProcessor());
        ImagePlus imgLocThk = locThk.getResultImage();
        imgLocThk.setCalibration(cal);
        Clean_Up_Local_Thickness cleanUp = new Clean_Up_Local_Thickness();
        cleanUp.runSilent = true;
        cleanUp.setup("", imgLocThk);
        cleanUp.run(imgLocThk.getProcessor());
        ImagePlus imgMap = cleanUp.getResultImage();
        // calibrate intensity to µm
        cal.setValueUnit("µm");
        imgMap.setCalibration(cal);
        IJ.run(imgMap, "Multiply...", "value="+cal.pixelWidth+" stack");
        img.close();
        imgEDT.close();
        imgLocThk.close();
        return(imgMap);
    }

    /*
    * Find objects in image
    */
    public static Objects3DPopulation getPopFromImage(ImagePlus img) {
        // label binary images first
        ImageLabeller labeller = new ImageLabeller();
        ImageInt labels = labeller.getLabels(ImageHandler.wrap(img));
        Objects3DPopulation pop = new Objects3DPopulation(labels);
        return pop;
    }
    

    /* count nucleus population
    /* 
    */
    public static Objects3DPopulation findNucleus(ImagePlus imgNuc, Roi roi, int first) throws IOException {
        Calibration cal = imgNuc.getCalibration();
        String imgTitle = imgNuc.getTitle();
        IJ.run(imgNuc, "Difference of Gaussians", "  sigma1=10 sigma2=5 stack");
    // threshold image                        
        imgNuc.setSlice(cytoDexCenterZ);
        imgNuc.updateAndDraw();        
        IJ.setAutoThreshold(imgNuc, "Moments dark stack");
        Prefs.blackBackground = false;
        IJ.run(imgNuc, "Convert to Mask", "method=Moments background=Dark black");
        if (roi != null) {
            imgNuc.setRoi(roi);
            IJ.run("Colors...", "foreground=white background=black selection=yellow");
            IJ.run(imgNuc, "Clear Outside","stack");
            imgNuc.deleteRoi();
        }
        removeSpheroid(imgNuc, 0);
//        FileSaver imgObjectsSave = new FileSaver(imgNuc);
//        imgObjectsSave.saveAsTiff(imgOutDir+fileNameWithOutExt+"_NucleusMask.tif");
        IJ.showStatus("Computing nucleus population ...");        
        Objects3DPopulation nucleus = getPopFromImage(imgNuc);
        if (nucleus == null)
            nbNucleus = 0;
        else {
            nbNucleus = nucleus.getNbObjects();
            float[] nucleusDist = new float[nbNucleus];
        // Save tagged image   
            ImageInt img = ImageInt.wrap(imgNuc);
            ImageHandler imgObjects = img.createSameDimensions();
            imgObjects.set332RGBLut();
            nucleus.draw(imgObjects, 30);
            imgObjects.setCalibration(cal);
            FileSaver imgSave = new FileSaver(imgObjects.getImagePlus());
            imgSave.saveAsTiff(imgOutDir+imgTitle+".tif");
            imgObjects.flush();
            imgObjects.closeImagePlus();
        // Compute distances  and save data  
            FileWriter fwNucleusDistances = new FileWriter(imgOutDir + imgTitle+"_distance.xls",false);
            BufferedWriter outputNucleusDistances = new BufferedWriter(fwNucleusDistances);
            outputNucleusDistances.flush();
            IJ.showStatus("Computing distance to Cytodex center");
            for (int i = 0; i < nucleus.getNbObjects(); i++) {
                nucleusDist[i] = (float) Math.abs(nucleus.getObject(i).distPixelCenter(cytoDexCenterX, cytoDexCenterY, cytoDexCenterZ));
            }
            ArrayUtil au = new ArrayUtil(nucleusDist);
            int distFreq[] = calcHistogram(nucleusDist, au.getMinimum(), au.getMaximum(), shollStep);
            //Save distances
            double min = au.getMinimum();
            outputNucleusDistances.write("Distance\t#Nucleus\n");
            for (int i = 0; i < distFreq.length; i++) {    
                outputNucleusDistances.write(min + "\t" + distFreq[i] + "\n");
                outputNucleusDistances.flush();
                min += shollStep;
            }
            outputNucleusDistances.close();
            imgNuc.changes = false;
            imgNuc.close();
            imgNuc.flush();
        }
        return(nucleus);
    }
    
    
    /**
	 * Prune end branches of a specific length
	 *
	 * @param stack input skeleton image
	 * @param results
	 * @param length limit length to prune the branches (in calibrated units)
	 *
	 */
    public static ImagePlus pruneEndBranches(ImageStack stack, SkeletonResult results, double length) {
        Graph[] graph = results.getGraph();
        ArrayList<Point> endPoints = results.getListOfEndPoints();
        for (int i = 0; i < graph.length; i++) {
            ArrayList<Edge> listEdges = graph[i].getEdges();
 
            // go through all branches and remove branches under threshold
            // in duplicate image
            for (Edge e : listEdges) {
                ArrayList<Point> p = e.getV1().getPoints();
                boolean v1End = endPoints.contains(p.get(0));
                ArrayList<Point> p2 = e.getV2().getPoints();
                boolean v2End = endPoints.contains(p2.get(0));
                // if any of the vertices is end-point 
                if (v1End || v2End) {
                    if (e.getLength() < length) {
                        if (v1End)
                            stack.setVoxel(p.get(0).x, p.get(0).y, p.get(0).z, 0);
                        if (v2End)
                            stack.setVoxel(p2.get(0).x, p2.get(0).y, p2.get(0).z, 0);
                        for (Point pt : e.getSlabs())
                            stack.setVoxel(pt.x, pt.y, pt.z, 0);
                    }
                }
            }
        }
        return(new ImagePlus("", stack));
    }



        
      
    // Calculate lenght of branches after skeletonize
    public static void analyzeSkel (ImagePlus img, BufferedWriter output, String outDir, double smallBranch, int iterPruning) {
	Calibration cal = img.getCalibration();
        String imgTitle = img.getTitle();
        AnalyzeSkeleton_ analyzeSkeleton = new AnalyzeSkeleton_();
        AnalyzeSkeleton_.calculateShortestPath = true;
        analyzeSkeleton.setup("",img);
        SkeletonResult skeletonResults = analyzeSkeleton.run(AnalyzeSkeleton_.NONE, false, true, null, true, true);
        // remove small branches
        IJ.showStatus("Removing small branches...");
        for (int i = 0; i < iterPruning; i++) {
            IJ.showStatus("Removing small branches step "+(i+1)+"/"+iterPruning);
            ImagePlus imgPrune = pruneEndBranches(img.getStack(), skeletonResults, smallBranch);
            analyzeSkeleton.setup("",imgPrune);
            skeletonResults = analyzeSkeleton.run(AnalyzeSkeleton_.NONE, false, false, null, true, false);
            //System.out.println("length = "+ skeletonResults.getBranches().length);
        }
        
        //  compute parameters for each skeleton
        ImageStack imgStackLab = analyzeSkeleton.getLabeledSkeletons();
        IJ.showStatus("Computing parameters for each skeleton ...");
        ImagePlus imgLab = new ImagePlus(imgTitle+"_LabelledSkel.tif", imgStackLab);
        ImagePlus imgLabProj = doZProjection(imgLab);
        IJ.run(imgLabProj, "3-3-2 RGB", "");
        imgLabProj.setCalibration(cal);
        imgLab.changes = false;
        imgLab.close();
        imgLab.flush();   
        int[] branchNumbers = skeletonResults.getBranches();
        double[] branchLengths = skeletonResults.getAverageBranchLength();
        double[] totalLengths = new double[branchLengths.length];
        for (int b = 0; b <  branchLengths.length; b++)
            totalLengths[b] += branchNumbers[b] * branchLengths[b];
        int[] nbEndPoints = skeletonResults.getEndPoints();
        int[] junctions = skeletonResults.getJunctions();
        for (int i = 0; i < skeletonResults.getGraph().length; i++) {
            try {
                if (branchNumbers[i] != 0) {
                    // write data
                    if (i == 0)
                        output.write(imgTitle + "\t" + (i+1) + "\t" + branchNumbers[i] + "\t" + junctions[i] + "\t" + nbEndPoints[i] + "\t" + totalLengths[i] + "\n");
                    else
                        output.write("\t" + (i+1) + "\t" + branchNumbers[i] + "\t" + junctions[i] + "\t" + nbEndPoints[i] + "\t" + totalLengths[i] + "\n");
                    output.flush();
                }
            } catch (IOException ex) {
                Logger.getLogger(Cytodex_Fluo3D.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

        FileSaver imgSave = new FileSaver(imgLabProj);
        imgSave.saveAsTiff(outDir+imgTitle.substring(0,imgTitle.indexOf(".tif"))+"_LabelledSkel.tif");
        imgLabProj.changes = false;
        imgLabProj.close();  
        imgLabProj.flush();  
    }
   


    public static ImagePlus doZProjection(ImagePlus img) {
        ZProjector zproject = new ZProjector();
        zproject.setMethod(ZProjector.MAX_METHOD);
        zproject.setStartSlice(1);
        zproject.setStopSlice(img.getNSlices());
        zproject.setImage(img);
        zproject.doProjection();
       return(zproject.getProjection());
    }
    
    public static void intersectionAnalysis(ImagePlus img) throws IOException {
        Calibration cal = img.getCalibration();
        String imgTitle = img.getTitle();
        final int n_cpus = Prefs.getThreads();
        int dx = (cytoDexCenterX <= img.getWidth()/2) ? (cytoDexCenterX-img.getWidth()) : cytoDexCenterX;
        int dy = (cytoDexCenterY <= img.getHeight()/2) ? (cytoDexCenterY-img.getHeight()) : cytoDexCenterY;
        int dz = (cytoDexCenterZ <= img.getNSlices()/2) ? (cytoDexCenterZ-img.getNSlices()) : cytoDexCenterZ;
        int maxEndRadius = (int) Math.sqrt(dx*dx + dy*dy + dz*dz);
        PointRoi shollCenter = new PointRoi(cytoDexCenterX, cytoDexCenterY);
        img.setZ(cytoDexCenterZ);
        img.setRoi(shollCenter, true);
        img.show();
        System.out.println(cytoDexRad*cal.pixelWidth + ", "+maxEndRadius+", "+shollStep);
        IJ.run("Metrics & Options...", "image spatial threshold center starting radius samples enclosing intersecting sum mean append=[] plots=[Only linear profile]"
                + " background=0 preferred=[Sampled values] parallel="+n_cpus+" file=.xls decimal=4");
        IJ.run(img, "Sholl Analysis...", "starting="+cytoDexRad*cal.pixelWidth+" ending="+maxEndRadius+" radius_step="+shollStep+" ignore enclosing=0 #_primary=10"+
                " fit linear polynomial=[Best fitting degree] most normalizer=Volume create overlay save directory="+imgOutDir+" do");    

        img.deleteRoi();
        //Add scale bar and calibration bar to saved sholl image
        ImagePlus imgSholl = new Opener().openImage(imgOutDir, imgTitle.substring(0,imgTitle.indexOf(".tif"))+"_ShollMask.tif");
        // Check if image map is saved may be problem occurs on Mac and PC ???
        if (imgSholl == null) {
            imgSholl = WindowManager.getCurrentImage();
        }
        IJ.run(imgSholl,  "Calibration Bar...", "location=[Upper Right] fill=White label=Black number=5 decimal=0 font=12 zoom=2 overlay");
        IJ.run(imgSholl, "Scale Bar...", "width=100 height=10 font=42 color=Yellow background=None location=[Lower Right]");
        FileSaver imgSave = new FileSaver(imgSholl);
        imgSave.saveAsTiff(imgOutDir+imgTitle.substring(0,imgTitle.indexOf(".tif"))+"_Skel_ShollMask.tif");
        imgSholl.changes = false;
        imgSholl.close();
        imgSholl.flush();
        img.hide();
        //IJ.run("Close");
    }
   

    // For each skeleton calculate mean diameter inside sphere from cytodex center using image map
    // compute number of intersections
    public static void diameterAnalysis (ImagePlus imgSkel, ImageHandler imgMap) throws IOException {
        Calibration cal = imgSkel.getCalibration();
        String imgTitle = imgSkel.getTitle();
        // parameters
        final int wdth = imgSkel.getWidth();
        final int hght = imgSkel.getHeight();
        final int depth = imgSkel.getNSlices();
        final double dx, dy, dz, maxEndRadius;
        final double vxWH = Math.sqrt(cal.pixelWidth * cal.pixelHeight);
	final double vxD = cal.pixelDepth;
        dx = (cytoDexCenterX <= wdth / 2) ? (cytoDexCenterX - wdth) * vxWH : cytoDexCenterX * vxWH;
        dy = (cytoDexCenterY <= hght / 2) ? (cytoDexCenterY - hght) * vxWH : cytoDexCenterY * vxWH;
        dz = (cytoDexCenterZ <= depth / 2) ? (cytoDexCenterZ - depth) * vxD : cytoDexCenterZ * vxD;
        
        maxEndRadius = Math.sqrt(dx * dx + dy * dy + dz * dz);
        // Calculate how many samples will be taken
        final int size = (int) ((maxEndRadius - cytoDexRad) / shollStep) + 1;
        ImagePlus imgDiameter = imgSkel.duplicate();
        // save projection Mask of diameters map
        ImagePlus imgDiameterProj = doZProjection(imgDiameter);
        imgDiameter.changes = false;
        imgDiameter.close();
        imgDiameter.flush();
        // Create arrays for radii (in physical units)
        final float[] radii = new float[size];
        for (int i = 0; i < size; i++) {
                radii[i] = (int)(cytoDexRad + cytoDexEnlarge) + i * shollStep;
        }        
        // Create an array to hold results
        double[] diameter = new double[radii.length];
        // create array list for intersections
        ArrayList<Integer> intersections = new ArrayList();
        ImageStack skelStack = imgSkel.getStack();
        Point3D point = new Point3D();
        // for each sphere find voxels
        // if voxel is on the skeleton read value on image map (diameter)
        IJ.showStatus("Computing diameters ...");
        for (int i = 0; i < radii.length; i++) {
            IJ.showStatus("Computing diameters and intersections for sphere "+i+"/"+radii.length+", rad="+radii[i]);
            ArrayList<Voxel3D> voxels = imgMap.getNeighborhoodLayerList(cytoDexCenterX, cytoDexCenterY, cytoDexCenterZ, radii[i], radii[i]+1);
            int index = 0;
            for (int n = 0; n < voxels.size(); n++) {
                double value = skelStack.getVoxel((int)voxels.get(n).x, (int)voxels.get(n).y, (int)voxels.get(n).z);
                if (value == 255) {
                    index++;
                    point.x = voxels.get(n).x;
                    point.y = voxels.get(n).y;
                    point.z = voxels.get(n).z;
                    float diam = imgMap.getPixel(point)*2;
                    diameter[i] += diam;   
                }
            }
            intersections.add(i, index);
            diameter[i] /= index;    
            // put pixel on skel image inside radius n and n+1
            if (i < radii.length -1) 
                drawInsideCircle(radii[i], radii[i+1],diameter[i],imgDiameterProj);
        }
        // Write shollAnalysis results for diameters
        // write headers for global results
        FileWriter fwAnalyze = new FileWriter(imgOutDir + imgTitle.substring(0,imgTitle.indexOf(".tif")) +"_Analyze_Diameter.xls",false);
        BufferedWriter outputAnalyze = new BufferedWriter(fwAnalyze);
        outputAnalyze.write("Distance\tMean diameter ("+cal.getUnit()+")\t# Intersections\n");
        outputAnalyze.flush();
        for (int i = 0; i < radii.length; i++) {
            outputAnalyze.write(radii[i]+"\t"+diameter[i]+"\t"+intersections.get(i)+"\n");
            outputAnalyze.flush();
        }
        fwAnalyze.close();
        // Apply LUT
        final double[] range = Tools.getMinMax(diameter);
        final boolean logMask = range[1] < 0;
        final int fcolor = -1;
        final int bcolor = 0;
        imgDiameterProj.getProcessor().setColorModel(matlabJetColorMap(bcolor, fcolor));
	imgDiameterProj.getProcessor().setMinAndMax(logMask ? range[0] : 0, range[1]);
        imgDiameterProj.setCalibration(cal);
        IJ.run(imgDiameterProj,  "Calibration Bar...", "location=[Upper Right] fill=White label=Black number=5 decimal=0 font=12 zoom=2 overlay");
        IJ.run(imgDiameterProj, "Scale Bar...", "width=100 height=10 font=42 color=Yellow background=None location=[Lower Right]");
        FileSaver imgSave = new FileSaver(imgDiameterProj);
        imgSave.saveAsTiff(imgOutDir+imgTitle.substring(0,imgTitle.indexOf(".tif"))+"_ShollDiameterSkel.tif");
        imgDiameterProj.changes = false;
        imgDiameterProj.close();
        imgDiameterProj.flush();       
    }
       
}
