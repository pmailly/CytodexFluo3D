package CytoDex;

//  This plugin Analyze nucleus spheroids repartition ...

import static CytoDex.Cytodex_Fluo3D_Lite.analyzeSkeleton;
import static Tools.Cytodex_tools.*;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.WaitForUserDialog;
import ij.io.FileSaver;
import ij.plugin.Duplicator;
import ij.plugin.PlugIn;
import ij.plugin.frame.ThresholdAdjuster;
import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom.Object3D;
import mcib3d.geom.ObjectCreator3D;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Point3D;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;
import mcib3d.utils.ArrayUtil;
import org.apache.commons.io.FilenameUtils;


public class Cytodex_Fluo3D_Nucleus implements PlugIn {
    
    private boolean canceled = false;
    private String imageDir = "";
    private BufferedWriter outPutAnalyze;
    private final double nucMinVol = 150;
    private final double nucMaxVol = Double.MAX_VALUE;
    private String thresholdMethod = "Default";
    private int shollStart = 100;
    private int shollStep = 50;
    private final double minBranchVol = 100;
    
   
    
    /**
     * Ask for image folder, Sholl start diameter and step
     */
    private boolean dialog() {
       boolean canceled = false;
       GenericDialogPlus gd = new GenericDialogPlus("Sholl parameters");
       gd.addDirectoryField("Choose Directory Containing CZI Files", imageDir); 
       gd.addMessage("Sholl parameters (µm)");
       gd.addNumericField("Start diameter", shollStart, 2);
       gd.addNumericField("Step diameter", shollStep, 2);
       gd.showDialog();
       if (gd.wasCanceled())
            canceled = true;
       imageDir = gd.getNextString();
       shollStart = (int)gd.getNextNumber();
       shollStep = (int)gd.getNextNumber();
       return(canceled);
    }
    
    
    /**
     * return objects population in an binary image
     * @param img
     * @return pop objects population
     */

    private static  Objects3DPopulation getPopFromImage(ImagePlus img) {
        // label binary images first
        ImageLabeller labeller = new ImageLabeller();
        ImageInt labels = labeller.getLabels(ImageHandler.wrap(img));
        Objects3DPopulation pop = new Objects3DPopulation(labels);
        return pop;
    }
    
    /**
     * Threshold dialog
     */
    private void dialogThreshod(ImagePlus img) {
        img.getProcessor().resetThreshold();
        img.setSlice(img.getNSlices()/2);
        img.updateAndDraw();
        img.show();
        IJ.setAutoThreshold(img, "Default dark");
        IJ.run("Threshold...");
        if (!WindowManager.getWindow("Threshold").isVisible());
            WindowManager.getWindow("Threshold").setVisible(true);
        new WaitForUserDialog("Choose threshold method and press OK").show();
        thresholdMethod = new ThresholdAdjuster().getMethod();
        WindowManager.getWindow("Threshold").setVisible(false);
        System.out.println("Threshold method :"+thresholdMethod);
        img.hide();
    }
    
        
        // Threshold images and fill holes
    private static void threshold(ImagePlus img, String thMed, boolean fill, boolean calcul) {
        //  Threshold and binarize
       String cal = "";
       img.setZ(img.getNSlices()/2);
       img.updateAndDraw();
       IJ.setAutoThreshold(img, thMed+" dark");
       Prefs.blackBackground = false;
       if (calcul)
           cal = " calculate";
        IJ.run(img, "Convert to Mask", "method="+thMed+cal+" background=Dark");
        if (fill) {
            IJ.run(img,"Fill Holes", "stack");
        }
    }
    
    
    /**
     * Find nucleus inside dist2-dist1
     * @param nucPop
     * @param center
     * @param rad
     * @param radOffset
     * @return nucleus included in a volume a distance rad + offset
     */
    private Objects3DPopulation nucleusInsideDistance(Objects3DPopulation nucPop, Point3D center, double dist1, double dist2) {
        Objects3DPopulation popInside = new Objects3DPopulation();
        for (int o = 0; o < nucPop.getNbObjects(); o++) {
            Object3D obj = nucPop.getObject(o);
            double dist = obj.distPixelBorderUnit(center.getX(), center.getY(), center.getZ());
            if (dist < dist2 && dist >= dist1) {
                obj.setName(String.valueOf(dist));
                popInside.addObject(obj);
                nucPop.removeObject(obj);
                o--;
            }
        }
        return(popInside);
    }
    
    
    /**
     * Calculate how many samples will be taken for sholl analyze
     * @param img
     * @param ShollStep
     * @return 
     */
    private int shollStep(ImagePlus img, Point3D centroid, int start, int shollStep) {
    // parameters
        final int wdth = img.getWidth();
        final int hght = img.getHeight();
        final int depth = img.getNSlices();
        final double dx, dy, dz, maxEndRadius;
        final double vxWH = Math.sqrt(cal.pixelWidth * cal.pixelHeight);
	final double vxD = cal.pixelDepth;
        dx = (centroid.getRoundX() <= wdth / 2) ? (centroid.getRoundX() - wdth) * vxWH : centroid.getRoundX() * vxWH;
        dy = (centroid.getRoundY()  <= hght / 2) ? (centroid.getRoundY()  - hght) * vxWH : centroid.getRoundY()  * vxWH;
        dz = (centroid.getRoundZ()  <= depth / 2) ? (centroid.getRoundZ()  - depth) * vxD : centroid.getRoundZ()  * vxD;
        maxEndRadius = Math.sqrt(dx * dx + dy * dy + dz * dz);
        // Calculate how many samples will be taken
        final int samples = (int) ((maxEndRadius - cytoDexRad) / shollStep) + 1;
        return(samples);
    }
    
    private double[] getMeanVolume(Objects3DPopulation pop) {
        ArrayUtil vol = new ArrayUtil(pop.getNbObjects());
        double[] volume = new double[2];
        for (int o = 0; o < pop.getNbObjects(); o++) {
            Object3D obj = pop.getObject(o);
            vol.addValue(o, obj.getVolumeUnit());
        }
        volume[0] = vol.getMean();
        volume[1] = vol.getStdDev();
        return(volume);
    }
 
    
    
    @Override
    public void run(String arg) {
        String imageExt = "czi";
        try {
            if (dialog() || imageDir == null) {
                IJ.showMessage("Plugin canceled");
                return;
            }
            File inDir = new File(imageDir);
            String[] imageFile = inDir.list();

            // create OME-XML metadata store of the latest schema version
            ServiceFactory factory;
            factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            Arrays.sort(imageFile);
            int imageNum = 0; 
            String outDirResults = null;
            for (int i = 0; i < imageFile.length; i++) {
                String fileExt = FilenameUtils.getExtension(imageFile[i]);
                if (fileExt.equals(imageExt)) {
                    int series = 0;
                    String imageName = inDir+ File.separator+imageFile[i];
                    String rootName = FilenameUtils.getBaseName(imageFile[i]);
                    imageNum++;
                    reader.setId(imageName);
                    // Check calibration
                    if (imageNum == 1) {
                        // create output folder
                        outDirResults = imageDir + File.separator+ "Results"+ File.separator;
                        File outDir = new File(outDirResults);
                        if (!Files.exists(Paths.get(outDirResults))) {
                            outDir.mkdir();
                        }
                        cal.pixelWidth = meta.getPixelsPhysicalSizeX(series).value().doubleValue();
                        cal.pixelHeight = cal.pixelWidth;
                        if (meta.getPixelsPhysicalSizeZ(series) != null)
                            cal.pixelDepth = meta.getPixelsPhysicalSizeZ(series).value().doubleValue();
                        else
                            cal.pixelDepth = 1;
                        cal.setUnit("microns");
                        System.out.println("x cal = " +cal.pixelWidth+", z cal=" + cal.pixelDepth);

                        /*
                         * Write headers results for results file
                         */
                        // Global file for mito results
                        FileWriter fwAnalyze = new FileWriter(outDirResults + "Analyze_skeleton_results.xls",false);
                        outPutAnalyze = new BufferedWriter(fwAnalyze);
                        outPutAnalyze.write("Image Name\tDistance form centroid\tNumber of nucleus\tMean closest distance\tStd closest distance\t"
                                + "Mean Volume\tStd volume\n");
                    }

                    reader.setSeries(0);
                    String seriesName = meta.getImageName(0);

                    ImporterOptions options = new ImporterOptions();
                    options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                    options.setId(imageName);
                    options.setSplitChannels(true);
                    options.setQuiet(true);

                     /*
                    * Open channels
                    */

                    // Actin
                    int actinCh = 1;
                    System.out.println("Opening actin channel");
                    options.setCBegin(0, actinCh);
                    options.setCEnd(0, actinCh);
                    ImagePlus imgActin = BF.openImagePlus(options)[0];
                    
                     // actine network
                    IJ.run(imgActin, "Subtract Background...", "rolling=50 stack");
                    IJ.run(imgActin, "Laplacian of Gaussian", "sigma=5 scale_normalised negate stack");
                    dialogThreshod(imgActin);
                    threshold(imgActin, thresholdMethod, false, false);

                    Objects3DPopulation actinPop = new Objects3DPopulation(getPopFromImage(imgActin).getObjectsWithinVolume(minBranchVol, Double.MAX_VALUE, false));
                    System.out.println("Actin population = "+ actinPop.getNbObjects());

                    // Draw actin population
                    ImageHandler imhActinObjects = ImageInt.wrap(imgActin).createSameDimensions();
                    for (int o = 0; o < actinPop.getNbObjects(); o++) {
                        actinPop.getObject(o).draw(imhActinObjects, o+100);
                    } 
                    imhActinObjects.set332RGBLut();
                    imhActinObjects.setTitle(rootName + "_" + seriesName + "_ActinObjects.tif");
                    imhActinObjects.save(outDirResults);
                    imhActinObjects.closeImagePlus();
                    
                    // Skeleton
                    ImageHandler imhActinPop = ImageInt.wrap(imgActin).createSameDimensions();
                    actinPop.draw(imhActinPop, 255);
                    ImagePlus imgSkel = new Duplicator().run(imhActinPop.getImagePlus());
                    imhActinPop.closeImagePlus();
                    
                    // Find actin network morphology
                    // Skeletonize
                    imgSkel.setCalibration(cal);
                    imgSkel.setTitle(imgActin.getTitle());
                    double[] skeletonParams = analyzeSkeleton(imgSkel, outDirResults);
                    imgSkel.close();
                    imgActin.close();
                    
                    // DAPI
                    int dapiCh = 0;                        
                    options.setCBegin(0, dapiCh);
                    options.setCEnd(0, dapiCh);
                    System.out.println("-- Series : "+ seriesName);
                    System.out.println("Opening Nucleus channel");
                    ImagePlus imgNuc= BF.openImagePlus(options)[0]; 

                    // Find nucleus
                    IJ.run(imgNuc, "Difference of Gaussians", "  sigma1=6 sigma2=4 stack");
                    threshold(imgNuc, "Moments", false, false);
                    IJ.run(imgNuc, "Watershed", "stack");
                    Objects3DPopulation nucPop = new Objects3DPopulation(getPopFromImage(imgNuc).getObjectsWithinVolume(nucMinVol, nucMaxVol, true));
                    int nucNumber = nucPop.getNbObjects();
                    System.out.println("Nucleus population = "+ nucNumber);
                    // get all nucleus as one object
                    Object3D nucObj = new Objects3DPopulation(imgNuc).getObject(0);
                    
                    // Find centroid from nucleus population
                    Point3D centroid = nucObj.getCenterAsPoint();
                    System.out.println(centroid.x+","+centroid.y+","+centroid.z);
                    
                    
                    // Find nucleus population at X distance from centroid
                    // Compute how many steps we need (nucleus populations)
                    
                    ArrayList<Objects3DPopulation> nucPopList = new ArrayList();
                    int shollSteps = shollStep(imgNuc, centroid, shollStart, shollStep);
                    int shollStop = shollStart + shollStep*shollSteps;
                    int index = 0;
                    for (int d = 0; d < shollStop; d+=shollStep) {
                        double dist = d+shollStep;
                        IJ.showStatus("Computing nucleus population between "+d+" and "+dist);
                        nucPopList.add(nucleusInsideDistance(nucPop, centroid, d, dist));
                        System.out.println("nucleus at "+d+ " "+nucPopList.get(index).getNbObjects());
                        index++;
                    }
                    
                    // Save nucleus population
                    System.out.println("Saving Object population ...");
                    IJ.run(imgNuc,"32-bit","");
                    ImageHandler imhNucObjects = ImageHandler.wrap(imgNuc).createSameDimensions();
                    for (Objects3DPopulation nucleus : nucPopList) {
                        if (nucleus.getNbObjects() > 0)
                            for (int o = 0; o < nucleus.getNbObjects(); o++) {
                                Object3D obj = nucleus.getObject(o);
                                float value = Float.valueOf(obj.getName());
                                obj.draw(imhNucObjects, value);
                            }
                    }
                    new ObjectCreator3D(imhNucObjects).createSphere(centroid.getX(), centroid.getY(), centroid.getZ(), 25, 255, false);
                    IJ.run(imhNucObjects.getImagePlus(), "Cyan Hot","");
                    imhNucObjects.setTitle(rootName + "_" + seriesName + "_NucleusObjects.tif");
                    imhNucObjects.save(outDirResults);
                    imhNucObjects.closeImagePlus();
                    imgNuc.close();
                    
                    // compute mean closest distances                    
                    // Write data
                    
                    int popIndex = 0;
                    for (int d = 0; d < shollStop; d+=shollStep) {
                        double dist = d+shollStep;
                        IJ.showStatus("Computing parameters diameter "+ dist);
                        Objects3DPopulation pop = nucPopList.get(popIndex);
                        if (pop.getNbObjects() > 1) {
                            ArrayUtil allDist = pop.distancesAllClosestCenter();
                            double meanDist = allDist.getMean();
                            double stdDist = allDist.getStdDev();
                            double[] vol = getMeanVolume(pop);
                            if (popIndex == 0)
                                outPutAnalyze.write(rootName+"\t"+dist+"\t"+pop.getNbObjects()+"\t"+meanDist+"\t"+stdDist+"\t"+vol[0]+"\t"+vol[1]+"\n");
                            else
                                outPutAnalyze.write("\t"+dist+"\t"+pop.getNbObjects()+"\t"+meanDist+"\t"+stdDist+"\t"+vol[0]+"\t"+vol[1]+"\n");
                        }
                        else
                            if (popIndex == 0)
                                outPutAnalyze.write(rootName+"\t"+dist+"\t0\t0\t0\t0\t0\n");
                            else
                                outPutAnalyze.write("\t"+dist+"\t0\t0\t0\t0\t0\n");
                        outPutAnalyze.flush();
                        popIndex++;
                    }
                    

                }                    
            } 
            outPutAnalyze.close();
            IJ.showStatus("End of process");
        } catch (IOException | DependencyException | ServiceException | FormatException ex) {
            Logger.getLogger(Cytodex_Fluo3D_Nucleus.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
