package CytoDex;

//  This plugin Analyze nucleus spheroids repartition ...

import Tools.Cytodex_Actin;
import static Tools.Cytodex_tools.*;
import Tools.Cytodex_Nucleus;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.WaitForUserDialog;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.Duplicator;
import ij.plugin.PlugIn;
import ij.plugin.RGBStackMerge;
import ij.process.AutoThresholder;
import ij.process.ImageProcessor;
import java.awt.Font;
import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Vector;
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
import mcib3d.geom.Vector3D;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.math3.exception.MathRuntimeException;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;


public class Cytodex_Fluo3D_Nucleus implements PlugIn {
    
    private boolean canceled = false;
    private String imageDir = "";
    private BufferedWriter outPutAnalyze, outPutGlobalNucleus, outPutNucleus, outPutActin;
    private double nucMinVol = 150;
    private final double nucMaxVol = Double.MAX_VALUE;
    private int shollStart = 100;
    private int shollStep = 50;
    private double actinMinVol = 1000;
    private double actinMinLength = 25;
    private int iterPruning = 5;
    public Calibration cal = new Calibration();
    
    
   
    
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
       gd.addNumericField("Actin minimum volume : ", actinMinVol);
       gd.addNumericField("Actin minimum branch length : ", actinMinLength);
       gd.addNumericField("Iteration number to prune : ", iterPruning);
       gd.addNumericField("Nucleus minimum volume : ", nucMinVol);
       gd.showDialog();
       if (gd.wasCanceled())
            canceled = true;
       imageDir = gd.getNextString();
       shollStart = (int)gd.getNextNumber();
       shollStep = (int)gd.getNextNumber();
       actinMinVol = gd.getNextNumber();
       actinMinLength = gd.getNextNumber();
       iterPruning = (int)gd.getNextNumber();
       nucMinVol = gd.getNextNumber();
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
     * Find nucleus inside dist2-dist1
     * @param nucPop
     * @param ActinPop
     * @param center
     * @param rad
     * @param radOffset
     * @return nucleus included inside dist2-dist1
     */
    private ArrayList<Cytodex_Nucleus> nucleusInsideDistance(ArrayList<Cytodex_Nucleus> nucList, double dist1, double dist2) {
        ArrayList<Cytodex_Nucleus> nucInsideList = new ArrayList<>();
        for (int i = 0; i < nucList.size(); i++) {
            double distObj = nucList.get(i).getCenterDistToCenter();
            //System.out.println("Nucleus ="+i+" d1 ="+dist1+" d2 ="+dist2+" dist ="+distObj);
            if (distObj >= dist1 && distObj < dist2) {
                nucInsideList.add(nucList.get(i));
            }
        }
        return(nucInsideList);
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
        final int samples = (int) ((maxEndRadius - start) / shollStep) + 1;
        return(samples);
    }
    
    
    private Object3D getAllObjectsInOne(Objects3DPopulation pop, ImageHandler imh) {
        ImageHandler imhPop = imh.duplicate();
        pop.draw(imhPop, 255);
        Objects3DPopulation objs = new Objects3DPopulation(imhPop, 1);
        imhPop.closeImagePlus();
        return(objs.getObject(0));
    }
    
    
    /**
     * 
     * 
    */
     private ArrayList<Cytodex_Actin> actinParameters (Objects3DPopulation actinPop, Point3D centroid, Vector3D allNucVector, ImagePlus img, String out, String name, String series) {
        ImageHandler imhActinObjects = ImageInt.wrap(img).createSameDimensions();
        ArrayList<Cytodex_Actin> actinParams = new ArrayList<>();
        int index = 1;
        for (int n = 0; n < actinPop.getNbObjects(); n++) {
            Object3D actinObj = actinPop.getObject(n);
            double distCenter = actinObj.distPixelCenter(centroid);
            double distBorder = actinObj.distPixelBorderUnit(centroid.x, centroid.y, centroid.z);
            double vol = actinObj.getVolumeUnit();
            Vector3D actinVec = actinObj.getMainAxis();
            double allAngle = angleBetweenVectors(actinVec, allNucVector);
            Vector3D actinVecCentroid = actinObj.vectorPixelUnitBorder(centroid.getX(), centroid.getY(), centroid.getZ());
            double centroidAngle = angleBetweenVectors(actinVec, actinVecCentroid);
            actinObj.draw(imhActinObjects, index);
            Cytodex_Actin actin = new Cytodex_Actin(index, vol, distCenter, distBorder, 0, centroidAngle, allAngle);
            actinParams.add(actin);
            index++;
        } 
        imhActinObjects.set332RGBLut();
        imhActinObjects.setTitle(name + "_" + series + "_ActinObjects.tif");
        imhActinObjects.save(out);
        imhActinObjects.closeImagePlus();
        return(actinParams);
     }
    
    
    /**
     * Find nucleus number in actin object
     * @param nucPop
     * @param actinPop
     * @return nucleusNumber
     */
    
    private void findNucleusInActin(Objects3DPopulation actinPop, Objects3DPopulation nucPop, ArrayList<Cytodex_Actin> actin) {
        for (int n = 0; n < actinPop.getNbObjects(); n++) {
            Object3D actinObj = actinPop.getObject(n);
            int nuc = 0;
            for (int i = 0; i < nucPop.getNbObjects(); i++) {
                Object3D nucObj = nucPop.getObject(i);
                if (actinObj.hasOneVoxelColoc(nucObj))
                    nuc++;
            }
            actin.get(n).setNucNumber(nuc);
        }
    }
    
    
    
    private int actinIndex(Object3D nucObj, Objects3DPopulation actinPop) {
        int actinNumber = -1;
        int totalActin = actinPop.getNbObjects();
        for (int n = 0; n < totalActin; n++) {
            Object3D actinObj = actinPop.getObject(n);
            if (nucObj.hasOneVoxelColoc(actinObj)) {
                actinNumber = (n+1);
                n = totalActin;
            }
        }
        return(actinNumber);
    }
    
    
    private double angleBetweenVectors(Vector3D v1, Vector3D v2) {
        double angle = (Math.acos((v1.x * v2.x + v1.y * v2.y)/((Math.sqrt(v1.x*v1.x + v1.y*v1.y) * Math.sqrt(v2.x*v2.x + v2.y*v2.y))))/Math.PI*180);
        return(angle);
    }
        
    /**
     * Find nucleus parameters
     */
    private ArrayList<Cytodex_Nucleus> nucleusParameters(Objects3DPopulation nucPop, Objects3DPopulation actinPop, Point3D centroid, Vector3D allNucVector, ImagePlus actinMap) {
        ArrayList<Cytodex_Nucleus> nucleusList = new ArrayList();
        ImageHandler imhActin = ImageHandler.wrap(actinMap);
        int index = 1;
        for (int n = 0; n < nucPop.getNbObjects(); n++) {
            Object3D nucObj = nucPop.getObject(n);
            double vol = nucObj.getVolumeUnit();
            double sph = nucObj.getSphericity(true);
            double comp = nucObj.getCompactness(true);
            double distCenter = nucObj.distPixelBorderUnit(centroid.getX(), centroid.getY(), centroid.getZ());
            Vector3D nucVec = nucObj.getMainAxis();
            double allAngle = angleBetweenVectors(nucVec, allNucVector);
            Vector3D nucVecCentroid = nucObj.vectorPixelUnitBorder(centroid.getX(), centroid.getY(), centroid.getZ());
            double centroidAngle = angleBetweenVectors(nucVec, nucVecCentroid);
            int actinNumber = actinIndex(nucObj, actinPop);
            double diam = 0;
            if (actinNumber != -1)
                diam = nucObj.getPixCenterValue(imhActin);
            Object3D closestObj = nucPop.closestCenter(nucObj, true);
            double closestDist = nucObj.distCenterUnit(closestObj);
            Cytodex_Nucleus nucleus = new Cytodex_Nucleus(index, vol, distCenter, sph, comp, centroidAngle, allAngle, actinNumber, closestDist, diam);
            nucleusList.add(nucleus);
            index++;
        }
        return(nucleusList);
    }
    
    
    /**
     * Write parameters for each nucleus
     * @param arg 
     */
    private void writeNucParameters(ArrayList<Cytodex_Nucleus> nucList, ImagePlus img, Point3D centroid, String name) throws IOException {
        int shollSteps = shollStep(img, centroid, shollStart, shollStep);
        int shollStop = shollStart + shollStep*shollSteps;
        for (int d = shollStart; d < shollStop; d+=shollStep) {
            int shollDiameter = d + shollStep;
            ArrayList<Cytodex_Nucleus> inSphere = nucleusInsideDistance(nucList, d, shollDiameter);
            if (inSphere.size() > 1) {
                for (int n = 0; n < inSphere.size(); n++) {
                    int index = inSphere.get(n).getIndex();
                    double distToCenter = inSphere.get(n).getCenterDistToCenter();
                    double closestDist = inSphere.get(n).getNucClosestDist();
                    double volume = inSphere.get(n).getVolume();
                    double sph = inSphere.get(n).getSphericity();
                    double comp = inSphere.get(n).getCompactness();
                    double allAngle = inSphere.get(n).getAllAngle();
                    double centroidAngle = inSphere.get(n).getCentroidAngle();
                    int actinObj = inSphere.get(n).getActinObjectIndex();
                    double ActinDiameter = 2 *  inSphere.get(n).getActinDiameter();
                    if (n == 0)
                        outPutNucleus.write(name+"\t"+d+"-"+shollDiameter+"\t"+index+"\t"+distToCenter+"\t"+volume+"\t"+closestDist+"\t"+actinObj+"\t"
                                +ActinDiameter+"\t"+sph+"\t"+comp+"\t"+allAngle+"\t"+centroidAngle+"\n");
                    else
                        outPutNucleus.write("\t\t"+index+"\t"+distToCenter+"\t"+volume+"\t"+closestDist+"\t"+actinObj+"\t"
                                +ActinDiameter+"\t"+sph+"\t"+comp+"\t"+allAngle+"\t"+centroidAngle+"\n");
                    outPutNucleus.flush();
                }
            }
        }
    }
    
    
    /**
     * 
     */
    private void writeNucleusGlobalParameters(ArrayList<Cytodex_Nucleus> nucleusList, Object3D nucObjs, ImagePlus imgNuc, Point3D centroid, String rootName) throws IOException {
        double globalSphericity = nucObjs.getSphericity(true);
        double globalCompactness = nucObjs.getCompactness(true);
        Vector3D mainAxe = nucObjs.getMainAxis();
        Vector3D center = new Vector3D(centroid);
        double globalAngle = mainAxe.angle(center);
        int shollSteps = shollStep(imgNuc, centroid, shollStart, shollStep);
        int shollStop = shollStart + shollStep*shollSteps;
        for (int d = shollStart; d < shollStop; d+=shollStep) {
            int shollDiameter = d + shollStep;
            DescriptiveStatistics allDist = new DescriptiveStatistics();
            DescriptiveStatistics allVol = new DescriptiveStatistics();
            DescriptiveStatistics allSph = new DescriptiveStatistics();
            DescriptiveStatistics allComp = new DescriptiveStatistics();
            DescriptiveStatistics allAngle = new DescriptiveStatistics();
            DescriptiveStatistics centroidAngle = new DescriptiveStatistics();
            ArrayList<Cytodex_Nucleus> inSphere = nucleusInsideDistance(nucleusList, d, shollDiameter);
            boolean findNuc = false;
            for (int n = 0; n < inSphere.size(); n++) {
                if (inSphere.get(n).getCenterDistToCenter() >= d && inSphere.get(n).getCenterDistToCenter() < shollDiameter) {
                    allDist.addValue(inSphere.get(n).getNucClosestDist());
                    allVol.addValue(inSphere.get(n).getVolume());
                    allSph.addValue(inSphere.get(n).getSphericity());
                    allComp.addValue(inSphere.get(n).getCompactness());
                    allAngle.addValue(inSphere.get(n).getAllAngle());
                    centroidAngle.addValue(inSphere.get(n).getCentroidAngle());
                    findNuc = true;
                }
            }
            if (findNuc) {
                if (d == shollStart)
                    outPutGlobalNucleus.write(rootName+"\t"+d+"-"+shollDiameter+"\t"+inSphere.size()+"\t"+allDist.getMean()+"\t"+allDist.getStandardDeviation()+"\t"
                            +allVol.getMean()+"\t"+allVol.getStandardDeviation()+"\t"+allSph.getMean()+"\t"+allSph.getStandardDeviation()+"\t"+allComp.getMean()+"\t"
                            +allComp.getStandardDeviation()+"\t"+allAngle.getMean()+"\t"+allAngle.getStandardDeviation()+"\t"+centroidAngle.getMean()+"\t"
                            +centroidAngle.getStandardDeviation()+"\t"+globalSphericity+"\t"+globalCompactness+"\t"+globalAngle+"\n");
                else
                    outPutGlobalNucleus.write("\t"+d+"-"+shollDiameter+"\t"+inSphere.size()+"\t"+allDist.getMean()+"\t"+allDist.getStandardDeviation()+"\t"
                            +allVol.getMean()+"\t"+allVol.getStandardDeviation()+"\t"+allSph.getMean()+"\t"+allSph.getStandardDeviation()+"\t"+allComp.getMean()+"\t"
                            +allComp.getStandardDeviation()+"\t"+allAngle.getMean()+"\t"+allAngle.getStandardDeviation()+"\t"+centroidAngle.getMean()+"\t"
                            +centroidAngle.getStandardDeviation()+"\t\t\t\n");
            }
            else
                outPutGlobalNucleus.write("\t"+d+"-"+shollDiameter+"\t0\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n");
            outPutGlobalNucleus.flush();
        }
    }
    
    /*
    * clear inside spheroid
    */
    public void removeSpheroid(ImagePlus img, int color, Point3D center) {
        IJ.showStatus("Removing spheroid center ...");
        Calibration cal = img.getCalibration();
        ObjectCreator3D sphereC0 = new ObjectCreator3D(img.getStack());
        double sphereRad = shollStart / cal.pixelWidth;
        //System.out.println("sphere rad = " + shollStart + shollStep);
        sphereC0.createEllipsoid(center.getRoundX(), center.getRoundY(), center.getRoundZ(), sphereRad, sphereRad, sphereRad/cal.pixelDepth, color, false);
        img.updateAndDraw();
    } 
    
     /**
     * Label object
     * @param popObj
     * @param img 
     */
    public void labelsObject (Object3D obj, int n, ImageHandler img) {
        int fontSize = Math.round(10f/(float)cal.pixelWidth);
        Font tagFont = new Font("SansSerif", Font.PLAIN, fontSize);
        int[] box = obj.getBoundingBox();
        int z = (int)obj.getCenterZ();
        int x = box[0] - 2;
        int y = box[2] - 2;
        img.getImagePlus().setSlice(z+1);
        ImageProcessor ip = img.getImagePlus().getProcessor();
        ip.setFont(tagFont);
        ip.setColor(255);
        ip.drawString(Integer.toString(n), x, y);
        img.getImagePlus().updateAndDraw();
    }
    
    /**
     * Save object Populations
     * @param nucPop
     * @param actinPop
     * @param imageName
     * @param img 
     */
    private void saveImgObjects(Objects3DPopulation nucPop, Objects3DPopulation actinPop, String imageName, ImagePlus img) {
        //create image objects population
        ImageHandler imgNucObj = ImageInt.wrap(img).createSameDimensions();
        ImageHandler imgActinObj = ImageInt.wrap(img).createSameDimensions();
        //nuc population green
        nucPop.draw(imgNucObj, 255);
        //vessel population red       
        actinPop.draw(imgActinObj, 255);
        
        // save image for objects population
        ImagePlus[] imgColors = {imgNucObj.getImagePlus(), imgActinObj.getImagePlus()};
        ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, false);
        imgObjects.setCalibration(cal);
        FileSaver ImgObjectsFile = new FileSaver(imgObjects);
        ImgObjectsFile.saveAsTiff(imageName + "_Objects.tif"); 
        imgNucObj.closeImagePlus();
        imgActinObj.closeImagePlus();
    }
    
    @Override
    public void run(String arg) {
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
                if (fileExt.equals("czi") || fileExt.equals("nd")) {
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
                        // Global file for analyze skeleton results
                        FileWriter fwAnalyze = new FileWriter(outDirResults + "Skeleton_results.xls",false);
                        outPutAnalyze = new BufferedWriter(fwAnalyze);
                        outPutAnalyze.write("Image\t#Skeletons\t#Branches\t#Junctions\t#End points\tTotal branch length\n");
                        
                        // Global file for nucleus results
                        FileWriter fwNucleus = new FileWriter(outDirResults + "Global_Nucleus_results.xls",false);
                        outPutGlobalNucleus = new BufferedWriter(fwNucleus);
                        outPutGlobalNucleus.write("Image Name\tSphere diameter\tNumber of nucleus\tMean closest distance\tStd closest distance"
                                + "\tNucleus mean Volume\tNucleus std volume\tMean sphericity\tStd sphericity\tMean Compactness\tStd compactness"
                                + "\tMean orientation to all nucleus\tStd orientation to all nucleus\tMean centroid orientation\tStd centroid orientation\tAll nucleus global Sphericity"
                                + "\tAll nucleus global compactness\tAll nucleus global Orientation\n");
                        
                        // file for analyze nucleus results
                        FileWriter fwEachNucleus = new FileWriter(outDirResults + "Nucleus_results.xls",false);
                        outPutNucleus = new BufferedWriter(fwEachNucleus);
                        outPutNucleus.write("Image Name\tSphere diameter\t#Nucleus\tDistance from centroid\tNucleus Volume\tDistance to closest nucleus"
                                + "\t#Actin\tActin diameter\tSphericity\tCompactness\tOrientation to all nucleus\tOrientation to centroid\n");
                        
                        // file for analyze actin results
                        FileWriter fwEachActin = new FileWriter(outDirResults + "Actin_results.xls",false);
                        outPutActin = new BufferedWriter(fwEachActin);
                        outPutActin.write("Image Name\t#Actin\tCenter distance to centroid\tBorder distance to centroid\tActine volume\tNucleus number\tOrientation to all nucleus"
                                + "\tOrientation to centroid nucleus\n");
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

                    // Actin find the tif file corresponding to Imaris mask
                    String actinFile = imageDir + File.separator + rootName + ".tif";
                    if (!new File(actinFile).exists()) {
                        IJ.showStatus("No actin file found !!!");
                        break;
                    }
                    else {    
                        System.out.println("Opening actin channel");
                        ImagePlus imgActin = IJ.openImage(actinFile);

                         // actine network
                        Objects3DPopulation actinPop = new Objects3DPopulation(getPopFromImage(imgActin).getObjectsWithinVolume(actinMinVol, Double.MAX_VALUE, true));
                        System.out.println("Actin population = "+ actinPop.getNbObjects());

                        
                        IJ.run(imgActin, "8-bit", "");
                        ImageHandler imhActinPop = ImageInt.wrap(imgActin).createSameDimensions();
                        imhActinPop.setTitle(rootName);
                        actinPop.draw(imhActinPop, 255);
                        
                        // get actin map
                        System.out.println("Find actin map ...");
                        ImagePlus actinMap = localThickness3D(imhActinPop.getImagePlus());
                                

                        // DAPI
                        int dapiCh = 0;                        
                        options.setCBegin(0, dapiCh);
                        options.setCEnd(0, dapiCh);
                        System.out.println("-- Series : "+ seriesName);
                        System.out.println("Opening Nucleus channel");
                        ImagePlus imgNuc= BF.openImagePlus(options)[0]; 

                        // Find nucleus
                        IJ.run(imgNuc, "Difference of Gaussians", "  sigma1=6 sigma2=4 stack");
                        threshold(imgNuc, AutoThresholder.Method.Moments, false, false);
                        IJ.run(imgNuc, "Watershed", "stack");

                        Objects3DPopulation nucPop = new Objects3DPopulation(getPopFromImage(imgNuc).getObjectsWithinVolume(nucMinVol, nucMaxVol, true));
                        int nucNumber = nucPop.getNbObjects();
                        System.out.println("Nucleus population = "+ nucNumber);
                        
                        // get all nucleus as one object
                        // Find centroid and fit elipse from nucleus population
                        ImageHandler imhNucObjects = ImageHandler.wrap(imgNuc).createSameDimensions();
                        Object3D nucObjTotal = getAllObjectsInOne(nucPop, imhNucObjects);
                        Point3D centroid = nucObjTotal.getCenterAsPoint();
                        Vector3D nucVector = nucObjTotal.getMainAxis();
                        //System.out.println(centroid.x+","+centroid.y+","+centroid.z);


                        // Nucleus parameters
                        ArrayList<Cytodex_Nucleus> nucleusList = nucleusParameters(nucPop, actinPop, centroid, nucVector, actinMap);
                        
                        // Save nucleus population
                        System.out.println("Saving Object population ...");
                        IJ.run(imgNuc,"32-bit","");
                        for (int n = 0; n < nucleusList.size(); n++) {
                            Object3D nucObj = nucPop.getObject(n);
                            float value = (float)nucleusList.get(n).getCenterDistToCenter();
                            nucObj.draw(imhNucObjects, value);
                            labelsObject(nucObj, nucleusList.get(n).getIndex(), imhNucObjects);
                        }
                        new ObjectCreator3D(imhNucObjects).createEllipsoidUnit(centroid.getX(), centroid.getY(), centroid.getZ(), 20, 20, 20/cal.pixelDepth, 255, false);
                        IJ.run(imhNucObjects.getImagePlus(), "Cyan Hot","");
                        imhNucObjects.setTitle(rootName + "_" + seriesName + "_NucleusObjects.tif");
                        imhNucObjects.save(outDirResults);
                        imhNucObjects.closeImagePlus();

                        // Write parameters for each nucleus
                        writeNucParameters(nucleusList, imgNuc, centroid, rootName);
                        imgNuc.close();
                        
//                        // Skeleton
//                        ImagePlus imgSkel = new Duplicator().run(imhActinPop.getImagePlus());
//                        imhActinPop.closeImagePlus();
//                        removeSpheroid(imgSkel, 0, centroid);
//                        // Find actin network morphology
//                        // Skeletonize
//                        imgSkel.setCalibration(cal);
//                        imgSkel.setTitle(imgActin.getTitle());
//                        IJ.run(imgSkel, "Skeletonize (2D/3D)", "");
//                        
//                        analyzeSkel(imgSkel, outPutAnalyze, outDirResults, actinMinLength, iterPruning);
//                        imgSkel.close();
//                        imgActin.close();
                        
                        // Get actin parameters
                        ArrayList<Cytodex_Actin> actinParams = actinParameters (actinPop, centroid, nucVector, imgActin, outDirResults, rootName, seriesName);
                        // find nucleus in actin
                        
                        // write actin results
                        //get nucleus number in actin object
                        findNucleusInActin(actinPop, nucPop, actinParams);       
                        for (int p = 0; p < actinParams.size(); p++) {
                            if (p == 0)
                                outPutActin.write(rootName+"\t"+actinParams.get(p).getIndex()+"\t"+actinParams.get(p).getCenterDistToCenter()+"\t"
                                        +actinParams.get(p).getBorderDistToCenter()+"\t"+actinParams.get(p).getVolume()+"\t"+actinParams.get(p).getNucNumber()+"\t"
                                        +actinParams.get(p).getAllAngle()+"\t"+actinParams.get(p).getCentroidAngle()+"\n");
                            else
                                outPutActin.write("\t"+actinParams.get(p).getIndex()+"\t"+actinParams.get(p).getCenterDistToCenter()+"\t"
                                        +actinParams.get(p).getBorderDistToCenter()+"\t"+actinParams.get(p).getVolume()+"\t"+actinParams.get(p).getNucNumber()+"\t"
                                        +actinParams.get(p).getAllAngle()+"\t"+actinParams.get(p).getCentroidAngle()+"\n");
                            outPutActin.flush();
                        }                        
                        // compute global nucleus parameters                    
                        // Write data
                        writeNucleusGlobalParameters(nucleusList, nucObjTotal, imgNuc, centroid, rootName);
                        
                        // save image with nucleus (red) and actin (green) objects
                        saveImgObjects(nucPop, actinPop, outDirResults+rootName, imgNuc);
                        
                        imgNuc.close();
                    }

                }                    
            } 
            outPutAnalyze.close();
            outPutActin.close();
            outPutGlobalNucleus.close();
            outPutNucleus.close();
            IJ.showStatus("End of process");
        } catch (IOException | DependencyException | ServiceException | FormatException ex) {
            Logger.getLogger(Cytodex_Fluo3D_Nucleus.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
