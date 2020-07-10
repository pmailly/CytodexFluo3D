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
import ij.measure.Calibration;
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
import mcib3d.geom.Object3DVoxels;
import mcib3d.geom.ObjectCreator3D;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Point3D;
import mcib3d.geom.Voxel3D;
import mcib3d.image3d.ImageFloat;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;
import mcib3d.image3d.distanceMap3d.EDT;
import mcib3d.utils.ArrayUtil;
import org.apache.commons.io.FilenameUtils;
import sc.fiji.localThickness.Clean_Up_Local_Thickness;
import sc.fiji.localThickness.EDT_S1D;
import sc.fiji.localThickness.Local_Thickness_Parallel;


public class Cytodex_Fluo3D_Nucleus implements PlugIn {
    
    private boolean canceled = false;
    private String imageDir = "";
    private BufferedWriter outPutAnalyze, outPutGlobalNucleus, outPutNucleus, outPutActin;
    private double nucMinVol = 150;
    private final double nucMaxVol = Double.MAX_VALUE;
    private String thresholdMethod = "Default";
    private int shollStart = 100;
    private int shollStep = 50;
    private double actinMinVol = 1000;
    private double actinMinLength = 25;
    private final double actinColoc = 50;
    
   
    
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
       gd.addNumericField("Nucleus minimum volume : ", nucMinVol);
       gd.showDialog();
       if (gd.wasCanceled())
            canceled = true;
       imageDir = gd.getNextString();
       shollStart = (int)gd.getNextNumber();
       shollStep = (int)gd.getNextNumber();
       actinMinVol = gd.getNextNumber();
       actinMinLength = gd.getNextNumber();
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
    
    private double[] findVolOfColoc(Object3D obj, Objects3DPopulation pop) {
        double vol = 0;
        double index = -1;
        for (int o = 0; o < pop.getNbObjects(); o++) {
            Object3D actinObj = pop.getObject(o);
            int coloc = obj.getColoc(actinObj);
            if (coloc >= actinColoc) {
                vol = actinObj.getVolumeUnit();
                index = o;
                break;
            }
        }
        double[] p = {index, vol};
        return(p);
    }
    
    
    
    
    /**
     * Find nucleus inside dist2-dist1
     * @param nucPop
     * @param ActinPop
     * @param center
     * @param rad
     * @param radOffset
     * @return nucleus included in a volume a distance rad + offset and add volume actin object
     */
    private Objects3DPopulation nucleusInsideDistance(Objects3DPopulation nucPop, Objects3DPopulation actinPop, Point3D center, double dist1, double dist2) {
        Objects3DPopulation popInside = new Objects3DPopulation();
        for (int o = 0; o < nucPop.getNbObjects(); o++) {
            Object3D obj = nucPop.getObject(o);
            double dist = obj.distPixelCenter(center);
            if (dist < dist2 && dist >= dist1) {
                double[] actinparams = findVolOfColoc(obj, actinPop);
                obj.setName(dist + "," + actinparams[0]);
                obj.setValue((int) actinparams[1]);
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
     private ArrayList<double[]> actinParameters (Objects3DPopulation actinPop, ImagePlus img, String out, String name, String series) {
        ImageHandler imhActinObjects = ImageInt.wrap(img).createSameDimensions();
        Object3D actinTotalObj = getAllObjectsInOne(actinPop, imhActinObjects);
        Point3D ptCenter = actinTotalObj.getCenterAsPoint();
        ArrayList<double[]> actinParam = new ArrayList<>();
        for (int o = 0; o < actinPop.getNbObjects(); o++) {
            Object3D actinObj = actinPop.getObject(o);
            double dist = actinObj.distPixelCenter(ptCenter);
            double vol = actinObj.getVolumeUnit();
            actinObj.draw(imhActinObjects, (o+1));
            double[] params = {dist, vol};
            actinParam.add(params);
        } 
        imhActinObjects.set332RGBLut();
        imhActinObjects.setTitle(name + "_" + series + "_ActinObjects.tif");
        imhActinObjects.save(out);
        imhActinObjects.closeImagePlus();
        return(actinParam);
     }
    
    
    /**
     * Find either volume object or volume of coloc object
     * @param pop
     * @return 
     */
    
    private double[] getMeanVolume(Objects3DPopulation pop, boolean coloc) {
        ArrayUtil vol = new ArrayUtil(pop.getNbObjects());
        double[] volume = new double[2];
        for (int o = 0; o < pop.getNbObjects(); o++) {
            Object3D obj = pop.getObject(o);
            if (coloc)
                vol.addValue(o, obj.getValue());
            else
                vol.addValue(o, obj.getVolumeUnit());
        }
        volume[0] = vol.getMean();
        volume[1] = vol.getStdDev();
        return(volume);
    }
 
    /**
     * compute local thickness
     * @param imgAstro
     * @return astroMap
    **/
    public ImagePlus localThickness3D (ImagePlus imgMask) {
        Calibration cal = imgMask.getCalibration();
        ImagePlus img = imgMask.duplicate();
        img.setCalibration(cal);
        threshold(img, "Li", false, false);
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
        ImagePlus astroMap = cleanUp.getResultImage();
        // calibrate intensity to µm
        cal.setValueUnit("µm");
        astroMap.setCalibration(cal);
        IJ.run(astroMap, "Multiply...", "value="+cal.pixelWidth+" stack");
        img.close();
        imgEDT.close();
        imgLocThk.close();
        return(astroMap);
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
                        // Global file for analyze skeleton results
                        FileWriter fwAnalyze = new FileWriter(outDirResults + "Skeleton_results.xls",false);
                        outPutAnalyze = new BufferedWriter(fwAnalyze);
                        outPutAnalyze.write("Image\t#Skeletons\t#Branches\t#Junctions\t#End points\tBranch length\n");
                        
                        FileWriter fwNucleus = new FileWriter(outDirResults + "Global_Nucleus_results.xls",false);
                        outPutGlobalNucleus = new BufferedWriter(fwNucleus);
                        outPutGlobalNucleus.write("Image Name\tSphere diameter\tNumber of nucleus\tMean closest distance\tStd closest distance"
                                + "\tNucleus mean Volume\tNucleus std volume\tActin mean volume\t Actin std volume\n");
                        
                        // file for analyze nucleus results
                        FileWriter fwEachNucleus = new FileWriter(outDirResults + "Nucleus_results.xls",false);
                        outPutNucleus = new BufferedWriter(fwEachNucleus);
                        outPutNucleus.write("Image Name\tSphere diameter\t#Nucleus\tDistance from centroid\tNucleus Volume\t#Actin\tActin coloc volume\tActin diameter\n");
                        
                        // file for analyze actin results
                        FileWriter fwEachActin = new FileWriter(outDirResults + "Actin_results.xls",false);
                        outPutActin = new BufferedWriter(fwEachActin);
                        outPutActin.write("Image Name\t#Actin\tDistance from centroid\tActin Volume\n");
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

                        Objects3DPopulation actinPop = new Objects3DPopulation(getPopFromImage(imgActin).getObjectsWithinVolume(actinMinVol, Double.MAX_VALUE, false));
                        System.out.println("Actin population = "+ actinPop.getNbObjects());

                        // Get actin parameters
                        ArrayList<double[]> actinParams = actinParameters (actinPop, imgActin, outDirResults, rootName, seriesName);
                        for (int p = 0; p < actinParams.size(); p++) {
                            if (p == 0)
                                outPutActin.write(rootName+"\t"+p+"\t"+actinParams.get(p)[0]+"\t"+actinParams.get(p)[1]+"\n");
                            else
                                outPutActin.write("\t\t"+p+"\t"+actinParams.get(p)[0]+"\t"+actinParams.get(p)[1]+"\n");
                            outPutActin.flush();
                        }
                        outPutActin.close();
                        IJ.run(imgActin, "8-bit", "");
                        ImageHandler imhActinPop = ImageInt.wrap(imgActin).createSameDimensions();
                        imhActinPop.setTitle(rootName);
                        actinPop.draw(imhActinPop, 255);
                        
                        // get actin map
                        ImagePlus actinMap = localThickness3D(imhActinPop.getImagePlus());
                                
                        // Skeleton
                        ImagePlus imgSkel = new Duplicator().run(imhActinPop.getImagePlus());
                        imhActinPop.closeImagePlus();
                        
                        // Find actin network morphology
                        // Skeletonize
                        imgSkel.setCalibration(cal);
                        imgSkel.setTitle(imgActin.getTitle());
                        IJ.run(imgSkel, "Skeletonize (2D/3D)", "");
                        analyzeSkel(imgSkel, outPutAnalyze, actinMinLength);
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
                        Object3D nucObj = nucPop.getMask();

                        // Find centroid from nucleus population
                        ImageHandler imhNucObjects = ImageHandler.wrap(imgNuc).createSameDimensions();
                        Object3D nucObjTotal = getAllObjectsInOne(nucPop, imhNucObjects);
                        Point3D centroid = nucObjTotal.getCenterAsPoint();
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
                            nucPopList.add(nucleusInsideDistance(nucPop, actinPop, centroid, d, dist));
                            //System.out.println("nucleus at "+d+ " "+nucPopList.get(index).getNbObjects());
                            index++;
                        }

                        // Save nucleus population
                        System.out.println("Saving Object population ...");
                        IJ.run(imgNuc,"32-bit","");
                        for (Objects3DPopulation nucleus : nucPopList) {
                            if (nucleus.getNbObjects() > 0)
                                for (int o = 0; o < nucleus.getNbObjects(); o++) {
                                    Object3D obj = nucleus.getObject(o);
                                    String[] objParam = obj.getName().split(",");
                                    float value = Float.valueOf(objParam[0]);
                                    obj.draw(imhNucObjects, value);
                                }
                        }
                        new ObjectCreator3D(imhNucObjects).createSphere(centroid.getX(), centroid.getY(), centroid.getZ(), 25, 255, false);
                        IJ.run(imhNucObjects.getImagePlus(), "Cyan Hot","");
                        imhNucObjects.setTitle(rootName + "_" + seriesName + "_NucleusObjects.tif");
                        imhNucObjects.save(outDirResults);
                        imhNucObjects.closeImagePlus();
                        imgNuc.close();

                        // Write parameters for each nucleus
                        int popIndex = 0;
                        ImageHandler imgMapActin = ImageHandler.wrap(actinMap);
                        for (int d = 0; d < shollStop; d+=shollStep) {
                            double diameter = d + shollStep;
                            Objects3DPopulation pop = nucPopList.get(popIndex);
                            int nuc = pop.getNbObjects();
                            if (nuc > 1) {
                                for (int o = 0; o < nuc; o++) {
                                    Object3D obj = pop.getObject(o);
                                    String[] objParam = obj.getName().split(",");
                                    if (o == 0)
                                        outPutNucleus.write(rootName+"\t"+diameter+"\t"+o+"\t"+objParam[0]+"\t"+obj.getVolumeUnit()+"\t"+objParam[1]+"\t"+obj.getValue()+"\t"+obj.getPixCenterValue(imgMapActin)+"\n");
                                    else
                                        outPutNucleus.write("\t\t"+o+"\t"+objParam[0]+"\t"+obj.getVolumeUnit()+"\t"+objParam[1]+"\t"+obj.getValue()+"\t"+obj.getPixCenterValue(imgMapActin)+"\n");
                                    outPutNucleus.flush();
                                }
                            }
                            popIndex++;
                        }

                        // compute mean closest distances                    
                        // Write data

                        popIndex = 0;
                        for (int d = 0; d < shollStop; d+=shollStep) {
                            double dist = d+shollStep;
                            IJ.showStatus("Computing parameters diameter "+ dist);
                            Objects3DPopulation pop = nucPopList.get(popIndex);
                            if (pop.getNbObjects() > 1) {
                                ArrayUtil allDist = pop.distancesAllClosestCenter();
                                double meanDist = allDist.getMean();
                                double stdDist = allDist.getStdDev();
                                double[] vol = getMeanVolume(pop, false);
                                double[] volActin = getMeanVolume(pop, true);
                                if (popIndex == 0)
                                    outPutGlobalNucleus.write(rootName+"\t"+dist+"\t"+pop.getNbObjects()+"\t"+meanDist+"\t"+stdDist+"\t"+vol[0]+"\t"+vol[1]+"\t"+volActin[0]+"\t"+volActin[1]+"\n");
                                else
                                    outPutGlobalNucleus.write("\t"+dist+"\t"+pop.getNbObjects()+"\t"+meanDist+"\t"+stdDist+"\t"+vol[0]+"\t"+vol[1]+"\t"+volActin[0]+"\t"+volActin[1]+"\n");
                            }
                            else
                                if (popIndex == 0)
                                    outPutGlobalNucleus.write(rootName+"\t"+dist+"\t0\t0\t0\t0\t0\n");
                                else
                                    outPutGlobalNucleus.write("\t"+dist+"\t0\t0\t0\t0\t0\n");
                            outPutNucleus.flush();
                            popIndex++;
                        }
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
