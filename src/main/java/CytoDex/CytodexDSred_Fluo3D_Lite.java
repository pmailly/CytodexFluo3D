package CytoDex;

//  This plugin extract branches of spherois and compute lengths branching ...

import static Tools.Cytodex_tools.*;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.EllipseRoi;
import ij.gui.Roi;
import ij.gui.WaitForUserDialog;
import ij.io.FileSaver;
import ij.plugin.Duplicator;
import ij.plugin.PlugIn;
import ij.plugin.RGBStackMerge;
import ij.plugin.frame.RoiManager;
import ij.process.ImageStatistics;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom.Object3D;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;


public class CytodexDSred_Fluo3D_Lite implements PlugIn {
    
    private final String positionName = "Experiment|AcquisitionBlock|RegionsSetup|SampleHolder|SingleTileRegion|Name #";
    public Roi roi = null;
    private boolean canceled = false;
    private String imageDir = "";
    int imgCh0, imgCh1, imgCh2;
    // remove the bead do not calculate skeleton
    private boolean beadOnly = false;
    private boolean cropZ = false;
    private String cropDir = "";
    private String autoThreshold_Method = "Default";
    
    
    
    private boolean dialog() {
        String[] channels = {"0", "1", "2"};
        String[] structures = {"Bead", "Dapi", "Cytodex"};
        String[] dir = {"Top", "Bottom"};
        GenericDialogPlus gd = new GenericDialogPlus("Channels");
        gd.addDirectoryField("Choose Directory Containing CZI or TIF Files", imageDir);
        for (int i = 0; i < structures.length; i++){
            gd.addChoice("Channel for "+structures[i]+ ": ", channels, channels[0]);
        }
        gd.addCheckbox(" Crop stack", cropZ);
        gd.addChoice("Direction : ", dir, dir[0]);
        gd.addCheckbox(" Do not calculate skeleton", beadOnly);
        gd.showDialog();
        if (gd.wasCanceled())
            canceled = true;
        imageDir = gd.getNextString()+ File.separator;
        imgCh0 = gd.getNextChoiceIndex();
        imgCh1 = gd.getNextChoiceIndex();
        imgCh2 = gd.getNextChoiceIndex();
        cropZ = gd.getNextBoolean();
        cropDir = dir[gd.getNextChoiceIndex()];
        beadOnly = gd.getNextBoolean();
        return(canceled);
    }
    
    
   /*
    * 
     * @param obj1
     * @param obj2
     * remove obj1 if 50% coloc with obj2
     * 
    */
    private void filterNucleus(Objects3DPopulation obj1, Objects3DPopulation obj2, ImagePlus img) {
        Objects3DPopulation nucleusInside = new Objects3DPopulation();
        double coloc = 0;
        int removedNucleus = 0;
        IJ.showStatus("Removing nucleus inside vessels ...");
        for (int ob1 = 0; ob1 < obj1.getNbObjects(); ob1++) {
            Object3D obj3D1 = obj1.getObject(ob1);
            for (int ob2 = 0; ob2 < obj2.getNbObjects(); ob2++) { 
                Object3D obj3D2 = obj2.getObject(ob2);
                // find coloc
                coloc = obj3D1.pcColoc(obj3D2);
                
                if (coloc > 50) {
                    nucleusInside.addObject(obj3D1);
                    obj1.removeObject(obj3D1);
                    ob1--;
                    removedNucleus++;
                }
            }
        }
        // create object image
        ImageHandler imgObjects = ImageInt.wrap(img).createSameDimensions();
        imgObjects.set332RGBLut();
        nucleusInside.draw(imgObjects, 225);
        imgObjects.setCalibration(cal);
        // save image
        FileSaver ObjectsFile = new FileSaver(imgObjects.getImagePlus());
        ObjectsFile.saveAsTiffStack(imgOutDir +File.separator+ fileNameWithOutExt + "_nucleusInside.tif");
        imgObjects.closeImagePlus();
        System.out.println(removedNucleus+" nucleus inside vessels");
    }
    
    /** Find min distance from nucleus centroid to vessel border
     *  Exclude nucleus if distance > 10 µm and < 0.7 µm
     * @param nucleus
     * @param vessels
     * @return 
     */
    private ArrayList<Double> distanceToBorder (Objects3DPopulation nucleus, Objects3DPopulation vessels, ImagePlus img) {
        ArrayList<Double> dist = new ArrayList<>();
        Object3D nucObj;
        Object3D vesselObj;
        Objects3DPopulation recructedNucleus = new Objects3DPopulation();
        double minDist = 0.7;
        double maxDist = 12;
        double distCenter;
        boolean limitOk;
        for (int ob1 = 0; ob1 < nucleus.getNbObjects(); ob1++) {
            maxDist = 12;
            limitOk = false;
            nucObj = nucleus.getObject(ob1);
            for (int ob2 = 0; ob2 < vessels.getNbObjects(); ob2++) {
                vesselObj = vessels.getObject(ob2);
                distCenter = nucObj.distCenterBorderUnit(vesselObj);
                if ((distCenter < maxDist) && (distCenter > minDist)) {
                    maxDist = distCenter;
                    limitOk = true;
                }
            }
            if (limitOk) {
                dist.add(maxDist);
                recructedNucleus.addObject(nucObj);
            }
        }            
        // create object image
        ImageHandler imgObjects = ImageInt.wrap(img).createSameDimensions();
        imgObjects.set332RGBLut();
        recructedNucleus.draw(imgObjects, 253);
        imgObjects.setCalibration(cal);
        // save image
        FileSaver ObjectsFile = new FileSaver(imgObjects.getImagePlus());
        ObjectsFile.saveAsTiffStack(imgOutDir +File.separator+ fileNameWithOutExt + "_recrutedNucleus.tif");
        imgObjects.closeImagePlus();
        return(dist);
    }
    
    // Find vessels with DOG
    public static ImagePlus vesselDOG(ImagePlus img, String seriesName, Roi roi, int first) {
        double minThreshold_Vessel = 0;
        double maxThreshold_Vessel = 0;
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
            new WaitForUserDialog("Select Threshold method and press OK").show();
            minThreshold_Vessel = img.getProcessor().getMinThreshold();
            maxThreshold_Vessel = img.getProcessor().getMaxThreshold();
            WindowManager.getWindow("Threshold").setVisible(false);
            EllipseRoi roi_Circle = new EllipseRoi(cytoDexCenterX - cytoDexRad, cytoDexCenterY - cytoDexRad, 
                    cytoDexCenterX + cytoDexRad, cytoDexCenterY + cytoDexRad, 1);
            IJ.setTool("oval");
            img.setRoi(roi_Circle);
            img.updateAndDraw();
            new WaitForUserDialog("Resize circle to define bead diameter and press OK").show();
            double cytoDexRadNew = roi_Circle.getFeretsDiameter()/2;
            cytoDexEnlarge = cytoDexRadNew - cytoDexRad;
            System.out.println("New cyto radius = "+cytoDexRadNew+ " enlarge thickness = "+cytoDexEnlarge);
            img.deleteRoi();   
            img.hide(); 
        }
        Prefs.blackBackground = false;
        IJ.setThreshold(img,minThreshold_Vessel, maxThreshold_Vessel);
        IJ.run(img, "Convert to Mask","background=Dark black");
            
        //IJ.run(img, "Options...", "iterations=20 count=1 do=Close stack");
        removeSpheroid(img, 0);
        if (roi != null) {
            img.setRoi(roi);
            IJ.run("Colors...", "foreground=white background=black selection=yellow");
            IJ.run(img, "Clear Outside","stack");
            img.deleteRoi();
        }
//        FileSaver imgObjectsSave = new FileSaver(img);
//        imgObjectsSave.saveAsTiff(imgOutDir+fileNameWithOutExt+"_VesselMask.tif");
        return(img);
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
    
    
    
    @Override
    public void run(String arg) {
        try {
            if (dialog()) {
                IJ.showMessage(" Pluging canceled");
                return;
            }
            FileWriter fwAnalyze = null;
            if (imageDir == null) return;
            File inDir = new File(imageDir);
            String [] imageFile = inDir.list();
            Arrays.sort(imageFile);
            if (imageFile == null) return;
            imgOutDir = imageDir+"Results/";
            File imgTmpDir = new File(imgOutDir);
            if (!imgTmpDir.isDirectory())
                imgTmpDir.mkdir();
            fwAnalyze = new FileWriter(imgOutDir + "Analyze_skeleton_results.xls",false);
            BufferedWriter outputAnalyze = new BufferedWriter(fwAnalyze);
            outputAnalyze.write("Image\t#Skeletons\t#Branches\t#Junctions\t#End points\tBranch length\n");
            outputAnalyze.flush();
            IJ.run("Colors...", "foreground=white background=black selection=yellow");
            // create OME-XML metadata store of the latest schema version
            ServiceFactory factory;
            factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta); 
            int imageNum = 0;
            for (int i = 0; i < imageFile.length; i++) {
                // for all czi or tif files
                if ((imageFile[i].endsWith(".czi")) || (imageFile[i].endsWith(".tif"))) {
                    String tilePosition;
                    String imageName = inDir+ File.separator+imageFile[i];
                    String seriesName = "";
                    
                    // read only first series
                    int series = 0; 
                    reader.setId(imageName);
                    reader.setSeries(series);
                    int width = reader.getSizeX();
                    int height = reader.getSizeY();
                    double sx = meta.getPixelsPhysicalSizeX(series).value().doubleValue();
                    double sy = meta.getPixelsPhysicalSizeY(series).value().doubleValue();
                    double sz = meta.getPixelsPhysicalSizeZ(series).value().doubleValue(); 
                    cal.pixelWidth = sx;
                    cal.pixelHeight = sy;
                    cal.pixelDepth = sz;
                    cal.setUnit("microns");
                    //System.out.println(cal.pixelWidth+", "+cal.pixelDepth);
                    // get file name without extension
                    // find number in parentheses to identify position name (well and cytodex number) in metadata
                    if (imageFile[i].endsWith(".czi")) {
                        fileNameWithOutExt = imageFile[i].substring(0,imageFile[i].indexOf(".czi"));
                        Matcher m = Pattern.compile("\\(([0-9]+)\\)").matcher(fileNameWithOutExt);
                        if (m.find()) {
                            tilePosition = m.group(1);  
                            int digit = 0;
                            // Find tile position name increase digit
                            Object metaScan;
                            do {
                                digit++;
                                metaScan = reader.getMetadataValue(positionName + IJ.pad(Integer.valueOf(tilePosition), digit));
                            } while (metaScan == null);   
                            seriesName = metaScan.toString();
                        }    
                    }                    
                    else {
                        fileNameWithOutExt = imageFile[i].substring(0,imageFile[i].indexOf(".tif"));
                    }
                    imageNum++;
                    // if no position name
                    if (seriesName.isEmpty()) 
                        seriesName = meta.getImageName(series);         // name of current series
                    
                    // check if roi file exist
                    if (new File(imageDir+fileNameWithOutExt+".roi").exists()) {
                        RoiManager rm = new RoiManager(false);
                        rm.runCommand("Open", imageDir+fileNameWithOutExt+".roi");
                        roi = rm.getRoi(0);
                        System.out.println("Roi found");
                        rm.close();
                    }
                    
                    // open red cytoDex bead  
                    ImagePlus imgCBead = readChannel(reader, width, height, imgCh0, fileNameWithOutExt, seriesName+"_CytodexBead");
                    // detect and remove bead
                    findSpheroid(imgCBead, roi);
                    // no bead
                    if (cytoDexCenterX == 0) {
                        imgCBead.changes = false;
                        imgCBead.close();
                        imgCBead.flush();
                        if (imageNum == 1)
                            imageNum = 0;
                    }
                    else {
                        imgCBead.changes = false;
                        imgCBead.close();
                        imgCBead.flush();
                        
                        // open green cytoDex capillaries 
                        ImagePlus imgCytodex = readChannel(reader, width, height, imgCh2, fileNameWithOutExt, seriesName+"_Vessel");
                        removeSpheroid(imgCytodex, 0);
                        
                        // crop stack in Z
                        if (cropZ)
                            cropStack(imgCytodex, cropDir);
                        
                       // open green cytoDex nucleus 
                        ImagePlus imgCNuc = readChannel(reader, width, height, imgCh1, fileNameWithOutExt, seriesName+"_Nucleus");
                        removeSpheroid(imgCNuc, 0);
                        // Crop stack in Z
                        if (cropZ)
                            cropStack(imgCNuc, cropDir);

                        ImagePlus[] images = new ImagePlus[2];
                        images[0] = imgCNuc.duplicate();
                        images[1] = imgCytodex.duplicate();
                        ImagePlus merged = new RGBStackMerge().mergeHyperstacks(images, false);
                        merged.setCalibration(cal);
                        new FileSaver(merged).saveAsTiff(imgOutDir+fileNameWithOutExt+"_"+seriesName+"_CytoDexBead.tif");
                        merged.close();
                        
                        IJ.run(imgCytodex, "Median...", "radius=2 stack");
                        
                        // Find vessel with DOG
                        ImagePlus imgCytodexDOG = vesselDOG(imgCytodex, seriesName, roi, imageNum);                      
                        imgCytodexDOG.setTitle(fileNameWithOutExt +"_"+seriesName+ "_Mask.tif");
                        IJ.saveAs(imgCytodexDOG, "Tiff", imgOutDir + fileNameWithOutExt + "_mask.tif");
                        // skeletonize
                        if (!beadOnly) {
                            ImagePlus imgSkelDOG = new Duplicator().run(imgCytodexDOG,1,imgCytodexDOG.getNSlices());
                            imgSkelDOG.setTitle(fileNameWithOutExt +"_"+seriesName+ "_Skel.tif");
                            IJ.run(imgSkelDOG, "Skeletonize (2D/3D)", "");                    
                            // Check if no branches
                            imgCytodexDOG.setSlice(cytoDexCenterZ);
                            ImageStatistics statsDOG = ImageStatistics.getStatistics(imgCytodexDOG.getProcessor(), ImageStatistics.MIN_MAX,imgCytodexDOG.getCalibration());

                            // no branches
                            if (statsDOG.max == statsDOG.min ) {
                                System.out.println("No vessel detected");
                                // write skeleton data with zeroSystem.out.println("No vessel detected");
                                outputAnalyze.write(imgOutDir+fileNameWithOutExt+"_"+seriesName+"\t0\t0\t0\t0\t0\n");
                                outputAnalyze.flush();
                                imgCytodexDOG.changes = false;
                                imgCytodexDOG.close();
                                imgCytodexDOG.flush();
                            } 
                            else {
                                // Analyze skeleton
                                analyzeSkel(imgSkelDOG,outputAnalyze, smallBranch);
                                // compute image map
                                //ImageFloat imgMapDOG = localThickness3D(imgC1DOG);                           
                                // compute mean diameter and intersections from concentric spheres
                                //diameterAnalysis(imgSkelDOG, imgMapDOG);

                                // find nucleus
                                Objects3DPopulation nucleus = findNucleus(imgCNuc, null, imageNum);
                                // convert vessel to 3DObjects
                                Objects3DPopulation vessels = getPopFromImage(imgCytodexDOG);
                                // remove nucleus with 50% coloc with vessels
                                filterNucleus(nucleus, vessels, imgCytodexDOG);
                                // compute distance between nucleus and vessel
                                ArrayList<Double> dist = distanceToBorder(nucleus, vessels, imgCytodexDOG);
                                fwAnalyze = new FileWriter(imgOutDir + fileNameWithOutExt+"_"+seriesName+"_Distresults.xls",false);
                                BufferedWriter outputDist = new BufferedWriter(fwAnalyze);
                                outputDist.write("Image\tDistance\n");
                                outputDist.flush();
                                for (int d = 0; d < dist.size(); d++) {
                                    outputDist.write(dist.get(d).toString()+"\n");
                                    outputDist.flush();
                                }
                                outputDist.close();
                                imgSkelDOG.changes=false;
                                imgSkelDOG.flush();
                                imgSkelDOG.close();
                            }
                            
                            imgCNuc.changes = false;
                            imgCNuc.close();
                            imgCNuc.flush();  
                            imgCytodexDOG.changes=false;
                            imgCytodexDOG.flush();
                            imgCytodexDOG.close();
                            
                        }
                        
                        imgCBead.changes = false;
                        imgCBead.close();
                        imgCBead.flush();
                    }                    
                }
                IJ.run("Close All");
            }   
            IJ.showStatus("End of process");
        } catch (IOException | DependencyException | ServiceException | FormatException ex) {
            Logger.getLogger(CytodexDSred_Fluo3D_Lite.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
