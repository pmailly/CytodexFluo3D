package CytoDex;

//  This plugin extract branches of spherois and compute lengths branching ...

import static Tools.Cytodex_tools.*;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.io.FileSaver;
import ij.plugin.Duplicator;
import ij.plugin.PlugIn;
import ij.plugin.frame.RoiManager;
import ij.process.ImageStatistics;
import java.io.*;
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
import mcib3d.image3d.ImageFloat;


public class CytodexDSred_Fluo3D implements PlugIn {
    
    private final String positionName = "Experiment|AcquisitionBlock|RegionsSetup|SampleHolder|SingleTileRegion|Name #";
    public Roi roi = null;
    
     
    @Override
    public void run(String arg) {
        try {
            FileWriter fwAnalyze = null;
            if (canceled) {
                IJ.showMessage(" Pluging canceled");
                return;
            }
            String imageDir = IJ.getDirectory("Choose Directory Containing CZI or TIF Files");
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
                    
                    // open red cytoDex bead C0 
                    ImagePlus imgC0 = readChannel(reader, width, height, 0, fileNameWithOutExt, seriesName+"_CytodexBead");
                    // 
                    // detect and remove bead
                    findSpheroid(imgC0, roi);
                    // no bead
                    if (cytoDexCenterX == 0) {
                        imgC0.changes = false;
                        imgC0.close();
                        imgC0.flush();
                        if (imageNum == 1)
                            imageNum = 0;
                    }
                    else {
                        removeSpheroid(imgC0, 60000);
                        FileSaver imgSave = new FileSaver(imgC0);
                        imgSave.saveAsTiff(imgOutDir+fileNameWithOutExt+"_"+seriesName+"_CytoDexBead.tif");
                        imgC0.changes = false;
                        imgC0.close();
                        imgC0.flush();
                        
                        

                        // open green cytoDex capillaries C1 
                        ImagePlus imgC1 = readChannel(reader, width, height, 1, fileNameWithOutExt, seriesName+"_Cytodex");
                        IJ.run(imgC1, "Median...", "radius=2 stack");
                        
                        // Find with DOG
                        ImagePlus imgC1DOG = vesselDOG(imgC1, seriesName, roi, imageNum);
                        // skeletonize
                        ImagePlus imgSkelDOG = new Duplicator().run(imgC1DOG,1,imgC1DOG.getNSlices());
                        imgSkelDOG.setTitle(fileNameWithOutExt +"_"+seriesName+ "_Skel.tif");
                        IJ.run(imgSkelDOG, "Skeletonize (2D/3D)", "");
                        
                        // Check if no branches
                        imgSkelDOG.setSlice(cytoDexCenterZ);
                        ImageStatistics statsDOG = ImageStatistics.getStatistics(imgSkelDOG.getProcessor(), ImageStatistics.MIN_MAX,imgSkelDOG.getCalibration());

                        // no branches
                        if (statsDOG.max == statsDOG.min ) {
                            // write skeleton data with zero
                            outputAnalyze.write(imgOutDir+fileNameWithOutExt+"_"+seriesName+"\t0\t0\t0\t0\t0\n");
                            outputAnalyze.flush();
                            imgC1DOG.changes = false;
                            imgC1DOG.close();
                            imgC1DOG.flush();

                        } 
                        else { 
                            // open blue nucleus C2
                            ImagePlus imgC2 = readChannel(reader, width, height, 2, fileNameWithOutExt, seriesName+"_Nucleus");
                            // find nucleus
                            findNucleus(imgC2, roi, imageNum);
                            imgC2.changes = false;
                            imgC2.close();
                            imgC2.flush(); 
                            
                            // compute intersections
                            //intersectionAnalysis(imgSkelDOG); 
                            
                            // Analyze skeleton
                            analyzeSkel(imgSkelDOG,outputAnalyze);
                            // compute image map
                            ImageFloat imgMapDOG = localThickness3D(imgC1DOG);                           
                            // compute mean diameter and intersections from concentric spheres
                            diameterAnalysis(imgSkelDOG, imgMapDOG);
                        }
                        imgSkelDOG.changes=false;
                        imgSkelDOG.flush();
                        imgSkelDOG.close();
                        imgC1.changes = false;
                        imgC1.close();
                        imgC1.flush();
                        
                    }                    
                }
                IJ.run("Close All");
            }   
            fwAnalyze.close();
            outputAnalyze.close();
            IJ.showStatus("End of process");
        } catch (IOException | DependencyException | ServiceException | FormatException ex) {
            Logger.getLogger(CytodexDSred_Fluo3D.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
