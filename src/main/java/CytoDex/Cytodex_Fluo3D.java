package CytoDex;

//  This plugin extract branches of spheroids and compute lengths branching ...


import static Tools.Cytodex_tools.*;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.io.Opener;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.ChannelSplitter;
import ij.plugin.Duplicator;
import ij.plugin.PlugIn;
import static ij.plugin.filter.Analyzer.getResultsTable;
import ij.process.AutoThresholder;
import ij.process.ImageStatistics;
import java.io.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import mcib3d.image3d.ImageHandler;


public class Cytodex_Fluo3D implements PlugIn {
    
    public Roi roi = null;
    private Calibration cal = new Calibration();
    private double smallBranch = 25;
    private int interPruning =5;
    

     /**
      * Read cytodex positions and size
      */
     
     void readCytodexPositions(String path, String image) {
         Opener.openResultsTable(path);
         ResultsTable table = getResultsTable();
         for (int n = 0; n < table.size(); n++) {
             if(table.getStringValue("Image", n).equals(image)) {
                 cytoDexCenterX = (int)table.getValue("Cx", n);
                 cytoDexCenterY = (int)table.getValue("Cy", n);
                 cytoDexCenterZ = (int)table.getValue("Cz", n);
                 cytoDexRad = (int)table.getValue("Cr", n);
             }
         }
         IJ.selectWindow("Results");
         IJ.run("Close");
     }
     
    @Override
    public void run(String arg) {
        FileWriter fwAnalyze = null;
        try {
            if (canceled) {
                IJ.showMessage(" Pluging canceled");
                return;
            }
            String imageDir = IJ.getDirectory("Choose Directory Containing Tif Files");
            if (imageDir == null) return;
            File inDir = new File(imageDir);
            String [] imageFile = inDir.list();
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
            int imageNum = 0;
            for (int i = 0; i < imageFile.length; i++) {
                // for all tif files
                if (imageFile[i].endsWith(".tif")) {
                    fileNameWithOutExt = imageFile[i].substring(0,imageFile[i].indexOf(".tif"));
                    // read image
                    Opener imgOpen = new Opener();
                    ImagePlus imgOrg = imgOpen.openImage(imageDir, imageFile[i]);
                    cal = imgOrg.getCalibration();
                    // Extract channels
                    ImagePlus[] imgC = ChannelSplitter.split(imgOrg);
                    imgOrg.close();
                    imgOrg.flush();
                    // Read cytodex positions and size for this image
                    readCytodexPositions(imageDir+"CytodexResults.xls",fileNameWithOutExt);
                    //median filter
                    IJ.run(imgC[0], "Median...", "radius=2 stack");
                    // check if roi file exist
                    if (new File(fileNameWithOutExt+".roi").exists()) {
                        IJ.open(inDir+fileNameWithOutExt+".roi"); 
                        roi = imgC[0].getRoi();
                    }
                    imageNum++;
                    // Find with DOG
                    ImagePlus imgC1DOG = vesselDOG(imgC[0], fileNameWithOutExt, roi, imageNum);
                    imgC[0].flush();
                    imgC[0].close();
                    // skeletonize
                    ImagePlus imgSkelDOG = new Duplicator().run(imgC1DOG,1,imgC1DOG.getNSlices());
                    imgSkelDOG.setTitle(fileNameWithOutExt + "_Skel.tif");
                    IJ.run(imgSkelDOG, "Skeletonize (2D/3D)", "");

                    // Check if no branches
                    imgSkelDOG.setSlice(cytoDexCenterZ);
                    ImageStatistics statsDOG = ImageStatistics.getStatistics(imgSkelDOG.getProcessor(), ImageStatistics.MIN_MAX,imgSkelDOG.getCalibration());

                    // no branches
                    if (statsDOG.max == statsDOG.min ) {
                        // write skeleton data with zero
                        outputAnalyze.write(imgOutDir+fileNameWithOutExt+"\t0\t0\t0\t0\t0\n");
                        outputAnalyze.flush();
                        imgC1DOG.changes = false;
                        imgC1DOG.close();
                        imgC1DOG.flush();

                    } 
                    else {  
                        // Analyze skeleton
                        analyzeSkel(imgSkelDOG,outputAnalyze, imgOutDir, smallBranch, interPruning);

                        // compute image map
                        ImagePlus imgMapDOG = localThickness3D(imgC1DOG);
                        
                        // compute intersections
                        intersectionAnalysis(imgSkelDOG);

                        // compute mean diameter and intersections from concentric spheres
                        diameterAnalysis(imgSkelDOG, ImageHandler.wrap(imgMapDOG));

                        // find nucleus
                        removeSpheroid(imgC[1], 0);
                        findNucleus(imgC[1], roi,imageNum);
                        imgC[1].changes = false;
                        imgC[1].close();
                        imgC[1].flush();  
                    }
                    imgSkelDOG.changes=false;
                    imgSkelDOG.flush();
                    imgSkelDOG.close();
                }
            }   fwAnalyze.close();
            outputAnalyze.close();
            IJ.showStatus("End of process");
        } catch (IOException ex) {
            Logger.getLogger(Cytodex_Fluo3D.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                fwAnalyze.close();
            } catch (IOException ex) {
                Logger.getLogger(Cytodex_Fluo3D.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }
    
}
