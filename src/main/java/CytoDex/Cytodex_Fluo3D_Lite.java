package CytoDex;

//  This plugin extract branches of spherois and compute lengths branching ...

import static Tools.Cytodex_tools.*;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.WaitForUserDialog;
import ij.io.FileSaver;
import ij.plugin.Duplicator;
import ij.plugin.PlugIn;
import ij.plugin.RGBStackMerge;
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
import mcib3d.geom.Objects3DPopulation;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;
import org.apache.commons.io.FilenameUtils;
import sc.fiji.analyzeSkeleton.AnalyzeSkeleton_;
import sc.fiji.analyzeSkeleton.Edge;
import sc.fiji.analyzeSkeleton.SkeletonResult;


public class Cytodex_Fluo3D_Lite implements PlugIn {
    
    private boolean canceled = false;
    private String imageDir = "";
    private BufferedWriter outPutAnalyze;
    private final double nucMinVol = 150;
    private final double nucMaxVol = Double.MAX_VALUE;
    private String thresholdMethod = "Default";
    private int minBranchVol = 300;
    
    
    
    // Calculate lenght of branches after skeletonize
    public static double[] analyzeSkeleton (ImagePlus img, String output) {
        IJ.run(img, "Skeletonize (2D/3D)", "");
	String imgTitle = img.getTitle();
        AnalyzeSkeleton_ analyzeSkeleton = new AnalyzeSkeleton_();
        AnalyzeSkeleton_.calculateShortestPath = true;
        analyzeSkeleton.setup("",img);
        IJ.showStatus("Analyze skeleton...");
        SkeletonResult skeletonResults = analyzeSkeleton.run(AnalyzeSkeleton_.NONE, false, true, null, true, false);
        ImageStack imgStackLab = analyzeSkeleton.getLabeledSkeletons();
        skeletonResults = pruneEndBranches(img.getStack(), imgStackLab, skeletonResults, smallBranch);
//        // remove small branches
//        IJ.showStatus("Removing small branches...");
//        for (int i = 0; i < 5; i++) {
//            removeSmallBranches(img, imgStackLab, skeletonResults);
//            analyzeSkeleton.setup("",img);
//            skeletonResults = analyzeSkeleton.run(AnalyzeSkeleton_.NONE, false, true, null, true, false);
//        }

        //  compute parameters for each skeleton
        IJ.showStatus("Computing parameters for each skeleton ...");
        int[] branchNumbers = skeletonResults.getBranches();
        double[] branchLengths = skeletonResults.getAverageBranchLength();
        int branches = 0;
        double totalBranchLength = 0;
        int skelNumber = skeletonResults.getGraph().length;
        for (int i = 0; i < skelNumber; i++) {
            totalBranchLength += branchLengths[i] * branchNumbers[i];
            branches += branchNumbers[i];
        }
        // other way to have total lenght
        double totalEdgeLength = 0;
        for (int i = 0; i < skelNumber; i++) {
            ArrayList<Edge> listEdges;
            listEdges = skeletonResults.getGraph()[i].getEdges();
            for (int e = 0; e < listEdges.size(); e++) {
                totalEdgeLength += listEdges.get(e).getLength();
            }
        }
        //System.out.println("branch length = "+branchLength+ ", total length = "+totalLength);
        double[] skeletonParams = {skelNumber, branches, totalBranchLength, totalEdgeLength};
        IJ.showStatus("Saving labeled skeleton ...");
        ImagePlus imgLab = new ImagePlus(imgTitle+"_LabelledSkel.tif", imgStackLab);
        ImagePlus imgLabProj = doZProjection(imgLab);
        IJ.run(imgLabProj, "3-3-2 RGB", "");
        imgLabProj.setCalibration(cal);
        imgLab.changes = false;
        imgLab.close();
        imgLab.flush();   
        FileSaver imgSave = new FileSaver(imgLabProj);
        imgSave.saveAsTiff(output+imgTitle+"_LabelledSkel.tif");
        imgLabProj.close();
        return(skeletonParams);
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
    
    
    
    private boolean dialog() {
        GenericDialogPlus gd = new GenericDialogPlus("Channels");
        gd.addDirectoryField("Choose Directory Containing CZI", imageDir);
        gd.addNumericField("Min branch volume (µm) :", minBranchVol, 1);
        gd.showDialog();
        if (gd.wasCanceled())
            canceled = true;
        imageDir = gd.getNextString()+ File.separator;
        minBranchVol = (int)gd.getNextNumber();
        return(canceled);
    }
        
        // Threshold images and fill holes
    public static void threshold(ImagePlus img, String thMed, boolean fill, boolean calcul) {
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
    
    @Override
    public void run(String arg) {
        String imageExt = "czi";
        try {
            if (dialog()) {
                IJ.showMessage(" Pluging canceled");
                return;
            }
            File inDir = new File(imageDir);
            String[] imageFile = inDir.list();
            if (imageFile == null) {
                System.out.println("No Image found in "+imageDir);
                return;
            }
            // create output folder
            String outDirResults = imageDir + File.separator+ "Results"+ File.separator;
            File outDir = new File(outDirResults);
            if (!Files.exists(Paths.get(outDirResults))) {
                outDir.mkdir();
            }

            // create OME-XML metadata store of the latest schema version
            ServiceFactory factory;
            factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            Arrays.sort(imageFile);
            int imageNum = 0; 
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
                        outPutAnalyze.write("Image Name\tSkeleton number\tBranch number\tTotal Branch Length\tTotal Edge Length\tNucleus number\n");
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
                    IJ.run(imgActin, "Subtract Background...", "rolling=25 stack");
                    median_filter(imgActin, 1.5);
                    IJ.run(imgActin, "Laplacian of Gaussian", "sigma=6 scale_normalised negate stack");
                    if (imageNum == 1)
                        dialogThreshod(imgActin);
                    threshold(imgActin, thresholdMethod, false, false);

                    Objects3DPopulation actinPop = new Objects3DPopulation(getPopFromImage(imgActin).getObjectsWithinVolume(minBranchVol, Double.MAX_VALUE, false));
                    System.out.println("Actin population = "+ actinPop.getNbObjects());

                    // Draw actin population
                    ImageHandler imhActinObjects = ImageInt.wrap(imgActin).createSameDimensions();
                    actinPop.draw(imhActinObjects, 255);
                    ImagePlus imgSkel = new Duplicator().run(imhActinObjects.getImagePlus());
                    
                    // Find actin network morphology
                    // Skeletonize
                    imgSkel.setCalibration(cal);
                    imgSkel.setTitle(imgActin.getTitle());
                    double[] skeletonParams = analyzeSkeleton(imgSkel, outDirResults);
                    imgSkel.close();
                    
                    
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
                    Objects3DPopulation nucPop = new Objects3DPopulation(getPopFromImage(imgNuc).getObjectsWithinVolume(nucMinVol, nucMaxVol, true));
                    int nucNumber = nucPop.getNbObjects();
                    System.out.println("Nucleus population = "+ nucNumber);
                    
                    outPutAnalyze.write(rootName+"\t"+skeletonParams[0]+"\t"+skeletonParams[1]+"\t"+skeletonParams[2]+"\t"+skeletonParams[3]+"\t"+nucNumber+"\n");
                    
                    // Save objects image
                    System.out.println("SavingObject population ...");
                    ImageHandler imhNucObjects = ImageInt.wrap(imgActin).createSameDimensions();
                    nucPop.draw(imhNucObjects, 255);
                    ImagePlus[] imgColors = {imhActinObjects.getImagePlus(), null, imhNucObjects.getImagePlus(), imgActin};
                    ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, false);
                    imgObjects.setCalibration(cal);
                    IJ.run(imgObjects, "Enhance Contrast", "saturated=0.35");
                    FileSaver ImgObjectsFile = new FileSaver(imgObjects);
                    ImgObjectsFile.saveAsTiff(outDirResults + rootName + "_" + seriesName + "_Objects.tif");
                    imhActinObjects.closeImagePlus();
                    imgObjects.close();
                    imhNucObjects.closeImagePlus();
                    imgActin.close();
                    imgNuc.close();
                    IJ.run("Collect Garbage", "");
                }                    
            } 
            outPutAnalyze.close();
            IJ.showStatus("End of process");
        } catch (IOException | DependencyException | ServiceException | FormatException ex) {
            Logger.getLogger(Cytodex_Fluo3D_Lite.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
