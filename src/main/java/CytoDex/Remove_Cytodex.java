package CytoDex;

//  This plugin extract branches of spherois and compute lengths branching ...


import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.*;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import ij.plugin.RGBStackMerge;
import ij.process.ImageProcessor;
import java.io.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.util.ImageProcessorReader;


public class Remove_Cytodex implements PlugIn {
	
    private final boolean canceled =false;
    private static String imgOutDir;
    private static String fileNameWithOutExt;
    private Calibration cal = new Calibration();
    private int cytoDexRad;
    private int cytoDexCenterZ;
    private int cytoDexCenterX;
    private int cytoDexCenterY;
    static int progressCounter;
    
    /**
     * Read channels
     * @param r
     * @param width
     * @param height
     * @param channel
     * @param name
     * @param series
     * @return 
     */
     private ImagePlus readChannel(ImageProcessorReader r, int width, int height, int channel, String name, String series) {
        ImageStack stack = new ImageStack(width, height);
        for (int n = channel; n < r.getImageCount(); n+=r.getSizeC()) {
            ImageProcessor ip = null;
            try {
                try {
                    ip = r.openProcessors(n)[0];
                } catch (FormatException ex) {
                    Logger.getLogger(Remove_Cytodex.class.getName()).log(Level.SEVERE, null, ex);
                }
            } catch (IOException ex) {
                Logger.getLogger(Remove_Cytodex.class.getName()).log(Level.SEVERE, null, ex);
            }
            IJ.showStatus("reading channel "+channel+" ...");
            IJ.showProgress(n, r.getImageCount());
            stack.addSlice("" + (n + 1), ip);
        }
        ImagePlus imgStack = new ImagePlus(name+" - "+series, stack);
        imgStack.setCalibration(cal);
        return imgStack;
    }   


    /*
    Compute cytoDex center positions and radius
    and shollStep
    */
    private void cytoDexGetCenter(Roi cytoDex) {
        // parameters for center, min, max and step radius 
        if (cytoDex.getBounds().width >= cytoDex.getBounds().height)
            cytoDexRad = cytoDex.getBounds().width/2;
        else
            cytoDexRad = cytoDex.getBounds().height/2;
        cytoDexCenterX = cytoDex.getBounds().x + cytoDex.getBounds().width/2;
        cytoDexCenterY = cytoDex.getBounds().y + cytoDex.getBounds().height/2;
    }
    
   
    
    @Override
    public void run(String arg) {
        try {
            if (canceled) {
                IJ.showMessage(" Pluging canceled");
                return;
            }           
            String imageDir = IJ.getDirectory("Choose Directory Containing CZI Files...");
            if (imageDir == null) return;
            File inDir = new File(imageDir);
            String [] imageFile = inDir.list();
            if (imageFile == null) return;
             // create directory to store images and data
            imgOutDir = imageDir+"Images/";
            File imgTmpDir = new File(imgOutDir);
            if (!imgTmpDir.isDirectory())
                imgTmpDir.mkdir();
            
            // write results headers for cytodex positions and size 
            FileWriter fwCytodex = new FileWriter(imgOutDir + "CytodexResults.xls",false);
            BufferedWriter outputCytodex = new BufferedWriter(fwCytodex);
            outputCytodex.write("Image\tCx\tCy\tCz\tCr\n");
            outputCytodex.flush();
            
            IJ.run("Colors...", "foreground=white background=black selection=yellow");
            // create OME-XML metadata store of the latest schema version
            ServiceFactory factory;
            factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta); 
            for (int i = 0; i < imageFile.length; i++) {
                // for all czi files
                if (imageFile[i].endsWith(".czi")) {
                    String imageName = inDir+ File.separator+imageFile[i];
                    // get file name without extension
                    fileNameWithOutExt = imageFile[i].substring(0,imageFile[i].indexOf(".czi"));
                    reader.setId(imageName);
                    // get number of series
                    int seriesCount = reader.getSeriesCount();
                    
                    // read only last series
                    int series =  seriesCount - 1;
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
                    String seriesName = meta.getImageName(series);         // name of current series
                    // open green cytoDex C0 
                    ImagePlus imgC0 = readChannel(reader, width, height, 0, fileNameWithOutExt, seriesName+"-Cytodex");
                    
                    // Select channel for cytodex auto contrast                    
                    for ( int z = 1; z <= imgC0.getNSlices(); z++) {
                        imgC0.setSlice(z);
                        IJ.run(imgC0,"Enhance Contrast", "saturated=0.35");
                    }
                    imgC0.setSlice(imgC0.getNSlices()/2);
                    imgC0.show();
                    
                    // ask for outline cytodex
                    IJ.setTool("oval");
                    new WaitForUserDialog("Go to the middle of the cytodex and outline the cytodex").show(); 			
                    Roi cytodexRoi = imgC0.getRoi();
                    cytoDexCenterZ = imgC0.getZ();
                    cytoDexGetCenter(cytodexRoi);
                    imgC0.deleteRoi();
                    
                    imgC0.hide();
                    IJ.resetMinAndMax(imgC0);
                    imgC0.updateAndDraw();
                    // open blue nucleus C1 
                    ImagePlus imgC1 = readChannel(reader, width, height, 1, fileNameWithOutExt, seriesName+"-Nucleus");
                    // Save images
                    ImageStack[] imgStack = {imgC0.getImageStack(), imgC1.getImageStack()};
                    ImagePlus stack = new RGBStackMerge().createComposite(imgC0.getWidth(), imgC0.getHeight(), imgC0.getStackSize(), imgStack, false);
                    stack.setCalibration(cal);
                    FileSaver imgC0Save = new FileSaver(stack);
                    imgC0Save.saveAsTiff(imgOutDir+fileNameWithOutExt+".tif");
                    imgC0.close();
                    imgC1.close();
                    stack.close();
                    outputCytodex.write(fileNameWithOutExt+"\t"+cytoDexCenterX+"\t"+cytoDexCenterY+"\t"+cytoDexCenterZ+"\t"+cytoDexRad+"\n");
                    outputCytodex.flush();
                }    
            }
            fwCytodex.close();
            outputCytodex.close();
        } catch (IOException | DependencyException | ServiceException | FormatException ex) {
            Logger.getLogger(Remove_Cytodex.class.getName()).log(Level.SEVERE, null, ex);
        } 
            IJ.showStatus("End of process");
    }
    
}
