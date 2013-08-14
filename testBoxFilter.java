import java.awt.image.BufferedImage;
import java.awt.image.Raster;
import java.awt.image.RenderedImage;
import java.awt.image.WritableRaster;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;


public class testBoxFilter {
    
    public static void main(String[] args) {
        
        //read in images TESTING
        RenderedImage NadirImage = null;
        try {
            NadirImage = ImageIO.read(new File ("test_nadir.png"));
            } catch(Exception e) {}
            
       //convert to arrays TESTING
       int imageWidth = NadirImage.getWidth();
       int imageHeight = NadirImage.getHeight(); 
       Raster NadirRaster = NadirImage.getData();
       
       int[][] nadirArray = new int[imageWidth][imageHeight];
       
       for (int i = 0; i < imageWidth; i++) {
           for (int j = 0; j < imageHeight; j++) {
               
               nadirArray[i][j] = NadirRaster.getSample(i, j, 0);
           }
       }
       
       int[] inPixels = new int[imageWidth*imageHeight];
       int[] outPixels = new int[imageWidth*imageHeight];
       int[] smoothedPixels = new int[imageHeight*imageWidth];
       int smoothingRadius = 2; 
       
       //convert pre cost array to 1d
       int count = 0;
       for (int x = 0; x < imageWidth; x++) {
           for (int y = 0; y < imageHeight; y++) {
               
               inPixels[count] = nadirArray[x][y];
               count++;
           }
       }
       
       //box filter average
       blur(inPixels, outPixels, imageWidth, imageHeight, smoothingRadius);
       blur(outPixels, smoothedPixels, imageHeight, imageWidth, smoothingRadius);
                      
       //loop over smoothed precost array and update costs xdis and ydis
       count = 0;
       for (int x = 0; x < imageWidth; x++) {
           for (int y = 0; y < imageHeight; y++) {
              
              nadirArray[x][y] = smoothedPixels[count];
              count++;
              
               } 
           } 
            
       
       //save output
       
       WritableRaster nadirOutput = Raster.createWritableRaster(NadirImage.getSampleModel(), null); 
       
       for (int x = 0; x < imageHeight; x++) {
           for (int y = 0; y < imageWidth; y++) {
               nadirOutput.setSample(x, y, 0, nadirArray[x][y]);
           
           }
       
       }
      
       //write out to BufferedImage TESTING
       BufferedImage yDisBI = new BufferedImage(NadirImage.getColorModel(), nadirOutput, false, null);
      
       File smoothedNadirFile = new File("smoothedNadir.png");
       try {
           ImageIO.write(yDisBI, "png", smoothedNadirFile);
       } catch (IOException e) {
           e.printStackTrace();
       }
            
    }
       
    
       public static void blur( int[] in, int[] out, int width, int height, int radius ) {
           int widthMinus1 = width-1; 
           int diameter = 2*radius + 1;
      
           int inIndex = 0;  
           
           //loop over the rows
           for (int y = 0; y < height; y++) {  
               int outIndex = y;   
               int currentSum = 0;   

               //initialise the moving window
               for (int i = -radius; i <= radius; i++) { 
                   int inputPixel = in[(int) (inIndex + clamp(i, 0, width-1))];    
                   currentSum += inputPixel;
               }
               
               //move along the columns
               for ( int x = 0; x < width; x++ ) {
                   out[outIndex] = currentSum/diameter;  
                   
                   //edge check
                   //selecting the next pixel to come into the window
                   int incoming = x+radius+1;   
                   if (incoming > widthMinus1)
                       incoming = widthMinus1;
                   //removing the last pixel in the window
                   int outgoing = x-radius;     
                   if (outgoing < 0)
                       outgoing = 0;
                   
                   //get incoming and outgoing pixels
                   int incomingPixel = in[inIndex+incoming];  
                   int outgoingPixel = in[inIndex+outgoing]; 
                   
                   //adjust the sum
                   currentSum += (incomingPixel)-(outgoingPixel); 
                   
                   //do the transpose
                   outIndex += height;
               }
               //do the transpose
               inIndex += width;
           }
       }
       
       //edge clamper
       public static float clamp(float x, float a, float b) {
           return (x < a) ? a : (x > b) ? b : x;
       }
       

}
