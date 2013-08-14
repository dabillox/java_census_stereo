import java.awt.image.BufferedImage;
import java.awt.image.Raster;
import java.awt.image.RenderedImage;
import java.awt.image.WritableRaster;
import java.awt.image.renderable.ParameterBlock;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import javax.imageio.ImageIO;
import javax.media.jai.JAI;
import javax.media.jai.KernelJAI;
import javax.media.jai.RenderedOp;


public class ModifiedCensusStereoSobel {
    
    public static void main(String[] args) {
         
        //define images TESTING
        RenderedImage NadirImage = null;
        RenderedImage ForwardImage = null;
        
        //read in images TESTING
        try {
            NadirImage = ImageIO.read(new File ("test_nadir.png"));
            } catch(Exception e) {}
        
        try {
            ForwardImage = ImageIO.read(new File ("test_forward.png"));
            } catch(Exception e) {}

        //convert to arrays TESTING
        int imageWidth = NadirImage.getWidth();
        int imageHeight = NadirImage.getHeight(); 
        Raster NadirRaster = NadirImage.getData();
        Raster ForwardRaster = ForwardImage.getData();
        
        int[][] nadirArray = new int[imageWidth][imageHeight];
        int[][] forwardArray = new int[imageWidth][imageHeight];
        
        for (int i = 0; i < imageWidth; i++) {
            for (int j = 0; j < imageHeight; j++) {
                
                nadirArray[i][j] = NadirRaster.getSample(i, j, 0);
                forwardArray[i][j] = ForwardRaster.getSample(i, j, 0);
            }
        }
        
        //create local variables
        int[][] PreCostArray = new int[imageWidth][imageHeight];
        int[][] CostArray = new int[imageWidth][imageHeight];
        int[][] xDisArray = new int[imageWidth][imageHeight];
        int[][] yDisArray = new int[imageWidth][imageHeight];
        int[][][] NadirCensusMap; 
        int[][][] ForwardCensusMap;
        int nBits = 8;
        int windowWidth = 9;
        int windowHeight = 9;
        int nBitNumbers = 2 * (windowWidth * windowHeight / 8) + 1;
        int halfWindowWidth = windowWidth/2;
        int halfWindowHeight = windowHeight/2;
        int imageWidthMinusHalfWindowWidth = imageWidth - halfWindowWidth;
        int imageHeightMinusHalfWindowHeight = imageHeight - halfWindowHeight;
        int xShift = 5;
        int yShift = 12;
        int NadirBitNumber;
        int ForwardBitNumber;
        int smoothingRadius = 1;
        HashMap<Integer, Integer> LookUpTable;
        int[] inPixels = new int[imageWidth*imageHeight];
        int[] outPixels = new int[imageWidth*imageHeight];
        int[] smoothedPixels = new int[imageHeight*imageWidth];
        int[] nadirInPixels = new int[imageWidth*imageHeight];
        int[] forwardInPixels = new int[imageWidth*imageHeight];
        int[] sobelNadirInPixels = new int[imageWidth*imageHeight];
        int[] sobelForwardInPixels = new int[imageWidth*imageHeight];
        int[] smoothedNadirPixels = new int[imageHeight*imageWidth];
        int[] smoothedForwardPixels = new int[imageHeight*imageWidth];
        int[] smoothedSobelNadirPixels = new int[imageHeight*imageWidth];
        int[] smoothedSobelForwardPixels = new int[imageHeight*imageWidth];
        ParameterBlock nadirPB = new ParameterBlock(); 
        ParameterBlock forwardPB = new ParameterBlock(); 
        int[][] nadirSobelArray = new int[imageWidth][imageHeight];
        int[][] forwardSobelArray = new int[imageWidth][imageHeight];
        int[][] smoothedNadirArray = new int[imageWidth][imageHeight];
        int[][] smoothedForwardArray = new int[imageWidth][imageHeight];
        int[][] smoothedNadirSobelArray = new int[imageWidth][imageHeight];
        int[][] smoothedForwardSobelArray = new int[imageWidth][imageHeight];
        KernelJAI sobelH = KernelJAI.GRADIENT_MASK_SOBEL_HORIZONTAL;
        KernelJAI sobelV = KernelJAI.GRADIENT_MASK_SOBEL_VERTICAL;
        ParameterBlock nadirSmoothPB = new ParameterBlock();
        ParameterBlock forwardSmoothPB = new ParameterBlock();
        
        //sobel filter (pre smoothing to get rid of noise effects) WILL NEED TO BE MODIFIED FOR THE MAIN ROUTINE TO CREATE THE RENDERED IMAGE INPUTS FROM THE INPUT ARRAYS
        nadirSmoothPB.addSource(NadirImage);
        nadirSmoothPB.add(3);
        nadirSmoothPB.add(3);
        nadirSmoothPB.add(3/2);
        nadirSmoothPB.add(3/2);
        RenderedOp smoothedNadirImage = JAI.create("boxfilter", nadirSmoothPB);
        
        nadirPB.addSource(smoothedNadirImage);
        nadirPB.add(sobelH);
        nadirPB.add(sobelV);
        RenderedOp NadirSobel = JAI.create( "gradientmagnitude", nadirPB);
        
        forwardSmoothPB.addSource(ForwardImage);
        forwardSmoothPB.add(3);
        forwardSmoothPB.add(3);
        forwardSmoothPB.add(3/2);
        forwardSmoothPB.add(3/2);
        RenderedOp smoothedForwardImage = JAI.create("boxfilter", forwardSmoothPB);
        
        forwardPB.addSource(smoothedForwardImage);
        forwardPB.add(sobelH);
        forwardPB.add(sobelV);
        RenderedOp ForwardSobel = JAI.create( "gradientmagnitude", forwardPB);
        
        Raster NadirSobelRaster = NadirSobel.getData();
        Raster ForwardSobelRaster = ForwardSobel.getData();
        
        //save an image of the gradienti mages for testing purposes   
        File nadirSobelFile = new File("nadirSobel.png");
        try {
            ImageIO.write(NadirSobel, "png", nadirSobelFile);
        } catch (IOException e) {
            e.printStackTrace();
        }
        
        File forwardSobelFile = new File("forwardSobel.png");
        try {
            ImageIO.write(ForwardSobel, "png", forwardSobelFile);
        } catch (IOException e) {
            e.printStackTrace();
        }
        
        
        for (int i = 0; i < imageWidth; i++) {
            for (int j = 0; j < imageHeight; j++) {
                
                nadirSobelArray[i][j] = NadirSobelRaster.getSample(i, j, 0);
                forwardSobelArray[i][j] = ForwardSobelRaster.getSample(i, j, 0);
            }
        }
                
        //mean filter each array to get central pixel
        int count = 0;
        for (int x = 0; x < imageWidth; x++) {
            for (int y = 0; y < imageHeight; y++) {
                
                nadirInPixels[count] = nadirArray[x][y];
                forwardInPixels[count] = forwardArray[x][y];
                sobelNadirInPixels[count] = nadirSobelArray[x][y];
                sobelForwardInPixels[count] = forwardSobelArray[x][y];
                count++;
            }
        }
        
        //box filter average
        blur(nadirInPixels, outPixels, imageWidth, imageHeight, windowWidth/2);
        blur(outPixels, smoothedNadirPixels, imageHeight, imageWidth, windowHeight/2);
        blur(forwardInPixels, outPixels, imageWidth, imageHeight, windowWidth/2);
        blur(outPixels, smoothedForwardPixels, imageHeight, imageWidth, windowHeight/2);
        blur(sobelNadirInPixels, outPixels, imageWidth, imageHeight, windowWidth/2);
        blur(outPixels, smoothedSobelNadirPixels, imageHeight, imageWidth, windowHeight/2);
        blur(sobelForwardInPixels, outPixels, imageWidth, imageHeight, windowWidth/2);
        blur(outPixels, smoothedSobelForwardPixels, imageHeight, imageWidth, windowHeight/2);
        
        //back to 2d
        count = 0;
        for (int x = 0; x < imageWidth; x++) {
            for (int y = 0; y < imageHeight; y++) {
                
                smoothedNadirArray[x][y] = smoothedNadirPixels[count];
                smoothedForwardArray[x][y] = smoothedForwardPixels[count];
                smoothedNadirSobelArray[x][y] = smoothedSobelNadirPixels[count];
                smoothedForwardSobelArray[x][y] = smoothedSobelForwardPixels[count];
                count++;
            }
        }
        
        //create census arrays
        NadirCensusMap = CensusImage(nadirArray, nadirSobelArray, smoothedNadirArray, smoothedNadirSobelArray, imageWidth, imageHeight, windowWidth, windowHeight); 
        ForwardCensusMap = CensusImage(forwardArray, forwardSobelArray, smoothedForwardArray, smoothedForwardSobelArray, imageWidth, imageHeight, windowWidth, windowHeight);
        
        //make 8bit cost look up
        LookUpTable = makeLookUpTable(nBits);
        
        //initialise cost rasters
        for (int x = halfWindowWidth; x < imageWidthMinusHalfWindowWidth; x++) {
            for (int y = halfWindowHeight; y < imageHeightMinusHalfWindowHeight; y++) {
 
                PreCostArray[x][y] = 999;
                CostArray[x][y] = 999;
            }           
        }
                 
        //loop over disparities
        for (int dx = -xShift; dx < xShift ; dx++) {
            for (int dy = 0; dy < yShift; dy++) {
                
                //loop over images
                for (int x = halfWindowWidth; x < imageWidthMinusHalfWindowWidth; x++) {
                    //System.out.println("x: " + x);
                    for (int y = halfWindowHeight; y < imageHeightMinusHalfWindowHeight; y++) {
                        
                        //create xdx and ydy
                        int xdx = x+dx;
                        int ydy = y+dy;
                        
                        //check if forward location in range
                        if (xdx >= halfWindowWidth && xdx < imageWidthMinusHalfWindowWidth
                           && ydy >= halfWindowHeight && ydy < imageHeightMinusHalfWindowHeight ) {
                            
                            //compute hamming distance
                            int hammingDistance = 0;
                            for (int d = 0; d < nBitNumbers; d++) {
                                
                                NadirBitNumber = NadirCensusMap[x][y][d]; 
                                ForwardBitNumber = ForwardCensusMap[xdx][ydy][d];
                                
                                //xor the bits
                                int xorBitNumber = NadirBitNumber ^ ForwardBitNumber;
                                       
                                //find associated cost in LUT
                                hammingDistance+=LookUpTable.get(xorBitNumber);
                                
                            }
                                  
                            //update precost array
                            PreCostArray[x][y] = hammingDistance;
                        
                        } //in range    
                    } //for y
                } //for x
                
                //convert pre cost array to 1d
                count = 0;
                for (int x = 0; x < imageWidth; x++) {
                    for (int y = 0; y < imageHeight; y++) {
                        
                        inPixels[count] = PreCostArray[x][y];
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
                       
                       int preCost = smoothedPixels[count];
                       count++;
                        //check and update
                        if (preCost < CostArray[x][y]) {
                            CostArray[x][y] = PreCostArray[x][y];
                            xDisArray[x][y] = dx;
                            yDisArray[x][y] = dy;
                        } 
                    } 
                }      
            } // for dy
        } //for dx
        
        //median filter output array
        //int[][] smoothedYDisArray = new int[imageWidth][imageHeight];
        //int minWidth = halfWindowWidth + smoothingRadius;
        //int maxWidth = imageWidth - minWidth;
        //int minHeight = halfWindowHeight + smoothingRadius;
        //int maxHeight = imageHeight - minHeight;
            
        //MedianFilter(yDisArray, smoothedYDisArray, minWidth, maxWidth, minHeight, maxHeight, smoothingRadius*2+1, 2);
        
        //Convert to writabelrasters TESTING
        WritableRaster yDisOutput = Raster.createWritableRaster(NadirImage.getSampleModel(), null); 
       
        for (int x = 0; x < imageHeight; x++) {
            for (int y = 0; y < imageWidth; y++) {
                yDisOutput.setSample(x, y, 0, yDisArray[x][y]);
            
            }
        
        }
       
        //write out to BufferedImage TESTING
        BufferedImage yDisBI = new BufferedImage(NadirImage.getColorModel(), yDisOutput, false, null);
       
        File yDisFile = new File("yDis.png");
        try {
            ImageIO.write(yDisBI, "png", yDisFile);
        } catch (IOException e) {
            e.printStackTrace();
        }
        
        
    }  //end of main
    
    
    public static void MedianFilter(int[][] inputImage, int[][] outputImage, int minWidth, int maxWidth, int minHeight, int maxHeight, int blockSize, double covThresh){
        
        int halfWindow = (blockSize)/2;
        int[] windowValues=new int[blockSize*blockSize];
        Quicksort Sorter = new Quicksort();
        
        //loop over image
        for (int x = minWidth; x < maxWidth; x++) {   
            for(int y = minHeight; y < maxHeight; y++) {   
                
                //reset local count
                int windowCount=0;
                int windowSum = 0;
                
                //loop over filter window and populate median array
                for (int i = -halfWindow; i <= halfWindow; i++){              
                    int ii=x+i;          
                    if (ii >= minWidth && ii < maxWidth){            
                        for (int j = -halfWindow;j <= halfWindow; j++){                                                  
                            int jj=y+j;                  
                            if (jj >= minHeight && jj < maxHeight){                      
                                if (i==0 || j==0 || i==j || i==-j){                                                 
                                    int value = inputImage[ii][jj];                          
                                    windowValues[windowCount] = value;
                                    windowSum += value; 
                                    windowCount++;                       
                                }                   
                            }               
                        }           
                    }       
                }
                
                //trim window values to the size of the window and calculate stdev
                int[] inboundWindowValues = new int[windowCount];
                double windowMean = windowSum/windowCount;
                double dispersion = 0;
                double std = 0;
                double cov = 0;
                for (int n = 0; n < windowCount; n++) {
                    dispersion += (windowValues[n] - windowMean) * (windowValues[n] - windowMean);
                    inboundWindowValues[n] = windowValues[n];                       
                }
                std = Math.sqrt(dispersion);
                
                //if coefficient of variation detects outlier then replace pixel with median
                //only works on ratio scale (positive numbers) so test using the absolute of
                //of the mean.  
                Math.abs(windowMean);
                cov = std / windowMean;
                
              
                
                if (cov > covThresh) {
                    //quicksort inbound values
                    Sorter.sort(inboundWindowValues);
                    
                    //determine median
                    int middle = windowCount/2-1;
                    if ((windowCount & 1) == 1 ) {
                        outputImage[x][y] = inboundWindowValues[middle]; 
                    } else {
                        outputImage[x][y] = (int) ((inboundWindowValues[middle] + inboundWindowValues[middle+1] + 0.5) / 2);
                    }
                    
                } else {
                    
                    outputImage[x][y] = inputImage[x][y];

                }
                
                       
            }
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
    
    
    public static HashMap<Integer, Integer> makeLookUpTable(int nBits) {
            
        //local variables
        int nCombinations =  (int) Math.pow(2, nBits);
        int value;
        HashMap<Integer, Integer> LookUpTable = new HashMap<Integer, Integer>(nCombinations);
            
        //create the binaryArray
        int[][] binaryArray = MakeBinaryArray(nBits);
            
        //loop over the 8 bit set
        for (int i = 0; i < nCombinations; i++) {
                           
            //loop over the number of bits
            value = 0;
            for (int n = 0; n < nBits; n++) {
                
                //get hamnming distance
                if (binaryArray[i][n] == 1) {
                    value+=1;
                }
                        
                }
                    
            //store bit Number and score 
            LookUpTable.put(i, value);
                       
        }
        return LookUpTable;
    } //end of MakeLookUpTable
        
    
    //routine for creating the binary array
    public static int[][] MakeBinaryArray(int N) {
            
        int nLength = N-1;
        int binaries = (int) Math.pow(2, N);
        int[][] codeArray = new int[binaries][N];
            
        //loop over all numbers start at 1 as binary 0 is all zeros
        for (int i = 1; i < binaries; i++) {
                
            //convert number to binary
            String binaryString = Integer.toBinaryString(i);
            int stringLength = binaryString.length();
         
            //loop over size of binary number and place into array at specified j        
            int counter = 1;
            for (int j = nLength; j > nLength - stringLength; j--) {
                    
                //extract part of the string
                Character bitChar = new Character(binaryString.charAt(stringLength - counter));
                String bitString = bitChar.toString();
                    
                //parse to integer
                int bitInt = Integer.parseInt(bitString);
                    
                //place in array
                codeArray[i][j] = bitInt;
                    
                //update counter
                counter++;    
            }
                    
        }
        return codeArray;
                
    }  //end of MakeBinaryArray
    
    public static int[][][] CensusImage(int[][] intensityImage, int[][] gradientImage, int[][] meanIntensityImage, int[][] meanGradientImage, int width, int height, int winWidth, int winHeight) {
        
        //variables
        int length = winWidth*winHeight;
        int pixel;
        int testPixel;
        int winHalfWidth = winWidth/2;
        int winHalfHeight = winHeight/2;
        StringBuilder bitString = new StringBuilder(length);
        String inputBitString;
        
        //depth is the number of binary number groups 
        int depth = ((winWidth * winHeight) / 8) + 1;
        int[][][] census = new int[width][height][depth*2];
   
        //loop over intensity image
        for (int x = winHalfWidth; x < width-winHalfWidth; x++) {
            for (int y = winHalfHeight; y < height-winHalfHeight; y++) {
                
                //get mean intensity value for comparison with intensity search window 
                pixel = meanIntensityImage[x][y];
                
                //loop over window and get bits 
                for (int i = x-winHalfWidth; i <= x+winHalfWidth ; i++) {
                    for (int j = y-winHalfHeight; j <= y+winHalfHeight; j++) {
                        testPixel = intensityImage[i][j];
                        
                        //check census
                        if (pixel < testPixel) {
                            
                            bitString.append(0);
                            
                        } else {
                            
                            bitString.append(1);
                        } 
                    }
                }
                
                //put value into array as a bit number
                inputBitString = bitString.toString();
                
                int count = 0;
                for (int d = 0; d <= inputBitString.length(); d+=8) {
                    int dend = d+8;
                    if (dend > length) {
                        dend = length;
                    }
                    String inputBitSubString = inputBitString.substring(d, dend);
                    census[x][y][count] = Integer.parseInt(inputBitSubString, 2);
                    count++;
                }
                
                
                //clear bit string
                bitString.delete(0,length);
            }
        }
        
        //loop over gradient image (easiest to do it with two loops)
        for (int x = winHalfWidth; x < width-winHalfWidth; x++) {
            for (int y = winHalfHeight; y < height-winHalfHeight; y++) {
                
                //get mean gradient value for comparison with gradient search window
                pixel = meanGradientImage[x][y];
                
                //loop over window and get bits 
                for (int i = x-winHalfWidth; i <= x+winHalfWidth ; i++) {
                    for (int j = y-winHalfHeight; j <= y+winHalfHeight; j++) {
                        testPixel = gradientImage[i][j];
                        
                        //check census
                        if (pixel < testPixel) {
                            
                            bitString.append(0);
                            
                        } else {
                            
                            bitString.append(1);
                        } 
                    }
                }
                
                //put value into array as a bit number
                inputBitString = bitString.toString();
                
                int count = depth;
                for (int d = 0; d <= inputBitString.length(); d+=8) {
                    int dend = d+8;
                    if (dend > length) {
                        dend = length;
                    }
                    String inputBitSubString = inputBitString.substring(d, dend);
                    census[x][y][count] = Integer.parseInt(inputBitSubString, 2);
                    count++;
                }
                
                
                //clear bit string
                bitString.delete(0,length);
            }
        }
        
        return census;
    }  //end of MakeCensusImage

    
    //edge clamper
    public static float clamp(float x, float a, float b) {
        return (x < a) ? a : (x > b) ? b : x;
    }
    
    


    
}
  