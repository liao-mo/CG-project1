///////////////////////////////////////////////////////////////////////////////
//
//      TargaImage.cpp                          Author:     Stephen Chenney
//                                              Modified:   Eric McDaniel
//                                              Date:       Fall 2004
//
//      Implementation of TargaImage methods.  You must implement the image
//  modification functions.
//
///////////////////////////////////////////////////////////////////////////////

#include "Globals.h"
#include "TargaImage.h"
#include "libtarga.h"
#include <stdlib.h>
#include <assert.h>
#include <memory.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include <queue>
#include <algorithm>

using namespace std;

// constants
const int           RED             = 0;                // red channel
const int           GREEN           = 1;                // green channel
const int           BLUE            = 2;                // blue channel
const unsigned char BACKGROUND[3]   = { 0, 0, 0 };      // background color


// Computes n choose s, efficiently
double Binomial(int n, int s)
{
    double        res;

    res = 1;
    for (int i = 1 ; i <= s ; i++)
        res = (n - i + 1) * res / i ;

    return res;
}// Binomial


///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage() : width(0), height(0), data(NULL), original_data(NULL), gray_data(NULL)
{}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h) : width(w), height(h)
{
   data = new unsigned char[width * height * 4];
   original_data = new unsigned char[width * height * 4];
   gray_data = new double[width * height];
   ClearToBlack();
}// TargaImage



///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables to values given.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h, unsigned char* d)
{
    int i;

    width = w;
    height = h;
    data = new unsigned char[width * height * 4];
    original_data = new unsigned char[width * height * 4];
    gray_data = new double[width * height];

    for (i = 0; i < width * height * 4; i++){
        data[i] = d[i];
        original_data[i] = d[i];
    }
}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Copy Constructor.  Initialize member to that of input
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(const TargaImage& image) 
{
   width = image.width;
   height = image.height;
   data = NULL; 
   original_data = NULL;
   gray_data = NULL;
   if (image.data != NULL) {
      data = new unsigned char[width * height * 4];
      original_data = new unsigned char[width * height * 4];
      gray_data = new double[width * height];
      memcpy(data, image.data, sizeof(unsigned char) * width * height * 4);
      memcpy(original_data, image.data, sizeof(unsigned char) * width * height * 4);
   }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Destructor.  Free image memory.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::~TargaImage()
{
    if (data)
        delete[] data;
}// ~TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Converts an image to RGB form, and returns the rgb pixel data - 24 
//  bits per pixel. The returned space should be deleted when no longer 
//  required.
//
///////////////////////////////////////////////////////////////////////////////
unsigned char* TargaImage::To_RGB(void)
{
    unsigned char* rgb = new unsigned char[width * height * 3];
    int i, j;

    if (! data)
	    return NULL;

    // Divide out the alpha
    for (i = 0 ; i < height ; i++)
    {
	    int in_offset = i * width * 4;
	    int out_offset = i * width * 3;

	    for (j = 0 ; j < width ; j++)
        {
	        RGBA_To_RGB(data + (in_offset + j*4), rgb + (out_offset + j*3));
	    }
    }

    return rgb;
}// TargaImage

//Get float array by the current data
float* TargaImage::get_float_data(void) {
    float* f_data = new float[width * height * 4];
    for (int i = 0; i < width * height * 4; ++i) {
        f_data[i] = (float)data[i];
    }
    return f_data;
}

///////////////////////////////////////////////////////////////////////////////
//
//      Save the image to a targa file. Returns 1 on success, 0 on failure.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Save_Image(const char *filename)
{
    TargaImage	*out_image = Reverse_Rows();

    if (! out_image)
	    return false;

    if (!tga_write_raw(filename, width, height, out_image->data, TGA_TRUECOLOR_32))
    {
	    cout << "TGA Save Error: %s\n", tga_error_string(tga_get_last_error());
	    return false;
    }

    delete out_image;

    return true;
}// Save_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Load a targa image from a file.  Return a new TargaImage object which 
//  must be deleted by caller.  Return NULL on failure.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Load_Image(char *filename)
{
    unsigned char   *temp_data;
    TargaImage	    *temp_image;
    TargaImage	    *result;
    int		        width, height;

    if (!filename)
    {
        cout << "No filename given." << endl;
        return NULL;
    }// if

    temp_data = (unsigned char*)tga_load(filename, &width, &height, TGA_TRUECOLOR_32);
    if (!temp_data)
    {
        cout << "TGA Error: %s\n", tga_error_string(tga_get_last_error());
	    width = height = 0;
	    return NULL;
    }
    temp_image = new TargaImage(width, height, temp_data);
    free(temp_data);

    result = temp_image->Reverse_Rows();

    delete temp_image;

    return result;
}// Load_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Convert image to grayscale.  Red, green, and blue channels should all 
//  contain grayscale value.  Alpha channel shoould be left unchanged.  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::To_Grayscale()
{
    if (!data)
        return false;

    for (int i = 0; i < height; i++)
    {
        int in_offset = i * width * 4;

        for (int j = 0; j < width; j++)
        {
            //iterate over the array, calculate the gray value and assign to each channel
            unsigned char gray = *(data + in_offset + j*4) * 0.299 + *(data + in_offset + j*4 + 1) * 0.587 + *(data + in_offset + j*4 + 2) * 0.114;
            gray_data[i * width + j] = gray;
            data[in_offset + j * 4 + 0] = gray;
            data[in_offset + j * 4 + 1] = gray;
            data[in_offset + j * 4 + 2] = gray;
        }
    }
    return true;
}// To_Grayscale


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using uniform quantization.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Uniform()
{
    if (!data)
        return false;

    for (int i = 0; i < height; i++)
    {
        int in_offset = i * width * 4;

        for (int j = 0; j < width; j++)
        {
            //find quantized values of each channel, and assign them back
            unsigned char R_val = (data[in_offset + j * 4 + 0] / 32) * 32;
            unsigned char G_val = (data[in_offset + j * 4 + 1] / 32) * 32;
            unsigned char B_val = (data[in_offset + j * 4 + 2] / 64) * 64;
            
            data[in_offset + j * 4 + 0] = R_val;
            data[in_offset + j * 4 + 1] = G_val;
            data[in_offset + j * 4 + 2] = B_val;
        }
    }
    

    return true;
}// Quant_Uniform


///////////////////////////////////////////////////////////////////////////////
//
//      Convert the image to an 8 bit image using populosity quantization.  
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Populosity()
{
    if (!data)
        return false;

    unsigned char* rgb_data = To_RGB();
    unsigned temp_rgb[3];
    
    vector<color_data> pop_colors;
    unordered_map<color_data, int, color_data_hasher> u_map;
    
    //put all the colors into the unordered map
    for (int i = 0; i < width * height; ++i) {
        color_data temp_color;
        int offset = i * 3;
        //use bit operation to do the first step of quantization
        temp_color.r_val = (rgb_data[offset + 0] >> 3) << 3;
        temp_color.g_val = (rgb_data[offset + 1] >> 3) << 3;
        temp_color.b_val = (rgb_data[offset + 2] >> 3) << 3;
        u_map[temp_color]++;
    }
    //find the most popular 256 colors from unorderd map by using the priority queue
    priority_queue<pair<int, color_data>> pq;
    for (auto it = u_map.begin(); it != u_map.end(); ++it) {
        pq.push(make_pair(it->second, it->first));
        if (pq.size() > u_map.size() - 256) {
            pop_colors.push_back(pq.top().second);
            pq.pop();
        }
    }

    //map each pixel to the 256 color table
    for (int i = 0; i < width * height; ++i) {
        int in_offset = i * 3;
        int out_offset = i * 4;
        double min_distance_square = DBL_MAX;
        color_data new_color;
        //find the closest color of this pixel
        for (color_data c : pop_colors) {
            double dR = rgb_data[in_offset + 0] - c.r_val;
            double dG = rgb_data[in_offset + 1] - c.g_val;
            double dB = rgb_data[in_offset + 2] - c.b_val;
            double distance_square = dR * dR + dG * dG + dB * dB;
            if (distance_square < min_distance_square) {
                min_distance_square = distance_square;
                new_color = c;
            }
        }
        //assign back to the current data
        data[out_offset + 0] = new_color.r_val;
        data[out_offset + 1] = new_color.g_val;
        data[out_offset + 2] = new_color.b_val;
    }
    return true;
}// Quant_Populosity


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image using a threshold of 1/2.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Threshold()
{
    if (!data)
        return false;

    To_Grayscale();
    for (int i = 0; i < width * height; ++i) {
        int offset = i * 4;
        if ((gray_data[i] / 255) < 0.5) {
            gray_data[i] = 0;
        }
        else {
            gray_data[i] = 255;
        }
        data[offset + 0] = gray_data[i];
        data[offset + 1] = gray_data[i];
        data[offset + 2] = gray_data[i];
    }

    return true;
}// Dither_Threshold


///////////////////////////////////////////////////////////////////////////////
//
//      Dither image using random dithering.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Random()
{
    if (!data)
        return false;
    To_Grayscale();
    for (int i = 0; i < width * height; ++i) {
        int offset = i * 4;
        double randNum = ((0.4 * ((double)rand() / (double)RAND_MAX)) - 0.2);
        double temp = (gray_data[i] / 255) + randNum;
        if (temp < 0.5) {
            gray_data[i] = 0;
        }
        else {
            gray_data[i] = 255;
        }
        data[offset + 0] = gray_data[i];
        data[offset + 1] = gray_data[i];
        data[offset + 2] = gray_data[i];
    }



    return true;
}// Dither_Random


///////////////////////////////////////////////////////////////////////////////
//
//      Perform Floyd-Steinberg dithering on the image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_FS()
{
    if (!data)
        return false;

    To_Grayscale();
    for (int i = 0; i < height; ++i) {
        //odd row or even row has different order
        if (i % 2 == 0) {
            //from right to left
            for (int j = width - 1; j >= 0; --j) {
                int current_index = i * width + j;
                int next_row_index = (i + 1) * width + j;

                //calculate the error with the threshold 0.5
                double error = 0;
                if (gray_data[current_index] / 255.0 < 0.5) {
                    error = gray_data[current_index];
                    gray_data[current_index] = 0;
                }
                else {
                    error = gray_data[current_index] - 255;
                    gray_data[current_index] = 255;
                }

                //distribute the error to the neighboring pixels
                //start
                if (j == width - 1) {
                    gray_data[current_index - 1] += error * (7.0 / 16.0);
                    if (i != height - 1) {
                        gray_data[next_row_index + 0] += error * (5.0 / 16.0);
                        gray_data[next_row_index - 1] += error * (1.0 / 16.0);
                    }
                }
                //end
                else if (j == 0) {
                    if (i != height - 1) {
                        gray_data[next_row_index + 0] += error * (5.0 / 16.0);
                        gray_data[next_row_index + 1] += error * (3.0 / 16.0);
                    }
                }
                //between
                else {
                    gray_data[current_index - 1] += error * (7.0 / 16.0);
                    if (i != height - 1) {
                        gray_data[next_row_index + 1] += error * (3.0 / 16.0);
                        gray_data[next_row_index + 0] += error * (5.0 / 16.0);
                        gray_data[next_row_index - 1] += error * (1.0 / 16.0);
                    }
                }
            }
        }
        else {
            //from left to right
            for (int j = 0; j < width; ++j) {
                int current_index = i * width + j;
                int next_row_index = (i + 1) * width + j;

                //calculate the error with the threshold 0.5
                double error = 0;
                if (gray_data[current_index] / 255.0 < 0.5) {
                    error = gray_data[current_index];
                    gray_data[current_index] = 0;
                }
                else {
                    error = gray_data[current_index] - 255;
                    gray_data[current_index] = 255;
                }

                //distribute the error to the neighboring pixels
                //start
                if (j == 0) {
                    gray_data[current_index + 1] += error * (7.0 / 16.0);
                    if (i != height - 1) {
                        gray_data[next_row_index + 0] += error * (5.0 / 16.0);
                        gray_data[next_row_index + 1] += error * (1.0 / 16.0);
                    }
                }
                //end
                else if (j == width - 1) {
                    if (i != height - 1) {
                        gray_data[next_row_index + 0] += error * (5.0 / 16.0);
                        gray_data[next_row_index - 1] += error * (3.0 / 16.0);
                    }
                }
                //between
                else {
                    gray_data[current_index + 1] += error * (7.0 / 16.0);
                    if (i != height - 1) {
                        gray_data[next_row_index - 1] += error * (3.0 / 16.0);
                        gray_data[next_row_index + 0] += error * (5.0 / 16.0);
                        gray_data[next_row_index + 1] += error * (1.0 / 16.0);
                    }
                }
            }
        }
    }
    
    for (int i = 0; i < width * height; ++i) {
        int offset = i * 4;
        data[offset + 0] = gray_data[i];
        data[offset + 1] = gray_data[i];
        data[offset + 2] = gray_data[i];
    }





    return true;
}// Dither_FS


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image while conserving the average brightness.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Bright()
{
    if (!data)
        return false;

    To_Grayscale();
    vector<int> intensities(255, 0);
    double sum = 0;
    //keep all the intensities into a vector
    for (int i = 0; i < width * height; ++i) {
        intensities[(int)gray_data[i]]++;
        sum += gray_data[i];
    }

    //find the average of intensity
    double average = sum / (width * height);
    //find the number of dark pixels
    int num_of_dark = (width * height) * (1 - (average / 255));
    //find the threshold by counting the number of dark pixels
    int intensity = 0;
    for (intensity = 0; intensity < 256; ++intensity) {
        num_of_dark -= intensities[intensity];
        if (num_of_dark <= 0) {
            break;
        }
    }
    double threshold = (double)intensity / 255.0;

    for (int i = 0; i < width * height; ++i) {
        int offset = i * 4;
        if ((gray_data[i] / 255) < threshold) {
            gray_data[i] = 0;
        }
        else {
            gray_data[i] = 255;
        }
        data[offset + 0] = gray_data[i];
        data[offset + 1] = gray_data[i];
        data[offset + 2] = gray_data[i];

    }
    return true;
}// Dither_Bright


///////////////////////////////////////////////////////////////////////////////
//
//      Perform clustered differing of the image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Cluster()
{
    if (!data)
        return false;

    To_Grayscale();
    vector<vector<double>> mask = { 
        {0.7059, 0.3529, 0.5882, 0.2353},
        {0.0588, 0.9412, 0.8235, 0.4118},
        {0.4706, 0.7647, 0.8824, 0.1176},
        {0.1765, 0.5294, 0.2941, 0.6471}
    };

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            int in_offset = y * width;
            int out_offset = y * width * 4;
            if ((gray_data[in_offset + x] / 255.0) < mask[x % 4][y % 4]) {
                gray_data[in_offset + x] = 0;
            }
            else {
                gray_data[in_offset + x] = 255;
            }
            data[out_offset + x * 4 + 0] = gray_data[in_offset + x];
            data[out_offset + x * 4 + 1] = gray_data[in_offset + x];
            data[out_offset + x * 4 + 2] = gray_data[in_offset + x];
        }
    }



    return true;
}// Dither_Cluster


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using Floyd-Steinberg dithering over
//  a uniform quantization - the same quantization as in Quant_Uniform.
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Color()
{
    if (!data)
        return false;

    //create vectors of each channel
    vector<double> r_data(width * height, 0);
    vector<double> g_data(width * height, 0);
    vector<double> b_data(width * height, 0);
    vector<vector<double>> rgb_channels = { r_data, g_data, b_data };
    vector<double> r_error(width * height, 0);
    vector<double> g_error(width * height, 0);
    vector<double> b_error(width * height, 0);
    vector<vector<double>> rgb_errors = { r_error, g_error, b_error };

    //do the quant-unif and calculate the error during this step
    for (int i = 0; i < height; i++)
    {
        int offset = i * width * 4;

        for (int j = 0; j < width; j++)
        {
            //find quantized values of each channel, and assign them back
            unsigned char R_val = (data[offset + j * 4 + 0] / 32) * 32;
            unsigned char G_val = (data[offset + j * 4 + 1] / 32) * 32;
            unsigned char B_val = (data[offset + j * 4 + 2] / 64) * 64;

            rgb_channels[0][i * width + j] = R_val;
            rgb_channels[1][i * width + j] = G_val;
            rgb_channels[2][i * width + j] = B_val;
            rgb_errors[0][i * width + j] = data[offset + j * 4 + 0] - R_val;
            rgb_errors[1][i * width + j] = data[offset + j * 4 + 1] - G_val;
            rgb_errors[2][i * width + j] = data[offset + j * 4 + 2] - B_val;
        }
    }


    //dither each channel
    for (int channel = 0; channel < 3; ++channel) {
        for (int i = 0; i < height; ++i) {
            //odd row or even row has different order
            if (i % 2 == 0) {
                //from right to left
                for (int j = width - 1; j >= 0; --j) {
                    int current_index = i * width + j;
                    int next_row_index = (i + 1) * width + j;

                    //take the error from previous step
                    double error = rgb_errors[channel][current_index];
                    //cout << error << endl;
                    //distribute the error to the neighboring pixels
                    //start
                    if (j == width - 1) {
                        rgb_channels[channel][current_index - 1] += error * (7.0 / 16.0);
                        if (i != height - 1) {
                            rgb_channels[channel][next_row_index + 0] += error * (5.0 / 16.0);
                            rgb_channels[channel][next_row_index - 1] += error * (1.0 / 16.0);
                        }
                    }
                    //end
                    else if (j == 0) {
                        if (i != height - 1) {
                            rgb_channels[channel][next_row_index + 0] += error * (5.0 / 16.0);
                            rgb_channels[channel][next_row_index + 1] += error * (3.0 / 16.0);
                        }
                    }
                    //between
                    else {
                        rgb_channels[channel][current_index - 1] += error * (7.0 / 16.0);
                        if (i != height - 1) {
                            rgb_channels[channel][next_row_index + 1] += error * (3.0 / 16.0);
                            rgb_channels[channel][next_row_index + 0] += error * (5.0 / 16.0);
                            rgb_channels[channel][next_row_index - 1] += error * (1.0 / 16.0);
                        }
                    }
                }
            }
            else {
                //from left to right
                for (int j = 0; j < width; ++j) {
                    int current_index = i * width + j;
                    int next_row_index = (i + 1) * width + j;

                    //take the error from previous step
                    double error = rgb_errors[channel][current_index];
                    //cout << error  << endl;

                    //distribute the error to the neighboring pixels
                    //start
                    if (j == 0) {
                        rgb_channels[channel][current_index + 1] += error * (7.0 / 16.0);
                        if (i != height - 1) {
                            rgb_channels[channel][next_row_index + 0] += error * (5.0 / 16.0);
                            rgb_channels[channel][next_row_index + 1] += error * (1.0 / 16.0);
                        }
                    }
                    //end
                    else if (j == width - 1) {
                        if (i != height - 1) {
                            rgb_channels[channel][next_row_index + 0] += error * (5.0 / 16.0);
                            rgb_channels[channel][next_row_index - 1] += error * (3.0 / 16.0);
                        }
                    }
                    //between
                    else {
                        gray_data[current_index + 1] += error * (7.0 / 16.0);
                        if (i != height - 1) {
                            rgb_channels[channel][next_row_index - 1] += error * (3.0 / 16.0);
                            rgb_channels[channel][next_row_index + 0] += error * (5.0 / 16.0);
                            rgb_channels[channel][next_row_index + 1] += error * (1.0 / 16.0);
                        }
                    }
                }
            }
        }
    }
    

    for (int i = 0; i < width * height; ++i) {
        int offset = i * 4;
        data[offset + 0] = (int)rgb_channels[0][i];
        data[offset + 1] = (int)rgb_channels[1][i];
        data[offset + 2] = (int)rgb_channels[2][i];

        
    }





    return true;
}// Dither_Color


///////////////////////////////////////////////////////////////////////////////
//
//      Composite the current image over the given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Over(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout <<  "Comp_Over: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_Over


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "in" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_In(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_In: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_In


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "out" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Out(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Out: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_Out


///////////////////////////////////////////////////////////////////////////////
//
//      Composite current image "atop" given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Atop(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Atop: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_Atop


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image with given image using exclusive or (XOR).  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Xor(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Xor: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_Xor


///////////////////////////////////////////////////////////////////////////////
//
//      Calculate the difference bewteen this imag and the given one.  Image 
//  dimensions must be equal.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Difference(TargaImage* pImage)
{
    if (!pImage)
        return false;

    if (width != pImage->width || height != pImage->height)
    {
        cout << "Difference: Images not the same size\n";
        return false;
    }// if

    for (int i = 0 ; i < width * height * 4 ; i += 4)
    {
        unsigned char        rgb1[3];
        unsigned char        rgb2[3];

        RGBA_To_RGB(data + i, rgb1);
        RGBA_To_RGB(pImage->data + i, rgb2);

        data[i] = abs(rgb1[0] - rgb2[0]);
        data[i+1] = abs(rgb1[1] - rgb2[1]);
        data[i+2] = abs(rgb1[2] - rgb2[2]);
        data[i+3] = 255;
    }

    return true;
}// Difference



double* TargaImage::find_filter_value(const std::vector<std::vector<double>>& m, const int x, const int y) {    
    //determine the size of the filter
    int filter_size = m.size();
    int filter_half = filter_size / 2;

    double output[3] = { 0, 0, 0 };
    
    //for each pixel, find their filtered number by calculating its nearby pixels and the filter
    for (int i = 0; i < filter_size; ++i) {
        int mask_x = x + i - filter_half;
        if (mask_x < 0) mask_x = -1 * mask_x;
        if (mask_x > width - 1) mask_x = (width - 1) - (mask_x - (width - 1));

        for (int j = 0; j < filter_size; ++j) {
            int mask_y = y + j - filter_half;
            if (mask_y < 0) mask_y = -1 * mask_y;
            if (mask_y > height - 1) mask_y = (height - 1) - (mask_y - (height - 1));

            int current_pos = mask_y * width + mask_x;
            int rgb_offset = current_pos * 4;

            output[0] += data[rgb_offset] * m[j][i];
            output[1] += data[rgb_offset + 1] * m[j][i];
            output[2] += data[rgb_offset + 2] * m[j][i];
        }
    }
    return output;
}


double* TargaImage::find_filter_value(const std::vector<std::vector<double>>& m, const int x1, const int x2, const int y1, const int y2) {
    //determine the size of the filter
    int filter_size_x = m[0].size();
    int filter_size_y = m.size();
    int filter_half_x = filter_size_x / 2;
    int filter_half_y = filter_size_y / 2;

    double output[3] = { 0, 0, 0 };

    //for each pixel, find their filtered number by calculating its nearby pixels and the filter
    for (int i = 0; i < filter_size_x; ++i) {
        int mask_x = x1 + i;
        if (mask_x < 0) mask_x = -1 * mask_x;
        if (mask_x > width - 1) mask_x = (width - 1) - (mask_x - (width - 1));

        for (int j = 0; j < filter_size_y; ++j) {
            int mask_y = y1 + j;
            if (mask_y < 0) mask_y = -1 * mask_y;
            if (mask_y > height - 1) mask_y = (height - 1) - (mask_y - (height - 1));

            int current_pos = mask_y * width + mask_x;
            int rgb_offset = current_pos * 4;

            output[0] += data[rgb_offset] * m[j][i];
            output[1] += data[rgb_offset + 1] * m[j][i];
            output[2] += data[rgb_offset + 2] * m[j][i];
        }
    }
    return output;
}

///////////////////////////////////////////////////////////////////////////////
//
//use this function for any filter function with a filter mask
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::Filter_template(const vector<vector<double>>& m) {
    if (m.empty()) {
        return;
    }

    
    //use this temperary array to store the filtered data
    unsigned char* rgb_data = To_RGB();
    
    //run through every pixels
    for (int pixel_offset = 0; pixel_offset < width * height; ++pixel_offset) {
        int x_pos = pixel_offset % width;
        int y_pos = pixel_offset / width;

        double* temp_rgb = find_filter_value(m, x_pos, y_pos);

        
        //put the filtered data back into the real data array
        int rgba_offset = pixel_offset * 4;
        data[rgba_offset] = temp_rgb[0];
        data[rgba_offset + 1] = temp_rgb[1];
        data[rgba_offset + 2] = temp_rgb[2];
    }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 box filter on this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Box()
{
    vector<vector<double>> mask(5, vector<double>(5, 1.0 / 25.0));
    Filter_template(mask);

    return true;
}// Filter_Box


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Bartlett filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Bartlett()
{
    vector<vector<double>> mask = { {0.0123, 0.0247, 0.0370, 0.0247, 0.0123},
                            {0.0247, 0.0494, 0.0741, 0.0494, 0.0247},
                            {0.0370, 0.0741, 0.1111, 0.0741, 0.0370},
                            {0.0247, 0.0494, 0.0741, 0.0494, 0.0247},
                            {0.0123, 0.0247, 0.0370, 0.0247, 0.0123} };
    Filter_template(mask);

    return true;
}// Filter_Bartlett

///////////////////////////////////////////////////////////////////////////////
//
//return a gaussian mask with given size
//
///////////////////////////////////////////////////////////////////////////////
vector<vector<double>> TargaImage::gaussian_mask_maker(int size) {
    vector<vector<double>> m(size, vector<double>(size, 0));
    vector<double> row(size, 0);
    for (int i = 0; i < size; ++i) {
        int n = size - 1;
        row[i] = Binomial(n, i);
    }

    double divider = (double)(pow(pow(2, (size / 2)), 4));

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            m[i][j] = row[i] * row[j] / divider;
        }
    }

    return m;
}

///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Gaussian()
{
    vector<vector<double>> mask = gaussian_mask_maker(5);
    Filter_template(mask);

    return true;
}// Filter_Gaussian

///////////////////////////////////////////////////////////////////////////////
//
//      Perform NxN Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////

bool TargaImage::Filter_Gaussian_N( unsigned int N )
{
    vector<vector<double>> mask = gaussian_mask_maker(N);
    Filter_template(mask);

    return true;
}// Filter_Gaussian_N


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 edge detect (high pass) filter on this image.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Edge()
{
    ClearToBlack();
    return false;
}// Filter_Edge


///////////////////////////////////////////////////////////////////////////////
//
//      Perform a 5x5 enhancement filter to this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Enhance()
{
    ClearToBlack();
    return false;
}// Filter_Enhance




///////////////////////////////////////////////////////////////////////////////
//
//      Run simplified version of Hertzmann's painterly image filter.
//      You probably will want to use the Draw_Stroke funciton and the
//      Stroke class to help.
// Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::NPR_Paint()
{
    ClearToBlack();
    return false;
}



///////////////////////////////////////////////////////////////////////////////
//
//      Halve the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Half_Size()
{
    vector<vector<double>> bartlett = { 
        {0.0625, 0.125, 0.0625}, 
        {0.125, 0.25, 0.125} ,
        {0.0625, 0.125, 0.0625} 
    };

    int new_width = width * 0.5;
    int new_height = height * 0.5;
    int new_size = new_width * new_height;
    unsigned char* output_pixels = new unsigned char[new_size * 4];

    for (int i = 0; i < new_height; ++i) {
        for (int j = 0; j < new_width; ++j) {
            int origin_y = i * 2;
            int origin_x = j * 2;

            double* tempRGB = find_filter_value(bartlett, origin_x, origin_y);

            int new_offset = (i * new_width + j) * 4;
            int origin_offset = (origin_y * width + origin_x) * 4;
            output_pixels[new_offset] = tempRGB[0];
            output_pixels[new_offset+1] = tempRGB[1];
            output_pixels[new_offset+2] = tempRGB[2];
            output_pixels[new_offset + 3] = data[origin_offset + 3];
        }
    }
    width = new_width;
    height = new_height;

    delete[] data;
    data = output_pixels;

    return true;
}// Half_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Double the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Double_Size()
{
    int new_width = width * 2;
    int new_height = height * 2;
    int new_size = new_width * new_height;
    unsigned char* output_pixels = new unsigned char[new_size * 4];



    for (int i = 0; i < new_height; ++i) {
        for (int j = 0; j < new_width; ++j) {
            int origin_y = (double)i / 2.0;
            int origin_x = (double)j / 2.0;

            vector<vector<double>> bartlett;
            double* tempRGB;



            //if (i % 2 == 0 && j % 2 == 0) {
            //    bartlett = {
            //        {0.0625, 0.125, 0.0625}, 
            //        {0.125, 0.25, 0.125} ,
            //        {0.0625, 0.125, 0.0625} 
            //    };
            //    tempRGB = find_filter_value(bartlett, origin_x, origin_y);
            //}
            //else if (i % 2 != 0 && j % 2 != 0) {
            //    bartlett = {
            //        {0.015625, 0.046875, 0.046875, 0.015625},
            //        {0.046875, 0.140625, 0.140625, 0.046875},
            //        {0.046875, 0.140625, 0.140625, 0.046875},
            //        {0.015625, 0.046875, 0.046875, 0.015625}
            //    };
            //    tempRGB = find_filter_value(bartlett, origin_x - 1, origin_y - 1, origin_x + 2, origin_y + 2);
            //}
            //else {
            //    bartlett = {
            //        {0.03125, 0.0625, 0.03125},
            //        {0.09375, 0.1875, 0.09375},
            //        {0.09375, 0.1875, 0.09375},
            //        {0.03125, 0.0625, 0.03125}
            //    };
            //    tempRGB = find_filter_value(bartlett, origin_x - 1, origin_y - 1, origin_x + 1, origin_y + 2);
            //}

            
            bartlett = {
                {0.0625, 0.125, 0.0625},
                {0.125, 0.25, 0.125} ,
                {0.0625, 0.125, 0.0625}
            };
            tempRGB = find_filter_value(bartlett, origin_x, origin_y);

            int new_offset = (i * new_width + j) * 4;
            int origin_offset = (origin_y * width + origin_x) * 4;
            output_pixels[new_offset] = tempRGB[0];
            output_pixels[new_offset + 1] = tempRGB[1];
            output_pixels[new_offset + 2] = tempRGB[2];
            output_pixels[new_offset + 3] = data[origin_offset + 3];
        }
    }

    width = new_width;
    height = new_height;

    delete[] data;
    data = output_pixels;

    return true;
}// Double_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Scale the image dimensions by the given factor.  The given factor is 
//  assumed to be greater than one.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Resize(float scale)
{
    vector<vector<double>> bartlett = {
    {0.0625, 0.125, 0.0625},
    {0.125, 0.25, 0.125} ,
    {0.0625, 0.125, 0.0625}
    };

    int new_width = width * scale;
    int new_height = height * scale;
    int new_size = new_width * new_height;
    unsigned char* output_pixels = new unsigned char[new_size * 4];

    for (int i = 0; i < new_height; ++i) {
        for (int j = 0; j < new_width; ++j) {
            int origin_y = i / scale;
            int origin_x = j / scale;

            double* tempRGB = find_filter_value(bartlett, origin_x, origin_y);

            int new_offset = (i * new_width + j) * 4;
            int origin_offset = (origin_y * width + origin_x) * 4;
            output_pixels[new_offset] = tempRGB[0];
            output_pixels[new_offset + 1] = tempRGB[1];
            output_pixels[new_offset + 2] = tempRGB[2];
            output_pixels[new_offset + 3] = data[origin_offset + 3];
        }
    }
    width = new_width;
    height = new_height;

    delete[] data;
    data = output_pixels;

    return true;
}// Resize


//////////////////////////////////////////////////////////////////////////////
//
//      Rotate the image clockwise by the given angle.  Do not resize the 
//  image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Rotate(float angleDegrees)
{
    vector<vector<double>> bartlett = {
        {0.0625, 0.125, 0.0625},
        {0.125, 0.25, 0.125} ,
        {0.0625, 0.125, 0.0625}
    };

    double center_x = (double)width / 2.0;
    double center_y = (double)height / 2.0;

   

    unsigned char* output_pixels = new unsigned char[width * height * 4];

    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {

            //translate the output pixel position to the center spot
            double translated_j = (double)j - center_x;
            double translated_i = (double)i - center_y;

            //find the rotated original position
            double rotated_x = translated_j * cos(DegreesToRadians(-angleDegrees)) - translated_i * sin(DegreesToRadians(-angleDegrees));
            double rotated_y = translated_j * sin(DegreesToRadians(-angleDegrees)) + translated_i * cos(DegreesToRadians(-angleDegrees));

            //translate again to find the original position

            int origin_x = rotated_x + center_x;
            int origin_y = rotated_y + center_y;



            double* tempRGB;
            double black[3] = { 0.0, 0.0, 0.0 };

            int new_offset = (i * width + j) * 4;
            int origin_offset;

            if (origin_x < 0 || origin_x >= width || origin_y < 0 || origin_y >= height) {
                tempRGB = black;
                origin_offset = 0;
            }
            else {
                tempRGB = find_filter_value(bartlett, origin_x, origin_y);
                origin_offset = (origin_y * width + origin_x) * 4;
            }
            

            output_pixels[new_offset] = tempRGB[0];
            output_pixels[new_offset + 1] = tempRGB[1];
            output_pixels[new_offset + 2] = tempRGB[2];
            output_pixels[new_offset + 3] = data[origin_offset + 3];
        }
    }


    delete[] data;
    data = output_pixels;

    return true;
}// Rotate


//////////////////////////////////////////////////////////////////////////////
//
//      Given a single RGBA pixel return, via the second argument, the RGB
//      equivalent composited with a black background.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::RGBA_To_RGB(unsigned char *rgba, unsigned char *rgb)
{
    const unsigned char	BACKGROUND[3] = { 0, 0, 0 };

    unsigned char  alpha = rgba[3];

    if (alpha == 0)
    {
        rgb[0] = BACKGROUND[0];
        rgb[1] = BACKGROUND[1];
        rgb[2] = BACKGROUND[2];
    }
    else
    {
	    float	alpha_scale = (float)255 / (float)alpha;
	    int	val;
	    int	i;

	    for (i = 0 ; i < 3 ; i++)
	    {
	        val = (int)floor(rgba[i] * alpha_scale);
	        if (val < 0)
		    rgb[i] = 0;
	        else if (val > 255)
		    rgb[i] = 255;
	        else
		    rgb[i] = val;
	    }
    }
}// RGA_To_RGB


///////////////////////////////////////////////////////////////////////////////
//
//      Copy this into a new image, reversing the rows as it goes. A pointer
//  to the new image is returned.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Reverse_Rows(void)
{
    unsigned char   *dest = new unsigned char[width * height * 4];
    TargaImage	    *result;
    int 	        i, j;

    if (! data)
    	return NULL;

    for (i = 0 ; i < height ; i++)
    {
	    int in_offset = (height - i - 1) * width * 4;
	    int out_offset = i * width * 4;

	    for (j = 0 ; j < width ; j++)
        {
	        dest[out_offset + j * 4] = data[in_offset + j * 4];
	        dest[out_offset + j * 4 + 1] = data[in_offset + j * 4 + 1];
	        dest[out_offset + j * 4 + 2] = data[in_offset + j * 4 + 2];
	        dest[out_offset + j * 4 + 3] = data[in_offset + j * 4 + 3];
        }
    }

    result = new TargaImage(width, height, dest);
    delete[] dest;
    return result;
}// Reverse_Rows


///////////////////////////////////////////////////////////////////////////////
//
//      Clear the image to all black.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::ClearToBlack()
{
    memset(data, 0, width * height * 4);
}// ClearToBlack


///////////////////////////////////////////////////////////////////////////////
//
//      Helper function for the painterly filter; paint a stroke at
// the given location
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::Paint_Stroke(const Stroke& s) {
   int radius_squared = (int)s.radius * (int)s.radius;
   for (int x_off = -((int)s.radius); x_off <= (int)s.radius; x_off++) {
      for (int y_off = -((int)s.radius); y_off <= (int)s.radius; y_off++) {
         int x_loc = (int)s.x + x_off;
         int y_loc = (int)s.y + y_off;
         // are we inside the circle, and inside the image?
         if ((x_loc >= 0 && x_loc < width && y_loc >= 0 && y_loc < height)) {
            int dist_squared = x_off * x_off + y_off * y_off;
            if (dist_squared <= radius_squared) {
               data[(y_loc * width + x_loc) * 4 + 0] = s.r;
               data[(y_loc * width + x_loc) * 4 + 1] = s.g;
               data[(y_loc * width + x_loc) * 4 + 2] = s.b;
               data[(y_loc * width + x_loc) * 4 + 3] = s.a;
            } else if (dist_squared == radius_squared + 1) {
               data[(y_loc * width + x_loc) * 4 + 0] = 
                  (data[(y_loc * width + x_loc) * 4 + 0] + s.r) / 2;
               data[(y_loc * width + x_loc) * 4 + 1] = 
                  (data[(y_loc * width + x_loc) * 4 + 1] + s.g) / 2;
               data[(y_loc * width + x_loc) * 4 + 2] = 
                  (data[(y_loc * width + x_loc) * 4 + 2] + s.b) / 2;
               data[(y_loc * width + x_loc) * 4 + 3] = 
                  (data[(y_loc * width + x_loc) * 4 + 3] + s.a) / 2;
            }
         }
      }
   }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke() {}

///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke(unsigned int iradius, unsigned int ix, unsigned int iy,
               unsigned char ir, unsigned char ig, unsigned char ib, unsigned char ia) :
   radius(iradius),x(ix),y(iy),r(ir),g(ig),b(ib),a(ia)
{
}

