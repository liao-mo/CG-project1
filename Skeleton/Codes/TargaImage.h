///////////////////////////////////////////////////////////////////////////////
//
//      TargaImage.h                            Author:     Stephen Chenney
//                                              Modified:   Eric McDaniel
//                                              Date:       Fall 2004
//
//      Class to manipulate targa images.  You must implement the image 
//  modification functions.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef _TARGA_IMAGE_H_
#define _TARGA_IMAGE_H_

#include <Fl/Fl.h>
#include <Fl/Fl_Widget.h>
#include <stdio.h>
#include <string>
#include <vector>

class Stroke;
class DistanceImage;

class TargaImage
{
    // methods
    public:
	    TargaImage(void);
        TargaImage(int w, int h);
	    TargaImage(int w, int h, unsigned char *d);
        TargaImage(const TargaImage& image);
	    ~TargaImage(void);

        unsigned char*	To_RGB(void);	            // Convert the image to RGB format,
        float* get_float_data(void);
        bool Save_Image(const char*);               // save the image to a file
        static TargaImage* Load_Image(char*);       // Load a file and return a pointer to a new TargaImage object.  Returns NULL on failure

        bool To_Grayscale();

        bool Quant_Uniform();
        bool Quant_Populosity();
        bool Quant_Median();

        bool Dither_Threshold();
        bool Dither_Random();
        bool Dither_FS();
        bool Dither_Bright();
        bool Dither_Cluster();
        bool Dither_Color();

        bool Comp_Over(TargaImage* pImage);
        bool Comp_In(TargaImage* pImage);
        bool Comp_Out(TargaImage* pImage);
        bool Comp_Atop(TargaImage* pImage);
        bool Comp_Xor(TargaImage* pImage);

        bool Difference(TargaImage* pImage);

        void Filter_template(const std::vector<std::vector<double> >& m);
        bool Filter_Box();
        bool Filter_Bartlett();
        std::vector<std::vector<double> > gaussian_mask_maker(int size);
        bool Filter_Gaussian();
        bool Filter_Gaussian_N(unsigned int N);
        bool Filter_Edge();
        bool Filter_Enhance();
        double* find_filter_value(const std::vector<std::vector<double>>& m, const int x, const int y);
        double* find_filter_value(const std::vector<std::vector<double>>& m, const int x1, const int x2, const int y1, const int y2);

        bool NPR_Paint();

        bool Half_Size();
        bool Double_Size();
        bool Resize(float scale);
        bool Rotate(float angleDegrees);

    private:
	// helper function for format conversion
        void RGBA_To_RGB(unsigned char *rgba, unsigned char *rgb);

        // reverse the rows of the image, some targas are stored bottom to top
	TargaImage* Reverse_Rows(void);

	// clear image to all black
        void ClearToBlack();

	// Draws a filled circle according to the stroke data
        void Paint_Stroke(const Stroke& s);

    // members
    public:
        int		width;	    // width of the image in pixels
        int		height;	    // height of the image in pixels
        unsigned char* data;	    // pixel data for the image, assumed to be in pre-multiplied RGBA format.
        unsigned char* original_data;
        double* gray_data;

};

class Stroke { // Data structure for holding painterly strokes.
public:
   Stroke(void);
   Stroke(unsigned int radius, unsigned int x, unsigned int y,
          unsigned char r, unsigned char g, unsigned char b, unsigned char a);
   
   // data
   unsigned int radius, x, y;	// Location for the stroke
   unsigned char r, g, b, a;	// Color
};

struct color_data {
    color_data() : r_val(0), g_val(0), b_val(0) {};
    unsigned char r_val;
    unsigned char g_val;
    unsigned char b_val;
    bool operator==(const color_data& c) const{
        return (r_val == c.r_val && g_val == c.g_val && b_val == c.b_val);
    }
    bool operator<(const color_data& c) const {
        if (r_val < c.r_val) return true;
        else return false;
    }
};

struct color_data_hasher {
    std::size_t operator()(const color_data& c) const
    {
        using std::size_t;
        using std::hash;
        using std::string;

        return ((hash<unsigned char>()(c.r_val)^ (hash<unsigned char>()(c.g_val) << 1)) >> 1)^ (hash<unsigned char>()(c.b_val) << 1);
    }
};

#endif


