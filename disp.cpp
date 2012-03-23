/*   CS580 HW   */
#include    "stdafx.h"  
#include	"Gz.h"
#include	"disp.h"

#include <climits>

#define GzClamp(intensity) (((intensity) < (0)) ? (0) : (((intensity) > (4095)) ? (4095) : (intensity)))

int GzNewFrameBuffer(char** framebuffer, int width, int height)
{
/* create a framebuffer:
 -- allocate memory for framebuffer : (sizeof)GzPixel x width x height
 -- pass back pointer 
*/

	*framebuffer = new char[3 * width * height];

	return GZ_SUCCESS;
}

int GzNewDisplay(GzDisplay **display, GzDisplayClass dispClass, int xRes, int yRes)
{
/* create a display:
  -- allocate memory for indicated class and resolution
  -- pass back pointer to GzDisplay object in display
*/

	*display = new GzDisplay;
	(*display)->xres = xRes;
	(*display)->yres = yRes;
	(*display)->dispClass = dispClass;
	(*display)->open = 1;
	(*display)->fbuf = new GzPixel[xRes * yRes];
	return GZ_SUCCESS;
}


int GzFreeDisplay(GzDisplay	*display)
{
/* clean up, free memory */

	delete [] display->fbuf;
	delete display;
	return GZ_SUCCESS;
}


int GzGetDisplayParams(GzDisplay *display, int *xRes, int *yRes, GzDisplayClass	*dispClass)
{
/* pass back values for an open display */

	*xRes = display->xres;
	*yRes = display->yres;
	*dispClass = display->dispClass;
	return GZ_SUCCESS;
}


int GzInitDisplay(GzDisplay	*display)
{
/* set everything to some default values - start a new frame */

	const int RESOLUTION = display->xres * display->yres;
	GzPixel initPixel = {4095, 4095, 4095, 255, INT_MAX};
	for(int i = 0; i < RESOLUTION; ++i)
	{
		display->fbuf[i] = initPixel;
	}

	return GZ_SUCCESS;
}


int GzPutDisplay(GzDisplay *display, int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z)
{
/* write pixel values into the display */

	GzPixel currPixel = {GzClamp(r), GzClamp(g), GzClamp(b), a, z};
	//Bound checking
	if(display == 0 || i < 0 || i >= display->xres || j < 0 || j >= display->yres) return GZ_FAILURE; //Ignore anything outside of the resolution
	display->fbuf[j * display->xres + i] = currPixel;
	return GZ_SUCCESS;
}


int GzGetDisplay(GzDisplay *display, int i, int j, GzIntensity *r, GzIntensity *g, GzIntensity *b, GzIntensity *a, GzDepth *z)
{
	/* pass back pixel value in the display */
	/* check display class to see what vars are valid */

	//Bound checking
	if(display == 0 || i < 0 || i >= display->xres || j < 0 || j >= display->yres) return GZ_FAILURE;
	GzPixel currPixel = display->fbuf[j * display->xres + i];
	*r = currPixel.red;
	*g = currPixel.green;
	*b = currPixel.blue;
	*a = currPixel.alpha;
	*z = currPixel.z;

	return GZ_SUCCESS;
}

int GzFlushDisplay2File(FILE* outfile, GzDisplay *display)
{
	/* write pixels to ppm file based on display class -- "P6 %d %d 255\r" */
	
	//Print header
	fprintf(outfile, "P6 %d %d 255\n", display->xres, display->yres); //"\n" for Linux, "\r" for Mac, "\r\n" for Windows
	//Print body
	const int BUFFERSIZE = display->xres * display->yres * 3;
	unsigned char *buffer = new unsigned char[BUFFERSIZE]; //Not a string, don't need +1 for '\0'

	GzIntensity r, g, b, a;
	GzDepth z;

	for(int i = 0; i < display->xres; ++i)
	{
		for(int j = 0; j < display->yres; ++j)
		{
			GzGetDisplay(display, i, j, &r, &g, &b, &a, &z);
			int currentIndex = (j * display->xres + i) * 3;
			buffer[currentIndex] = static_cast<unsigned char>(r >> 4);
			buffer[currentIndex + 1] = static_cast<unsigned char>(g >> 4);
			buffer[currentIndex + 2] = static_cast<unsigned char>(b >> 4);
		}
	}

	fwrite(buffer, sizeof(char), BUFFERSIZE, outfile);
	delete [] buffer;
	
	return GZ_SUCCESS;
}

int GzFlushDisplay2FrameBuffer(char* framebuffer, GzDisplay *display)
{
	/* write pixels to framebuffer: 
		- Put the pixels into the frame buffer
		- Caution: store the pixel to the frame buffer as the order of blue, green, and red 
		- Not red, green, and blue !!!
	*/

	delete [] framebuffer;
	framebuffer = new char[display->xres * display->yres * 3];

	GzIntensity r, g, b, a;
	GzDepth z;
	
	for(int i = 0; i < display->xres; ++i)
	{
		for(int j = 0; j < display->yres; ++j)
		{
			GzGetDisplay(display, i, j, &r, &g, &b, &a, &z);
			int currentIndex = (j * display->xres + i) * 3;
			framebuffer[currentIndex] = static_cast<unsigned char>(b >> 4);
			framebuffer[currentIndex + 1] = static_cast<unsigned char>(g >> 4);
			framebuffer[currentIndex + 2] = static_cast<unsigned char>(r >> 4);
		}
	}

	return GZ_SUCCESS;
}