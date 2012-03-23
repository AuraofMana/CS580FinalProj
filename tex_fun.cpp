/* Texture functions for cs580 GzLib	*/
#include    "stdafx.h" 
#include	"stdio.h"
#include	"Gz.h"

GzColor	*image;
int xs, ys;
int reset = 1;

/* Image texture function */
int tex_fun(float u, float v, GzColor color)
{
  unsigned char		pixel[3];
  unsigned char     dummy;
  char  		foo[8];
  int   		i, j;
  FILE			*fd;

  if (reset) {          /* open and load texture file */
    fd = fopen ("texture", "rb");
    if (fd == NULL) {
      fprintf (stderr, "texture file not found\n");
      exit(-1);
    }
    fscanf (fd, "%s %d %d %c", foo, &xs, &ys, &dummy);
    image = (GzColor*)malloc(sizeof(GzColor)*(xs+1)*(ys+1));
    if (image == NULL) {
      fprintf (stderr, "malloc for texture image failed\n");
      exit(-1);
    }

    for (i = 0; i < xs*ys; i++) {	/* create array of GzColor values */
      fread(pixel, sizeof(pixel), 1, fd);
      image[i][RED] = (float)((int)pixel[RED]) * (1.0 / 255.0);
      image[i][GREEN] = (float)((int)pixel[GREEN]) * (1.0 / 255.0);
      image[i][BLUE] = (float)((int)pixel[BLUE]) * (1.0 / 255.0);
      }

    reset = 0;          /* init is done */
	fclose(fd);
  }

/* bounds-test u,v to make sure nothing will overflow image array bounds */
/* determine texture cell corner values and perform bilinear interpolation */
/* set color to interpolated GzColor value and return */

  if(u < 0.0f) u = 0.0f;
  if(v < 0.0f) v = 0.0f;
  if(u > 1.0f) u = 1.0f;
  if(v > 1.0f) v = 1.0f;

  int u0 = u * (xs - 1);
  int v0 = v * (ys - 1);

  int u1, v1;

  if(u0 == xs - 1) u1 = u0 - 1;
  else u1 = u0 + 1;
  if(v0 == ys - 1) v1 = v0 - 1;
  else v1 = v0 + 1;

  float s = (u * (xs - 1)) - float(u0);
  float t = (v * (ys - 1)) - float(v0);

  GzColor *topLeft, *topRight, *bottomRight, *bottomLeft;
  topLeft = &image[v0 * xs + u0];
  topRight = &image[v0 * xs + u1];
  bottomRight = &image[v1 * xs + u1];
  bottomLeft = &image[v1 * xs + u0];

  float topLeftRatio = (1 - s) * (1 - t);
  float topRightRatio = s * (1 - t);
  float bottomRightRatio = s * t;
  float bottomLeftRatio = (1 - s) * t;

  color[0] = bottomRightRatio * (*bottomRight)[0] + bottomLeftRatio * (*bottomLeft)[0] + topRightRatio * (*topRight)[0] + topLeftRatio * (*topLeft)[0];
  color[1] = bottomRightRatio * (*bottomRight)[1] + bottomLeftRatio * (*bottomLeft)[1] + topRightRatio * (*topRight)[1] + topLeftRatio * (*topLeft)[1];
  color[2] = bottomRightRatio * (*bottomRight)[2] + bottomLeftRatio * (*bottomLeft)[2] + topRightRatio * (*topRight)[2] + topLeftRatio * (*topLeft)[2];

  return GZ_SUCCESS;
}


/* Procedural texture function */
int ptex_fun(float u, float v, GzColor color)
{
  unsigned char		pixel[3];
  unsigned char     dummy;
  char  		foo[8];
  int   		i, j, c;
  FILE			*fd;

  xs = ys = 64;

  if (reset) { 
	  image = (GzColor*)malloc(sizeof(GzColor)*xs*ys);
	  if (image == NULL) {
		  fprintf (stderr, "malloc for texture image failed\n");
		  exit(-1);
	  }

	  //Source: http://glprogramming.com/red/chapter09.html
	  for (i = 0; i < ys; i++) {
		  for (j = 0; j < xs; j++) {
			  c = ((((i&0x8)==0)^((j&0x8))==0));
			  int imgIndex = j * (xs) + i;
			  image[imgIndex][0] = c;
			  image[imgIndex][1] = c;
			  image[imgIndex][2] = c;
		  }
	  }
	  reset = 0;          /* init is done */
  }

/* bounds-test u,v to make sure nothing will overflow image array bounds */
/* determine texture cell corner values and perform bilinear interpolation */
/* set color to interpolated GzColor value and return */

  if(u < 0.0f || u > 1.0f || v < 0.0f || v > 1.0f)
  {
	  return GZ_FAILURE;
  }

  int u0 = u * (xs - 1);
  int v0 = v * (ys - 1);

  int u1, v1;

  if(u0 == xs - 1) u1 = u0 - 1;
  else u1 = u0 + 1;
  if(v0 == ys - 1) v1 = v0 - 1;
  else v1 = v0 + 1;

  float s = (u * (xs - 1)) - float(u0);
  float t = (v * (ys - 1)) - float(v0);

  GzColor *topLeft, *topRight, *bottomRight, *bottomLeft;
  topLeft = &image[v0 * xs + u0];
  topRight = &image[v0 * xs + u1];
  bottomRight = &image[v1 * xs + u1];
  bottomLeft = &image[v1 * xs + u0];

  float topLeftRatio = (1 - s) * (1 - t);
  float topRightRatio = s * (1 - t);
  float bottomRightRatio = s * t;
  float bottomLeftRatio = (1 - s) * t;

  color[0] = bottomRightRatio * (*bottomRight)[0] + bottomLeftRatio * (*bottomLeft)[0] + topRightRatio * (*topRight)[0] + topLeftRatio * (*topLeft)[0];
  color[1] = bottomRightRatio * (*bottomRight)[1] + bottomLeftRatio * (*bottomLeft)[1] + topRightRatio * (*topRight)[1] + topLeftRatio * (*topLeft)[1];
  color[2] = bottomRightRatio * (*bottomRight)[2] + bottomLeftRatio * (*bottomLeft)[2] + topRightRatio * (*topRight)[2] + topLeftRatio * (*topLeft)[2];

  return GZ_SUCCESS;
}

