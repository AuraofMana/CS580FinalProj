#pragma once

#include "disp.h" /* include your own disp.h file (e.g. hw1)*/

/* Camera defaults */
#define	DEFAULT_FOV		35.0
#define	DEFAULT_IM_Z	(-10.0)  /* world coords for image plane origin */
#define	DEFAULT_IM_Y	(5.0)    /* default look-at point = 0,0,0 */
#define	DEFAULT_IM_X	(-10.0)

#define	DEFAULT_AMBIENT	{0.1, 0.1, 0.1}
#define	DEFAULT_DIFFUSE	{0.7, 0.6, 0.5}
#define	DEFAULT_SPECULAR	{0.2, 0.3, 0.4}
#define	DEFAULT_SPEC		32

#define	MATLEVELS	100		/* how many matrix pushes allowed */
#define	MAX_LIGHTS	10		/* how many lights allowed */

#define BB_MINX -5;
#define BB_MAXX 5;
#define BB_MINY -5;
#define BB_MAXY 5;
#define BB_MINZ -5;
#define BB_MAXZ 5;

//Stereo Camera Display Index defines
#define ACTUALDISPLAY 0
#define STEREOLEFT 1
#define STEREORIGHT 2

struct CubeMap
{
	int xSize, ySize;
	GzColor *posX;
	GzColor *negX;
	GzColor *posY;
	GzColor *negY;
	GzColor *posZ;
	GzColor *negZ;
};

enum CUBEMAPSIDE
{
	LEFT,
	RIGHT,
	UP,
	DOWN,
	FRONT,
	BACK
};

#ifndef GZRENDER
#define GZRENDER
typedef struct {			/* define a renderer */
  GzRenderClass	renderClass;
  GzDisplay		*display[3]; //Actual display + left and right display
  short		    open;
  GzCamera		camera, leftCamera, rightCamera;
  short		    matlevel;	        /* top of stack - current xform */
  GzMatrix		Ximage[MATLEVELS];	/* stack of xforms (Xsm) */
  GzMatrix		Xnorm[MATLEVELS];	/* xforms for norms (Xim) */
  GzMatrix		Xsp;		        /* NDC to screen (pers-to-screen) */
  GzColor		flatcolor;          /* color state for flat shaded triangles */
  int			interp_mode;
  int			numlights;
  GzLight		lights[MAX_LIGHTS];
  GzLight		ambientlight;
  GzColor		Ka, Kd, Ks;
  float		    spec;		/* specular power */
  GzTexture		tex_fun;    /* tex_fun(float u, float v, GzColor color) */
  CubeMap		cmap;
}  GzRender;
#endif

// Function declaration
// HW2
int GzNewRender(GzRender **render, GzRenderClass renderClass, GzDisplay *display);
int GzFreeRender(GzRender *render);
int GzBeginRender(GzRender	*render);
int GzPutAttribute(GzRender	*render, int numAttributes, GzToken	*nameList, 
	GzPointer *valueList);
int GzPutTriangle(GzRender *render, int	numParts, GzToken *nameList,
	GzPointer *valueList);

//HW2 Added functions
short ctoi(float color); //Defined by original code but not declared

//Bilinear interpolation based on: http://ezekiel.vancouver.wsu.edu/~cs442/archive/lectures/gouraud/gouraud.pdf
#ifndef GZEDGE
#define GZEDGE
typedef struct
{
	GzCoord start;
	GzCoord end;
	GzCoord current;
	float slopeX;
	float slopeZ;

	GzCoord dataStart;
	GzCoord dataEnd;
	float ratio;

	GzTextureIndex texStart;
	GzTextureIndex texEnd;
} GzEdge;
#endif

#ifndef GZSPAN
#define GZSPAN
typedef struct 
{
	float start[2];
	float end[2];
	float current[2];
	float slopeZ;

	GzCoord dataStart;
	GzCoord dataEnd;
	float ratio;

	GzTextureIndex texStart;
	GzTextureIndex texEnd;
} GzSpan;
#endif

int GzVertSorting(GzCoord *vertices, GzCoord *normals, GzTextureIndex *textures);
int GzCeiling(float fvar);

// HW3
int GzPutCamera(GzRender *render, GzCamera *camera);
int GzPushMatrix(GzRender *render, GzMatrix	matrix);
int GzPopMatrix(GzRender *render);

// Object Translation
int GzRotXMat(float degree, GzMatrix mat);
int GzRotYMat(float degree, GzMatrix mat);
int GzRotZMat(float degree, GzMatrix mat);
int GzTrxMat(GzCoord translate, GzMatrix mat);
int GzScaleMat(GzCoord scale, GzMatrix mat);

//HW3 Added functions
typedef float   GzVector[4];

void GzVectorDiff(const GzCoord &v1, const GzCoord &v2, GzCoord &diff);
float GzDotProduct(const GzCoord &v1, const GzCoord &v2);
void GzCrossProduct(const GzCoord &v1, const GzCoord &v2, GzCoord &res);
float GzVectorMagnitude(const GzCoord &vec);
void GzNormalizeVector(const GzCoord &vec, GzCoord &res);
void GzSetVector(const GzCoord &vec, GzCoord &res);
void GzMultiplyVector(const GzCoord &vec, const float scalar, GzCoord &res);
void GzAddVector(const GzCoord &v1, const GzCoord &v2, GzCoord &res);
void GzSubtractVector(const GzCoord &v1, const GzCoord &v2, GzCoord &res);
void GzMatrixMultiplication(const GzMatrix &m1, const GzMatrix &m2, GzMatrix &res);
void GzMatrixTimesVector(const GzMatrix &m, const GzVector &vec, GzVector &res);
void GzCoordToGzVector(const GzCoord &vec, GzVector &res);
void GzVectorToGzCoord(const GzVector &vec, GzCoord &res);
void GzConcatMatrix(GzRender *render, GzMatrix &res);

//HW4 Added functions
void GzCalculateColor(const GzRender *render, const GzCoord &normal, GzColor &color, boolean gtexture);
void GzColorMultiply(const GzColor &color1, const GzColor &color2, GzColor &res);
//Can't use const & because API code sucks
typedef float   Gz3x3Matrix[3][3];
void GzMatrixTo3x3(GzMatrix orig, Gz3x3Matrix &res);
void GzMatrix3x3TimesScalar(const Gz3x3Matrix &m, float scalar, Gz3x3Matrix &res);
void GzMatrix3x3Transpose(const Gz3x3Matrix &m, Gz3x3Matrix &res);
void GzConcatMatrixNormal(GzRender *render, GzMatrix &res);

//HW5 Added functions
float GzNewVz(float currZ);
void GzXformToPerspective(const GzTextureIndex &texi, float Vz, GzTextureIndex &res);
void GzXformToAffine(const GzTextureIndex &texi, float Vz, GzTextureIndex &res);

//Final Project Added functions
//Cube Mapping
void GzLoadCubeMaps(GzRender *render);
void GzGetCubeMapColor(GzRender *render, const GzCoord &normal, GzColor &color);
void GzGetCubeMapTexture(GzRender *render, CUBEMAPSIDE cmEnum, float u, float v, GzColor &color);

//Stereoscopic 3D
void GzCopyCamera(const GzCamera &cameraSrc, GzCamera &cameraDest);
void GzLoadXiw(GzCamera &camera);
void GzStereoInit(GzRender *render, const GzCoord &leftPos, const GzCoord &rightPos);
void GzInsertXiw(GzRender *render, GzMatrix matrix);
int GzStereoPutTriangle(GzRender *render, int numParts, GzToken *nameList, GzPointer *valueList);
int GzStereoPutTriangleHelper(GzRender *render, GzCoord vertices[3], GzCoord normals[3], GzTextureIndex textures[3], bool leftCamera);
void GzCombineDisplays(GzRender *render);