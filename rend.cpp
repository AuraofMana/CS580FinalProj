/* CS580 Homework 3 */

#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"

#include <algorithm>
#include <cmath>
#include <string>

using std::copy;
using std::string;

#define GZDEGREETORADIAN(GzDegree) (GzDegree * 3.14159265 / 180.0)

int GzRotXMat(float degree, GzMatrix mat)
{
// Create rotate matrix : rotate along x axis
// Pass back the matrix using mat value

	double radian = GZDEGREETORADIAN(degree);
	float sinRes = sin(radian);
	float cosRes = cos(radian);
	mat[0][0] = 1.0;
	mat[0][1] = 0.0;
	mat[0][2] = 0.0;
	mat[0][3] = 0.0;

	mat[1][0] = 0.0;
	mat[1][1] = cosRes;
	mat[1][2] = -sinRes;
	mat[1][3] = 0.0;

	mat[2][0] = 0.0;
	mat[2][1] = sinRes;
	mat[2][2] = cosRes;
	mat[2][3] = 0.0;

	mat[3][0] = 0.0;
	mat[3][1] = 0.0;
	mat[3][2] = 0.0;
	mat[3][3] = 1.0;

	return GZ_SUCCESS;
}


int GzRotYMat(float degree, GzMatrix mat)
{
// Create rotate matrix : rotate along y axis
// Pass back the matrix using mat value

	double radian = GZDEGREETORADIAN(degree);
	float sinRes = sin(radian);
	float cosRes = cos(radian);
	mat[0][0] = cosRes;
	mat[0][1] = 0.0;
	mat[0][2] = sinRes;
	mat[0][3] = 0.0;

	mat[1][0] = 0.0;
	mat[1][1] = 1.0;
	mat[1][2] = 0.0;
	mat[1][3] = 0.0;

	mat[2][0] = -sinRes;
	mat[2][1] = 0.0;
	mat[2][2] = cosRes;
	mat[2][3] = 0.0;

	mat[3][0] = 0.0;
	mat[3][1] = 0.0;
	mat[3][2] = 0.0;
	mat[3][3] = 1.0;

	return GZ_SUCCESS;
}


int GzRotZMat(float degree, GzMatrix mat)
{
// Create rotate matrix : rotate along z axis
// Pass back the matrix using mat value

	double radian = GZDEGREETORADIAN(degree);
	float sinRes = sin(radian);
	float cosRes = cos(radian);
	mat[0][0] = cosRes;
	mat[0][1] = -sinRes;
	mat[0][2] = 0.0;
	mat[0][3] = 0.0;

	mat[1][0] = sinRes;
	mat[1][1] = cosRes;
	mat[1][2] = 0.0;
	mat[1][3] = 0.0;

	mat[2][0] = 0.0;
	mat[2][1] = 0.0;
	mat[2][2] = 1.0;
	mat[2][3] = 0.0;

	mat[3][0] = 0.0;
	mat[3][1] = 0.0;
	mat[3][2] = 0.0;
	mat[3][3] = 1.0;

	return GZ_SUCCESS;
}


int GzTrxMat(GzCoord translate, GzMatrix mat)
{
// Create translation matrix
// Pass back the matrix using mat value

	mat[0][0] = 1.0;
	mat[0][1] = 0.0;
	mat[0][2] = 0.0;
	mat[0][3] = translate[0];

	mat[1][0] = 0.0;
	mat[1][1] = 1.0;
	mat[1][2] = 0.0;
	mat[1][3] = translate[1];

	mat[2][0] = 0.0;
	mat[2][1] = 0.0;
	mat[2][2] = 1.0;
	mat[2][3] = translate[2];

	mat[3][0] = 0.0;
	mat[3][1] = 0.0;
	mat[3][2] = 0.0;
	mat[3][3] = 1.0;

	return GZ_SUCCESS;
}


int GzScaleMat(GzCoord scale, GzMatrix mat)
{
// Create scaling matrix
// Pass back the matrix using mat value

	mat[0][0] = scale[0];
	mat[0][1] = 0.0;
	mat[0][2] = 0.0;
	mat[0][3] = 0.0;

	mat[1][0] = 0.0;
	mat[1][1] = scale[1];
	mat[1][2] = 0.0;
	mat[1][3] = 0.0;

	mat[2][0] = 0.0;
	mat[2][1] = 0.0;
	mat[2][2] = scale[2];
	mat[2][3] = 0.0;

	mat[3][0] = 0.0;
	mat[3][1] = 0.0;
	mat[3][2] = 0.0;
	mat[3][3] = 1.0;

	return GZ_SUCCESS;
}


//----------------------------------------------------------
// Begin main functions

int GzNewRender(GzRender **render, GzRenderClass renderClass, GzDisplay	*display)
{
/*  
- malloc a renderer struct 
- keep closed until all inits are done 
- setup Xsp and anything only done once 
- span interpolator needs pointer to display 
- check for legal class GZ_Z_BUFFER_RENDER 
- init default camera 
*/

	//Set up Xsp
	float XsOverTwo = float(display->xres / 2);
	float YsOverTwo = float(display->yres / 2);
	float ZmaxOverD = INT_MAX * tan(GZDEGREETORADIAN(DEFAULT_FOV / 2.0));

	*render = new GzRender;
	GzRender tempRender =
	{
		renderClass, 
		display, 
		1,
		{
			{0.0}, //Transform from world to image space
			{0.0}, //Perspective projection transform
			{DEFAULT_IM_X, DEFAULT_IM_Y, DEFAULT_IM_Z}, //Position of image plane origin
			{0.0, 0.0, 0.0}, //Position of look-at-point
			{0.0, 1.0, 0.0}, //World up-vector (almost screen up)
			DEFAULT_FOV //Horizontal field of view
		}, //Camera
		-1, //-1 to denote that there is nothing on the stack
		{0.0}, 
		{0.0}, 
		{
			{XsOverTwo, 0.0, 0.0, XsOverTwo},
			{0.0, -YsOverTwo, 0.0, YsOverTwo},
			{0.0, 0.0, ZmaxOverD, 0.0},
			{0.0, 0.0, 0.0, 1.0}
		}, 
		{GZ_RGB_COLOR, GZ_RGB_COLOR, GZ_RGB_COLOR}, 
		0, 
		0, 
		{0}, 
		{0}, 
		DEFAULT_AMBIENT, 
		DEFAULT_DIFFUSE, 
		DEFAULT_SPECULAR,
		DEFAULT_SPEC, 
		0
	};
	**render = tempRender;

	LoadCubeMaps(*render);

	return GZ_SUCCESS;
}


int GzFreeRender(GzRender *render)
{
/* 
-free all renderer resources
*/
	free(render->cmap.posX);
	free(render->cmap.negX);
	free(render->cmap.posY);
	free(render->cmap.negY);
	free(render->cmap.posZ);
	free(render->cmap.negZ);

	delete render;
	return GZ_SUCCESS;
}


int GzBeginRender(GzRender *render)
{
/*  
- set up for start of each frame - clear frame buffer 
- compute Xiw and projection xform Xpi from camera definition 
- init Ximage - put Xsp at base of stack, push on Xpi and Xiw 
- now stack contains Xsw and app can push model Xforms if it want to. 
*/
	
	//Clear frame buffer
	GzPixel initPixel = {2048, 1792, 1536, 255, INT_MAX};
	for(int i = 0; i < render->display->xres * render->display->yres; ++i)
	{
		render->display->fbuf[i] = initPixel;
	}

	//Xiw
	GzCoord cl;
	GzVectorDiff(render->camera.lookat, render->camera.position, cl);
	
	//New Z
	GzCoord newZ;
	GzNormalizeVector(cl, newZ);

	//New Up
	float newZdotUp = GzDotProduct(render->camera.worldup, newZ);

	GzCoord temp;
	GzMultiplyVector(newZ, newZdotUp, temp);

	GzSubtractVector(render->camera.worldup, temp, temp);

	GzSetVector(temp, render->camera.worldup);

	//New Y
	GzCoord newY;
	GzNormalizeVector(temp, newY);

	//New X
	GzCoord newX;
	GzCrossProduct(newY, newZ, newX);

	//-X, -Y, -Z dot C
	GzCoord negNewX, negNewY, negNewZ;
	GzMultiplyVector(newX, -1.0, negNewX);
	GzMultiplyVector(newY, -1.0, negNewY);
	GzMultiplyVector(newZ, -1.0, negNewZ);

	float negXdotC = GzDotProduct(negNewX, render->camera.position);
	float negYdotC = GzDotProduct(negNewY, render->camera.position);
	float negZdotC = GzDotProduct(negNewZ, render->camera.position);

	render->camera.Xiw[0][0] = newX[0];
	render->camera.Xiw[0][1] = newX[1];
	render->camera.Xiw[0][2] = newX[2];
	render->camera.Xiw[0][3] = negXdotC;

	render->camera.Xiw[1][0] = newY[0];
	render->camera.Xiw[1][1] = newY[1];
	render->camera.Xiw[1][2] = newY[2];
	render->camera.Xiw[1][3] = negYdotC;

	render->camera.Xiw[2][0] = newZ[0];
	render->camera.Xiw[2][1] = newZ[1];
	render->camera.Xiw[2][2] = newZ[2];
	render->camera.Xiw[2][3] = negZdotC;

	render->camera.Xiw[3][0] = 0.0;
	render->camera.Xiw[3][1] = 0.0;
	render->camera.Xiw[3][2] = 0.0;
	render->camera.Xiw[3][3] = 1.0;


	//Xpi
	float oneOverD = tan(GZDEGREETORADIAN(render->camera.FOV / 2.0));

	render->camera.Xpi[0][0] = 1.0;
	render->camera.Xpi[0][1] = 0.0;
	render->camera.Xpi[0][2] = 0.0;
	render->camera.Xpi[0][3] = 0.0;

	render->camera.Xpi[1][0] = 0.0;
	render->camera.Xpi[1][1] = 1.0;
	render->camera.Xpi[1][2] = 0.0;
	render->camera.Xpi[1][3] = 0.0;

	render->camera.Xpi[2][0] = 0.0;
	render->camera.Xpi[2][1] = 0.0;
	render->camera.Xpi[2][2] = 1.0;
	render->camera.Xpi[2][3] = 0.0;

	render->camera.Xpi[3][0] = 0.0;
	render->camera.Xpi[3][1] = 0.0;
	render->camera.Xpi[3][2] = oneOverD;
	render->camera.Xpi[3][3] = 1.0;


	//init Ximage - put Xsp at base of stack, push on Xpi and Xiw 
	int status = 0;
	
	status |= GzPushMatrix(render, render->Xsp);
	status |= GzPushMatrix(render, render->camera.Xpi);
	status |= GzPushMatrix(render, render->camera.Xiw);
	if(status) return GZ_FAILURE;
	else return GZ_SUCCESS;
}

int GzPutCamera(GzRender *render, GzCamera *camera)
{
/*
- overwrite renderer camera structure with new camera definition
*/

	copy(&camera->Xiw[0][0], &camera->Xiw[0][0] + 16, &render->camera.Xiw[0][0]);
	copy(&camera->Xpi[0][0], &camera->Xpi[0][0] + 16, &render->camera.Xpi[0][0]);
	copy(&camera->position[0], &camera->position[0] + 3, &render->camera.position[0]);
	copy(&camera->lookat[0], &camera->lookat[0] + 3, &render->camera.lookat[0]);
	copy(&camera->worldup[0], &camera->worldup[0] + 3, &render->camera.worldup[0]);
	GzNormalizeVector(render->camera.worldup, render->camera.worldup); //Normalize up just in case
	render->camera.FOV = camera->FOV;

	render->Xsp[2][2] = INT_MAX * tan(GZDEGREETORADIAN(render->camera.FOV / 2.0));

	return GZ_SUCCESS;
}

int GzPushMatrix(GzRender *render, GzMatrix	matrix)
{
/*
- push a matrix onto the Ximage stack
- check for stack overflow
*/

	if(render->matlevel == MATLEVELS - 1) return GZ_FAILURE;
	++(render->matlevel);
	copy(&matrix[0][0], &matrix[0][0] + 16, &render->Ximage[render->matlevel][0][0]);

	GzMatrix concat = {0.0};
	copy(&matrix[0][0], &matrix[0][0] + 16, &concat[0][0]);

	Gz3x3Matrix newMatrix = {0.0};
	GzMatrixTo3x3(concat, newMatrix);

	if(render->matlevel < 2)
	{
		concat[0][0] = 1.0;
		concat[0][1] = 0.0;
		concat[0][2] = 0.0;
		concat[0][3] = 0.0;

		concat[1][0] = 0.0;
		concat[1][1] = 1.0;
		concat[1][2] = 0.0;
		concat[1][3] = 0.0;

		concat[2][0] = 0.0;
		concat[2][1] = 0.0;
		concat[2][2] = 1.0;
		concat[2][3] = 0.0;

		concat[3][0] = 0.0;
		concat[3][1] = 0.0;
		concat[3][2] = 0.0;
		concat[3][3] = 1.0;

		copy(&concat[0][0], &concat[0][0] + 16, &render->Xnorm[render->matlevel][0][0]);

		return GZ_SUCCESS;
	}
	else if(render->matlevel > 2)
	{
		//Form unitary rotations
		float URScale = 1.0 / sqrt(newMatrix[0][0] * newMatrix[0][0] + newMatrix[0][1] * newMatrix[0][1] + newMatrix[0][2] * newMatrix[0][2]);
		GzMatrix3x3TimesScalar(newMatrix, URScale, newMatrix);
	}

	//Invert the matrix
	GzCoord x0 = {0.0};
	GzCoord x1 = {0.0};
	GzCoord x2 = {0.0};
	GzCoord x1crossx2 = {0.0};
	GzCoord x2crossx0 = {0.0};
	GzCoord x0crossx1 = {0.0};

	x0[0] = newMatrix[0][0];
	x0[1] = newMatrix[0][1];
	x0[2] = newMatrix[0][2];

	x1[0] = newMatrix[1][0];
	x1[1] = newMatrix[1][1];
	x1[2] = newMatrix[1][2];

	x2[0] = newMatrix[2][0];
	x2[1] = newMatrix[2][1];
	x2[2] = newMatrix[2][2];

	GzCrossProduct(x1, x2, x1crossx2);
	GzCrossProduct(x2, x0, x2crossx0);
	GzCrossProduct(x0, x1, x0crossx1);

	float det = GzDotProduct(x0, x1crossx2);
	det = 1.0 / det;

	Gz3x3Matrix invTemp = {0.0};
	invTemp[0][0] = x1crossx2[0];
	invTemp[1][0] = x1crossx2[1];
	invTemp[2][0] = x1crossx2[2];
	invTemp[0][1] = x2crossx0[0];
	invTemp[1][1] = x2crossx0[1];
	invTemp[2][1] = x2crossx0[2];
	invTemp[0][2] = x0crossx1[0];
	invTemp[1][2] = x0crossx1[1];
	invTemp[2][2] = x0crossx1[2];

	GzMatrix3x3TimesScalar(invTemp, det, newMatrix);

	GzMatrix3x3Transpose(newMatrix, newMatrix);

	//Convert 3x3 back to GzMatrix
	concat[0][0] = newMatrix[0][0];
	concat[0][1] = newMatrix[0][1];
	concat[0][2] = newMatrix[0][2];
	concat[0][3] = 0.0;

	concat[1][0] = newMatrix[1][0];
	concat[1][1] = newMatrix[1][1];
	concat[1][2] = newMatrix[1][2];
	concat[1][3] = 0.0;

	concat[2][0] = newMatrix[2][0];
	concat[2][1] = newMatrix[2][1];
	concat[2][2] = newMatrix[2][2];
	concat[2][3] = 0.0;

	concat[3][0] = 0.0;
	concat[3][1] = 0.0;
	concat[3][2] = 0.0;
	concat[3][3] = 1.0;

	copy(&concat[0][0], &concat[0][0] + 16, &render->Xnorm[render->matlevel][0][0]);

	return GZ_SUCCESS;
}

int GzPopMatrix(GzRender *render)
{
/*
- pop a matrix off the Ximage stack
- check for stack underflow
*/

	if(render->matlevel < 0) return GZ_FAILURE;
	--(render->matlevel);
	return GZ_SUCCESS;
}


int GzPutAttribute(GzRender	*render, int numAttributes, GzToken	*nameList, 
	GzPointer	*valueList) /* void** valuelist */
{
/*
- set renderer attribute states (e.g.: GZ_RGB_COLOR default color)
- later set shaders, interpolaters, texture maps, and lights
*/
	for(int i = 0; i < numAttributes; ++i)
	{
		switch(nameList[i])
		{
		case GZ_DIRECTIONAL_LIGHT:
			{
				render->lights[i] = *((GzLight *) valueList[i]);
				break;
			}
		case GZ_AMBIENT_LIGHT:
			{
				render->ambientlight = *((GzLight *) valueList[i]);
				break;
			}
		case GZ_AMBIENT_COEFFICIENT:
			{
				GzColor *valueListColor = (GzColor *) valueList[i];
				render->Ka[0] = (*valueListColor)[0];
				render->Ka[1] = (*valueListColor)[1];
				render->Ka[2] = (*valueListColor)[2];
				break;
			}
		case GZ_SPECULAR_COEFFICIENT:
			{
				GzColor *valueListColor = (GzColor *) valueList[i];
				render->Ks[0] = (*valueListColor)[0];
				render->Ks[1] = (*valueListColor)[1];
				render->Ks[2] = (*valueListColor)[2];
				break;
			}
		case GZ_DIFFUSE_COEFFICIENT:
			{
				GzColor *valueListColor = (GzColor *) valueList[i];
				render->Kd[0] = (*valueListColor)[0];
				render->Kd[1] = (*valueListColor)[1];
				render->Kd[2] = (*valueListColor)[2];
				break;
			}
		case GZ_DISTRIBUTION_COEFFICIENT:
			{
				render->spec = *((float *) valueList[i]);
				break;
			}
		case GZ_INTERPOLATE:
			{
				render->interp_mode = *((int *) valueList[i]);
				break;
			}
		case GZ_TEXTURE_MAP:
			{
				render->tex_fun = *((GzTexture) valueList[i]);
				break;
			}
		}
	}

	return GZ_SUCCESS;
}

int GzPutTriangle(GzRender	*render, int numParts, GzToken *nameList, 
				  GzPointer	*valueList)
/* numParts : how many names and values */
{
/*  
- pass in a triangle description with tokens and values corresponding to 
      GZ_POSITION:3 vert positions in model space 
- Xform positions of verts  
- Clip - just discard any triangle with verts behind view plane 
       - test for triangles with all three verts off-screen 
- invoke triangle rasterizer  
*/

	GzCoord vertices[3] = {0};
	GzCoord normals[3] = {0};
	GzTextureIndex textures[3] = {0};

	for(int i = 0; i < numParts; ++i)
	{
		switch(nameList[i])
		{
		case GZ_POSITION:
			{
				GzCoord *valueListTriangles = (GzCoord *) valueList[i];
				vertices[0][0] = valueListTriangles[0][0];
				vertices[0][1] = valueListTriangles[0][1];
				vertices[0][2] = valueListTriangles[0][2];
				vertices[1][0] = valueListTriangles[1][0];
				vertices[1][1] = valueListTriangles[1][1];
				vertices[1][2] = valueListTriangles[1][2];
				vertices[2][0] = valueListTriangles[2][0];
				vertices[2][1] = valueListTriangles[2][1];
				vertices[2][2] = valueListTriangles[2][2];
				break;
			}
		case GZ_NORMAL:
			{
				GzCoord *valueListNormals = (GzCoord *) valueList[i];
				normals[0][0] = valueListNormals[0][0];
				normals[0][1] = valueListNormals[0][1];
				normals[0][2] = valueListNormals[0][2];
				normals[1][0] = valueListNormals[1][0];
				normals[1][1] = valueListNormals[1][1];
				normals[1][2] = valueListNormals[1][2];
				normals[2][0] = valueListNormals[2][0];
				normals[2][1] = valueListNormals[2][1];
				normals[2][2] = valueListNormals[2][2];
				break;
			}
		case GZ_TEXTURE_INDEX:
			{
				GzTextureIndex *valueListTextures = (GzTextureIndex *) valueList[i];
				textures[0][0] = valueListTextures[0][0];
				textures[0][1] = valueListTextures[0][1];
				textures[1][0] = valueListTextures[1][0];
				textures[1][1] = valueListTextures[1][1];
				textures[2][0] = valueListTextures[2][0];
				textures[2][1] = valueListTextures[2][1];
				break;
			}
		}
	}


	//Grab the concatenation matrix for vertices
	GzMatrix concat = {0.0};
	GzConcatMatrix(render, concat);

	//Multiply the vertices with the concatenated matrix
	GzVector triVec;
	//Vertex 0
	GzCoordToGzVector(vertices[0], triVec);
	GzMatrixTimesVector(concat, triVec, triVec);
	GzVectorToGzCoord(triVec, vertices[0]);
	GzMultiplyVector(vertices[0], 1.0 / triVec[3], vertices[0]); //Perspective

	//Vertex 1
	GzCoordToGzVector(vertices[1], triVec);
	GzMatrixTimesVector(concat, triVec, triVec);
	GzVectorToGzCoord(triVec, vertices[1]);
	GzMultiplyVector(vertices[1], 1.0 / triVec[3], vertices[1]); //Perspective

	//Vertex 2
	GzCoordToGzVector(vertices[2], triVec);
	GzMatrixTimesVector(concat, triVec, triVec);
	GzVectorToGzCoord(triVec, vertices[2]);
	GzMultiplyVector(vertices[2], 1.0 / triVec[3], vertices[2]); //Perspective

	//Clip
	if(vertices[0][2] < 0.0 || vertices[1][2] < 0.0 || vertices[2][2] < 0.0)
	{
		return GZ_SUCCESS; //Done
	}


	//Grab the concatenation matrix for normals
	GzConcatMatrixNormal(render, concat);

	GzVector normVec;

	//Normal 0
	GzCoordToGzVector(normals[0], normVec);
	GzMatrixTimesVector(concat, normVec, normVec);
	GzVectorToGzCoord(normVec, normals[0]);

	//Normal 1
	GzCoordToGzVector(normals[1], normVec);
	GzMatrixTimesVector(concat, normVec, normVec);
	GzVectorToGzCoord(normVec, normals[1]);

	//Normal 2
	GzCoordToGzVector(normals[2], normVec);
	GzMatrixTimesVector(concat, normVec, normVec);
	GzVectorToGzCoord(normVec, normals[2]);


	//Transform the texture indices to perspective space

	//Tex ID 0
	GzXformToPerspective(textures[0], GzNewVz(vertices[0][2]), textures[0]);

	//Tex ID 1
	GzXformToPerspective(textures[1], GzNewVz(vertices[1][2]), textures[1]);

	//Tex ID 2
	GzXformToPerspective(textures[2], GzNewVz(vertices[2][2]), textures[2]);


	//Step 1: Sort vertices by Y
	int result = GzVertSorting(vertices, normals, textures);

	//Step 2: Setup edge DDAs
	/*
	Edge0: 1 -> 0; Edge1: 2 -> 1; Edge 2: 2 -> 0
	*/
	GzEdge edge10 = 
	{
		{vertices[1][0], vertices[1][1], vertices[1][2]},
		{vertices[0][0], vertices[0][1], vertices[0][2]},
		{vertices[1][0], vertices[1][1], vertices[1][2]},
		(vertices[0][0] - vertices[1][0]) / (vertices[0][1] - vertices[1][1]),
		(vertices[0][2] - vertices[1][2]) / (vertices[0][1] - vertices[1][1]),
		{0.0},
		{0.0},
		0.0,
		{textures[1][0], textures[1][1]},
		{textures[0][0], textures[0][1]}
	};

	GzEdge edge21 =
	{
		{vertices[2][0], vertices[2][1], vertices[2][2]},
		{vertices[1][0], vertices[1][1], vertices[1][2]},
		{vertices[2][0], vertices[2][1], vertices[2][2]},
		(vertices[1][0] - vertices[2][0]) / (vertices[1][1] - vertices[2][1]),
		(vertices[1][2] - vertices[2][2]) / (vertices[1][1] - vertices[2][1]),
		{0.0},
		{0.0},
		0.0,
		{textures[2][0], textures[2][1]},
		{textures[1][0], textures[1][1]}
	};

	GzEdge edge20 =
	{
		{vertices[2][0], vertices[2][1], vertices[2][2]},
		{vertices[0][0], vertices[0][1], vertices[0][2]},
		{vertices[2][0], vertices[2][1], vertices[2][2]},
		(vertices[0][0] - vertices[2][0]) / (vertices[0][1] - vertices[2][1]),
		(vertices[0][2] - vertices[2][2]) / (vertices[0][1] - vertices[2][1]),
		{0.0},
		{0.0},
		0.0,
		{textures[2][0], textures[2][1]},
		{textures[0][0], textures[0][1]}
	};

	if(render->interp_mode == GZ_COLOR)
	{
		GzNormalizeVector(normals[0], normals[0]);
		GzNormalizeVector(normals[1], normals[1]);
		GzNormalizeVector(normals[2], normals[2]);

		GzColor n1 = {0.0}, n0 = {0.0}, n2 = {0.0};
		if(render->tex_fun != 0)
		{
			GzCalculateColor(render, normals[0], n0, true);
			GzCalculateColor(render, normals[1], n1, true);
			GzCalculateColor(render, normals[2], n2, true);
		}
		else
		{
			GzCalculateColor(render, normals[0], n0, false);
			GzCalculateColor(render, normals[1], n1, false);
			GzCalculateColor(render, normals[2], n2, false);
		}

		copy(&n0[0], &n0[0] + 3, &edge10.dataEnd[0]);
		copy(&n0[0], &n0[0] + 3, &edge20.dataEnd[0]);

		copy(&n1[0], &n1[0] + 3, &edge10.dataStart[0]);
		copy(&n1[0], &n1[0] + 3, &edge21.dataEnd[0]);

		copy(&n2[0], &n2[0] + 3, &edge20.dataStart[0]);
		copy(&n2[0], &n2[0] + 3, &edge21.dataStart[0]);
	}
	else if(render->interp_mode == GZ_NORMAL)
	{
		copy(&normals[1][0], &normals[1][0] + 3, &(edge10.dataStart[0]));
		copy(&normals[0][0], &normals[0][0] + 3, &(edge10.dataEnd[0]));

		copy(&normals[2][0], &normals[2][0] + 3, &(edge21.dataStart[0]));
		copy(&normals[1][0], &normals[1][0] + 3, &(edge21.dataEnd[0]));

		copy(&normals[2][0], &normals[2][0] + 3, &(edge20.dataStart[0]));
		copy(&normals[0][0], &normals[0][0] + 3, &(edge20.dataEnd[0]));
	}
	else
	{
		GzCoord avgNormal = 
		{
			(normals[0][0] + normals[1][0] + normals[2][0]) / 3.0,
			(normals[0][1] + normals[1][1] + normals[2][1]) / 3.0,
			(normals[0][2] + normals[1][2] + normals[2][2]) / 3.0
		};

		GzCalculateColor(render, avgNormal, render->flatcolor, false);
	}

	//Step 3 & 4: Sort edges by L and R
	/*
	if result = 0 then L & R will be E1 & E2 (E1 may not be left)
	result = 0 has a special case where it's just a line
	if result = 1 then L & R will be E0 & E2 (E0 may not be left)
	if result = 2 then first L & R (starting from top) will be E1 & E2 (E1 may not be left)
	and second L & R will be E0 & E2 (E0 may not be left, but E2 continues to be left or right)

	if result = 0 then L slope < R slope
	if result = 1 then L slope > R slope
	if result = 2 then L slope < R slope
	*/
	GzEdge *leftEdge, *rightEdge;
	if(result != 1) //Y0 = Y1 != Y2 OR Y0 != Y1 != Y2 OR Y0 = Y1 = Y2
	{
		if(edge21.slopeX < edge20.slopeX)
		{
			leftEdge = &edge21;
			rightEdge = &edge20;
		}
		else
		{
			leftEdge = &edge20;
			rightEdge = &edge21;
		}

		//Step 5: Advance DDA current positions to top y-scan line (ceiling)
		//You know the top is vertex 2
		float deltaY = GzCeiling(vertices[2][1]) - vertices[2][1];
		leftEdge->current[0] += (leftEdge->slopeX * deltaY);
		leftEdge->current[1] += deltaY;
		leftEdge->current[2] += (leftEdge->slopeZ * deltaY);
		rightEdge->current[0] += (rightEdge->slopeX * deltaY);
		rightEdge->current[1] += deltaY;
		rightEdge->current[2] += (rightEdge->slopeZ * deltaY);

		leftEdge->ratio = deltaY / (leftEdge->end[1] - leftEdge->start[1]);
		rightEdge->ratio = deltaY / (rightEdge->end[1] - rightEdge->start[1]);

		if(result == 2)
		{
			float deltaYEdge10 = GzCeiling(edge10.start[1]) - edge10.start[1];
			edge10.current[0] += (edge10.slopeX * deltaYEdge10);
			edge10.current[1] += deltaYEdge10;
			edge10.current[2] += (edge10.slopeZ * deltaYEdge10);
			
			edge10.ratio = deltaYEdge10 / (edge10.end[1] - edge10.start[1]);
		}
	}
	else//Y0 != Y1 = Y2
	{
		if(edge10.slopeX > edge20.slopeX)
		{
			leftEdge = &edge10;
			rightEdge = &edge20;
		}
		else
		{
			leftEdge = &edge20;
			rightEdge = &edge10;
		}

		//Step 5: Advance DDA current positions to top y-scan line (ceiling)
		float deltaYLeft = GzCeiling(leftEdge->start[1]) - leftEdge->start[1];
		float deltaYRight = GzCeiling(rightEdge->start[1]) - rightEdge->start[1];
		leftEdge->current[0] += (leftEdge->slopeX * deltaYLeft);
		leftEdge->current[1] += deltaYLeft;
		leftEdge->current[2] += (leftEdge->slopeZ * deltaYLeft);
		rightEdge->current[0] += (rightEdge->slopeX * deltaYRight);
		rightEdge->current[1] += deltaYRight;
		rightEdge->current[2] += (rightEdge->slopeZ * deltaYRight);

		leftEdge->ratio = deltaYLeft / (leftEdge->end[1] - leftEdge->start[1]);
		rightEdge->ratio = deltaYRight / (rightEdge->end[1] - rightEdge->start[1]);
	}

	//result = 3 => left edge ends first
	//result = 4 => right edge ends first
	//result = 5 => no edge to switch to
	if(result == 2) //Triangles with no horizontal edge must be tested to see which edge ends first
	{
		if(leftEdge->end[1] < rightEdge->end[1]) result = 3;
		else if(leftEdge->end[1] > rightEdge->end[1]) result = 4;
		else result = 5; //Passed in just a line then no need to switch edge
	}
	else result = 5; //No edge to switch because there is a horizontal edge

	GzSpan currSpan = {0};

	//Step 13: Continue spans until vert 3 is passed
	while(1)
	{
		//Step 11: Test for ending edge and switch if need be
		if(leftEdge->current[1] > leftEdge->end[1] || rightEdge->current[1] > rightEdge->end[1])
		{
			if(result == 5) break;
			else if(result == 3)
			{
				leftEdge = &edge10;
				//Check if edges change over an integer y
				if(leftEdge->start[1] == int(leftEdge->start[1]))
				{
					leftEdge->current[0] += leftEdge->slopeX;
					++leftEdge->current[1];
					leftEdge->current[2] += leftEdge->slopeZ;
					rightEdge->current[0] += rightEdge->slopeX;
					++rightEdge->current[1];
					rightEdge->current[2] += rightEdge->slopeZ;

					leftEdge->ratio = (leftEdge->current[1] - leftEdge->start[1]) / (leftEdge->end[1] - leftEdge->start[1]);
					rightEdge->ratio = (rightEdge->current[1] - rightEdge->start[1]) / (rightEdge->end[1] - rightEdge->start[1]);
				}
			}
			else
			{
				rightEdge = &edge10;
				//Check if edges change over an integer y
				if(rightEdge->start[1] == int(rightEdge->start[1]))
				{
					leftEdge->current[0] += leftEdge->slopeX;
					++leftEdge->current[1];
					leftEdge->current[2] += leftEdge->slopeZ;
					rightEdge->current[0] += rightEdge->slopeX;
					++rightEdge->current[1];
					rightEdge->current[2] += rightEdge->slopeZ;

					leftEdge->ratio = (leftEdge->current[1] - leftEdge->start[1]) / (leftEdge->end[1] - leftEdge->start[1]);
					rightEdge->ratio = (rightEdge->current[1] - rightEdge->start[1]) / (rightEdge->end[1] - rightEdge->start[1]);
				}
			}
			result = 5;
			continue;
		}

		//Step 6 & 7: Setup span DDA and set its current and end positions to left and right values, respectively
		//Or update current span
		currSpan.start[0] = leftEdge->current[0];
		currSpan.start[1] = leftEdge->current[2];
		currSpan.end[0] = rightEdge->current[0];
		currSpan.end[1] = rightEdge->current[2];
		currSpan.current[0] = leftEdge->current[0];
		currSpan.current[1] = leftEdge->current[2];
		currSpan.slopeZ = (rightEdge->current[2] - leftEdge->current[2]) / (rightEdge->current[0] - leftEdge->current[0]);

		float oneMinusRatio = 1.0 - leftEdge->ratio;;
		currSpan.dataStart[0] = leftEdge->dataStart[0] * oneMinusRatio + leftEdge->dataEnd[0] * leftEdge->ratio;
		currSpan.dataStart[1] = leftEdge->dataStart[1] * oneMinusRatio + leftEdge->dataEnd[1] * leftEdge->ratio;
		currSpan.dataStart[2] = leftEdge->dataStart[2] * oneMinusRatio + leftEdge->dataEnd[2] * leftEdge->ratio;
		currSpan.texStart[0] = leftEdge->texStart[0] * oneMinusRatio + leftEdge->texEnd[0] * leftEdge->ratio;
		currSpan.texStart[1] = leftEdge->texStart[1] * oneMinusRatio + leftEdge->texEnd[1] * leftEdge->ratio;

		oneMinusRatio = 1.0 - rightEdge->ratio;
		currSpan.dataEnd[0] = rightEdge->dataStart[0] * oneMinusRatio + rightEdge->dataEnd[0] * rightEdge->ratio;
		currSpan.dataEnd[1] = rightEdge->dataStart[1] * oneMinusRatio + rightEdge->dataEnd[1] * rightEdge->ratio;
		currSpan.dataEnd[2] = rightEdge->dataStart[2] * oneMinusRatio + rightEdge->dataEnd[2] * rightEdge->ratio;
		currSpan.texEnd[0] = rightEdge->texStart[0] * oneMinusRatio + rightEdge->texEnd[0] * rightEdge->ratio;
		currSpan.texEnd[1] = rightEdge->texStart[1] * oneMinusRatio + rightEdge->texEnd[1] * rightEdge->ratio;

		currSpan.ratio = 0.0;

		//Step 8: Advance span current position to leftmost covered pixel (ceiling)
		float deltaX = GzCeiling(leftEdge->current[0]) - leftEdge->current[0];
		currSpan.current[0] += deltaX;
		currSpan.current[1] += (currSpan.slopeZ * deltaX);
		currSpan.ratio = deltaX / (currSpan.end[0] - currSpan.start[0]);

		//Step 9: Interpolate span position and parameters (Z) until current position > end
		while(currSpan.current[0] <= currSpan.end[0])
		{
			//Step 10: Test interpolated-Z against FB-Z for each pixel - low Z wins
			GzIntensity r, g, b, a;
			GzDepth z;
			GzGetDisplay(render->display, currSpan.current[0], leftEdge->current[1], &r, &g, &b, &a, &z);

			if(currSpan.current[1] < z)
			{
				GzColor currColor = {0.0}, texColor = {0.0};
				GzTextureIndex currTex = {0.0};
				float oneMinusRatio = 1.0 - currSpan.ratio;;

				//Calculate texture color
				currTex[0] = currSpan.texStart[0] * oneMinusRatio + currSpan.texEnd[0] * currSpan.ratio;
				currTex[1] = currSpan.texStart[1] * oneMinusRatio + currSpan.texEnd[1] * currSpan.ratio;
				GzXformToAffine(currTex, GzNewVz(currSpan.current[1]), currTex);
				if(render->tex_fun != 0) render->tex_fun(currTex[0], currTex[1], texColor);

				if(render->interp_mode == GZ_COLOR)
				{
					currColor[0] = currSpan.dataStart[0] * oneMinusRatio + currSpan.dataEnd[0] * currSpan.ratio;
					currColor[1] = currSpan.dataStart[1] * oneMinusRatio + currSpan.dataEnd[1] * currSpan.ratio;
					currColor[2] = currSpan.dataStart[2] * oneMinusRatio + currSpan.dataEnd[2] * currSpan.ratio;

					if(render->tex_fun != 0) GzColorMultiply(texColor, currColor, currColor);
					GzPutDisplay(render->display, currSpan.current[0], leftEdge->current[1], ctoi(currColor[0]), ctoi(currColor[1]), ctoi(currColor[2]), a, currSpan.current[1]);
				}
				else if(render->interp_mode == GZ_NORMAL)
				{
					GzCoord currNormal = {0.0};

					currNormal[0] = currSpan.dataStart[0] * oneMinusRatio + currSpan.dataEnd[0] * currSpan.ratio;
					currNormal[1] = currSpan.dataStart[1] * oneMinusRatio + currSpan.dataEnd[1] * currSpan.ratio;
					currNormal[2] = currSpan.dataStart[2] * oneMinusRatio + currSpan.dataEnd[2] * currSpan.ratio;
					GzNormalizeVector(currNormal, currNormal);

					GetCubeMapColor(render, currNormal, texColor);

					if(render->tex_fun != 0)
					{
						copy(&texColor[0], &texColor[0] + 3, &(render->Kd[0]));
						copy(&texColor[0], &texColor[0] + 3, &(render->Ka[0]));
					}

					GzCalculateColor(render, currNormal, currColor, false);
					GzPutDisplay(render->display, currSpan.current[0], leftEdge->current[1], ctoi(currColor[0]), ctoi(currColor[1]), ctoi(currColor[2]), a, currSpan.current[1]);
				}
				else
				{
					GzPutDisplay(render->display, currSpan.current[0], leftEdge->current[1], ctoi(render->flatcolor[0]), ctoi(render->flatcolor[1]), ctoi(render->flatcolor[2]), a, currSpan.current[1]);
				}
			}
			++currSpan.current[0];
			currSpan.current[1] += currSpan.slopeZ;
			currSpan.ratio = (currSpan.current[0] - currSpan.start[0]) / (currSpan.end[0] - currSpan.start[0]);
		}

		//Advance current for left and right edges with deltaY = 1
		leftEdge->current[0] += leftEdge->slopeX;
		++leftEdge->current[1];
		leftEdge->current[2] += leftEdge->slopeZ;
		rightEdge->current[0] += rightEdge->slopeX;
		++rightEdge->current[1];
		rightEdge->current[2] += rightEdge->slopeZ;

		leftEdge->ratio = (leftEdge->current[1] - leftEdge->start[1]) / (leftEdge->end[1] - leftEdge->start[1]);
		rightEdge->ratio = (rightEdge->current[1] - rightEdge->start[1]) / (rightEdge->end[1] - rightEdge->start[1]);
	}
	return GZ_SUCCESS;
}

/* NOT part of API - just for general assistance */

short	ctoi(float color)		/* convert float color to GzIntensity short */
{
  return(short)((int)(color * ((1 << 12) - 1)));
}

int GzVertSorting(GzCoord *vertices, GzCoord *normals, GzTextureIndex *textures)
{
	using std::swap;
	if(vertices[0][1] < vertices[1][1])
	{
		swap(vertices[0], vertices[1]);
		swap(normals[0], normals[1]);
		swap(textures[0], textures[1]);
	}
	if(vertices[0][1] < vertices[2][1])
	{
		swap(vertices[0], vertices[2]);
		swap(normals[0], normals[2]);
		swap(textures[0], textures[2]);
	}
	if(vertices[1][1] < vertices[2][1])
	{
		swap(vertices[1], vertices[2]);
		swap(normals[1], normals[2]);
		swap(textures[1], textures[2]);
	}

	if(vertices[0][1] == vertices[2][1]) //Y0 = Y1 = Y2
	{
		return -1; //Do nothing
	}

	if(vertices[0][1] == vertices[1][1]) //Y0 = Y1 != Y2
	{
		return 0; //Y0 and Y1 are the smallest
	}
	else if(vertices[1][1] == vertices[2][1]) //Y0 != Y1 = Y2
	{
		return 1; //Y0 and Y1 are the smallest
	}
	else //Y0 != Y1 != Y2
	{
		return 2;
	}
}

int GzCeiling(float fvar)
{
	return fvar - int(fvar) > 0 ? int(fvar) + 1 : fvar;
}

void GzVectorDiff(const GzCoord &v1, const GzCoord &v2, GzCoord &diff)
{
	diff[0] = v1[0] - v2[0];
	diff[1] = v1[1] - v2[1];
	diff[2] = v1[2] - v2[2];
}

float GzDotProduct(const GzCoord &v1, const GzCoord &v2)
{
	return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

void GzCrossProduct(const GzCoord &v1, const GzCoord &v2, GzCoord &res)
{
	GzCoord temp;
	temp[0] = v1[1] * v2[2] - v2[1] * v1[2];
	temp[1] = v2[0] * v1[2] - v1[0] * v2[2];
	temp[2] = v1[0] * v2[1] - v2[0] * v1[1];
	res[0] = temp[0];
	res[1] = temp[1];
	res[2] = temp[2];
}

float GzVectorMagnitude(const GzCoord &vec)
{
	return sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}

void GzNormalizeVector(const GzCoord &vec, GzCoord &res)
{
	float length = GzVectorMagnitude(vec);
	res[0] = vec[0] / length;
	res[1] = vec[1] / length;
	res[2] = vec[2] / length;
}

void GzSetVector(const GzCoord &vec, GzCoord &res)
{
	if(&vec == &res) return;
	copy(&vec[0], &vec[0] + 3, &res[0]);
}

void GzMultiplyVector(const GzCoord &vec, const float scalar, GzCoord &res)
{
	res[0] = vec[0] * scalar;
	res[1] = vec[1] * scalar;
	res[2] = vec[2] * scalar;
}

void GzAddVector(const GzCoord &v1, const GzCoord &v2, GzCoord &res)
{
	res[0] = v1[0] + v2[0];
	res[1] = v1[1] + v2[1];
	res[2] = v1[2] + v2[2];
}

void GzSubtractVector(const GzCoord &v1, const GzCoord &v2, GzCoord &res)
{
	res[0] = v1[0] - v2[0];
	res[1] = v1[1] - v2[1];
	res[2] = v1[2] - v2[2];
}

void GzMatrixMultiplication(const GzMatrix &m1, const GzMatrix &m2, GzMatrix &res)
{
	GzMatrix temp = {0.0};
	for(int i = 0; i < 4; ++i)
	{
		for(int j = 0; j < 4; ++j)
		{
			for(int k = 0; k < 4; ++k)
			{
				temp[i][j] += m1[i][k] * m2[k][j];
			}
		}
	}
	copy(&temp[0][0], &temp[0][0] + 16, &res[0][0]);
}

void GzMatrixTimesVector(const GzMatrix &m, const GzVector &vec, GzVector &res)
{
	GzVector temp = {0.0};
	for(int i = 0; i < 4; ++i)
	{
		for(int j = 0; j < 4; ++j)
		{
			temp[i] += m[i][j] * vec[j];
		}
	}
	copy(&temp[0], &temp[0] + 4, &res[0]);
}

void GzCoordToGzVector(const GzCoord &vec, GzVector &res)
{
	res[0] = vec[0];
	res[1] = vec[1];
	res[2] = vec[2];
	res[3] = 1.0;
}

void GzVectorToGzCoord(const GzVector &vec, GzCoord &res)
{
	res[0] = vec[0];
	res[1] = vec[1];
	res[2] = vec[2];
}

void GzConcatMatrix(GzRender *render, GzMatrix &res)
{
	if(render->matlevel < 0) return;

	//Grab the first matrix
	copy(&(render->Ximage[render->matlevel][0][0]), &(render->Ximage[render->matlevel][0][0]) + 16, &(res[0][0]));

	int savedMatLevel = render->matlevel; //Save the current mat level

	GzPopMatrix(render);
	
	if(render->matlevel < 0) return;

	//Keep on grabbing
	while(render->matlevel != -1)
	{
		GzMatrixMultiplication(render->Ximage[render->matlevel], res, res);
		GzPopMatrix(render);
	}

	render->matlevel = savedMatLevel; //Recall old matlevel
}

void GzCalculateColor(const GzRender *render, const GzCoord &normal, GzColor &color, boolean gtexture)
{
	GzColor specularColor = {0.0};
	GzColor diffuseColor = {0.0};
	GzColor ambientColor = {0.0};

	for(int i = 0; i < MAX_LIGHTS; ++i)
	{
		//Specular
		GzCoord specularR;
		float NdotL = GzDotProduct(normal, render->lights[i].direction);
		GzMultiplyVector(normal, 2 * NdotL, specularR);
		GzSubtractVector(specularR, render->lights[i].direction, specularR);

		GzCoord specularE = {0.0, 0.0, -1.0};

		float RdotE = GzDotProduct(specularR, specularE);

		float NdotE = GzDotProduct(normal, specularE);
		if((NdotL < 0 && NdotE > 0) || (NdotL > 0 && NdotE < 0)) continue;
		else if(NdotL < 0 && NdotE < 0)
		{
			//New Normal
			GzCoord newNormal;
			GzMultiplyVector(normal, -1.0, newNormal);

			//New Specular
			NdotL = GzDotProduct(newNormal, render->lights[i].direction);
			GzMultiplyVector(newNormal, 2 * NdotL, specularR);
			GzSubtractVector(specularR, render->lights[i].direction, specularR);

			//New RdotE
			RdotE = GzDotProduct(specularR, specularE);
		}

		if(RdotE < 0) RdotE = 0.0;
		RdotE = pow(RdotE, render->spec);

		GzColor currentSpecularColor;
		GzMultiplyVector(render->lights[i].color, RdotE, currentSpecularColor);

		GzAddVector(specularColor, currentSpecularColor, specularColor);


		//Diffuse
		GzColor currentDiffuseColor;
		GzMultiplyVector(render->lights[i].color, NdotL, currentDiffuseColor);

		GzAddVector(diffuseColor, currentDiffuseColor, diffuseColor);
	}

	if(!gtexture)
	{
		GzColorMultiply(specularColor, render->Ks, specularColor);
		GzColorMultiply(diffuseColor, render->Kd, diffuseColor);
		GzColorMultiply(render->ambientlight.color, render->Ka, ambientColor);
	}
	else
	{
		ambientColor[0] = render->ambientlight.color[0];
		ambientColor[1] = render->ambientlight.color[1];
		ambientColor[2] = render->ambientlight.color[2];
	}

	color[0] = specularColor[0] + diffuseColor[0] + ambientColor[0];
	color[1] = specularColor[1] + diffuseColor[1] + ambientColor[1];
	color[2] = specularColor[2] + diffuseColor[2] + ambientColor[2];

	//Prevent color overflow
	if(color[0] > 1.0) color[0] = 1.0;
	if(color[1] > 1.0) color[1] = 1.0;
	if(color[2] > 1.0) color[2] = 1.0;
}

void GzColorMultiply(const GzColor &color1, const GzColor &color2, GzColor &res)
{
	res[0] = color1[0] * color2[0];
	res[1] = color1[1] * color2[1];
	res[2] = color1[2] * color2[2];
}

void GzMatrixTo3x3(GzMatrix orig, Gz3x3Matrix &res)
{
	res[0][0] = orig[0][0];
	res[0][1] = orig[0][1];
	res[0][2] = orig[0][2];

	res[1][0] = orig[1][0];
	res[1][1] = orig[1][1];
	res[1][2] = orig[1][2];

	res[2][0] = orig[2][0];
	res[2][1] = orig[2][1];
	res[2][2] = orig[2][2];
}

void GzMatrix3x3TimesScalar(const Gz3x3Matrix &m, float scalar, Gz3x3Matrix &res)
{
	res[0][0] = m[0][0] * scalar;
	res[0][1] = m[0][1] * scalar;
	res[0][2] = m[0][2] * scalar;

	res[1][0] = m[1][0] * scalar;
	res[1][1] = m[1][1] * scalar;
	res[1][2] = m[1][2] * scalar;

	res[2][0] = m[2][0] * scalar;
	res[2][1] = m[2][1] * scalar;
	res[2][2] = m[2][2] * scalar;
}

void GzMatrix3x3Transpose(const Gz3x3Matrix &m, Gz3x3Matrix &res)
{
	float temp;
	for(int i = 0; i < 3; ++i)
	{
		for(int j = i + 1; j < 3; ++j)
		{
			temp = m[i][j];
			res[i][j] = m[j][i];
			res[j][i] = temp;
		}
	}
}

void GzConcatMatrixNormal(GzRender *render, GzMatrix &res)
{
	if(render->matlevel < 0) return;

	//Grab the first matrix
	copy(&(render->Xnorm[render->matlevel][0][0]), &(render->Xnorm[render->matlevel][0][0]) + 16, &(res[0][0]));

	int savedMatLevel = render->matlevel; //Save the current mat level

	GzPopMatrix(render);

	if(render->matlevel < 0) return;

	//Keep on grabbing
	while(render->matlevel != -1)
	{
		GzMatrixMultiplication(render->Xnorm[render->matlevel], res, res);
		GzPopMatrix(render);
	}

	render->matlevel = savedMatLevel; //Recall old matlevel
}

float GzNewVz(float currZ)
{
	return (currZ / (INT_MAX - currZ));
}

void GzXformToPerspective(const GzTextureIndex &texi, float Vz, GzTextureIndex &res)
{
	res[0] = texi[0] / (Vz + 1);
	res[1] = texi[1] / (Vz + 1);
}

void GzXformToAffine(const GzTextureIndex &texi, float Vz, GzTextureIndex &res)
{
	res[0] = texi[0] * (Vz + 1);
	res[1] = texi[1] * (Vz + 1);
}

void LoadCubeMaps(GzRender *render)
{
	unsigned char pixel[3];
	unsigned char dummy;
	char foo[8];
	FILE *fd;
	GzColor *image;

	for(int i = 0; i < 6; ++i)
	{
		string filePath = "CubeMapping\\";
		switch(i)
		{
		case 0:
			{
				image = render->cmap.posX;
				filePath += "posx.ppm";
				break;
			}
		case 1:
			{
				image = render->cmap.negX;
				filePath += "negx.ppm";
				break;
			}
		case 2:
			{
				image = render->cmap.posY;
				filePath += "posy.ppm";
				break;
			}
		case 3:
			{
				image = render->cmap.negY;
				filePath += "negy.ppm";
				break;
			}
		case 4:
			{
				image = render->cmap.posZ;
				filePath += "posz.ppm";
				break;
			}
		case 5:
			{
				image = render->cmap.negZ;
				filePath += "negz.ppm";
				break;
			}
		}
		fd = fopen(filePath.c_str(), "rb");
		if(fd == NULL)
		{
			fprintf(stderr, "texture file not found\n");
			exit(-1);
		}
		fscanf(fd, "%s %d %d %c", foo, &render->cmap.xSize, &render->cmap.ySize, &dummy);
		image = (GzColor*)malloc(sizeof(GzColor) * (render->cmap.xSize + 1) * (render->cmap.ySize + 1));
		if(image == NULL)
		{
			fprintf(stderr, "malloc for texture image failed\n");
			exit(-1);
		}

		for(int j = 0; j < render->cmap.xSize * render->cmap.ySize; ++j) /* create array of GzColor values */
		{
			fread(pixel, sizeof(pixel), 1, fd);
			image[j][RED] = (float)((int)pixel[RED]) * (1.0 / 255.0);
			image[j][GREEN] = (float)((int)pixel[GREEN]) * (1.0 / 255.0);
			image[j][BLUE] = (float)((int)pixel[BLUE]) * (1.0 / 255.0);
		}

		fclose(fd);
	}
}

void GetCubeMapColor(GzRender *render, const GzCoord &normal, GzColor &color)
{
	if((fabsf(normal[0]) >= fabsf(normal[1])) && (fabsf(normal[0]) >= fabsf(normal[2])))
	{
		if(normal[0] > 0.0f)
		{
			GetCubeMapTexture(render, CUBEMAPSIDE::RIGHT, 1.0f - (normal[2] / normal[0] + 1.0f) * 0.5f, (normal[1] / normal[0] + 1.0f) * 0.5f, color);
		}
		else if(normal[0] < 0.0f)
		{
			GetCubeMapTexture(render, CUBEMAPSIDE::LEFT, 1.0f - (normal[2] / normal[0] + 1.0f) * 0.5f, 1.0f - (normal[1] / normal[0] + 1.0f) * 0.5f, color);
		}
	}
	else if ((fabsf(normal[1]) >= fabsf(normal[0])) && (fabsf(normal[1]) >= fabsf(normal[2])))
	{
		if(normal[1] > 0.0f)
		{
			GetCubeMapTexture(render, CUBEMAPSIDE::UP, (normal[0] / normal[1] + 1.0f) * 0.5f, 1.0f - (normal[2]/ normal[1] + 1.0f) * 0.5f, color);
		}
		else if(normal[1] < 0.0f)
		{
			GetCubeMapTexture(render, CUBEMAPSIDE::DOWN, 1.0f - (normal[0] / normal[1] + 1.0f) * 0.5f, (normal[2]/normal[1] + 1.0f) * 0.5f, color);
		}
	}
	else if((fabsf(normal[2]) >= fabsf(normal[0])) && (fabsf(normal[2]) >= fabsf(normal[1])))
	{
		if(normal[2] > 0.0f)
		{
			GetCubeMapTexture(render, CUBEMAPSIDE::FRONT, (normal[0] / normal[2] + 1.0f) * 0.5f, (normal[1]/normal[2] + 1.0f) * 0.5f, color);
		}
		else if(normal[2] < 0.0f)
		{
			GetCubeMapTexture(render, CUBEMAPSIDE::BACK, (normal[0] / normal[2] + 1.0f) * 0.5f, 1.0f - (normal[1] /normal[2]+1) * 0.5f, color);
		}
	}
}

void GetCubeMapTexture(GzRender *render, CUBEMAPSIDE cmEnum, float u, float v, GzColor &color)
{
	GzColor *cmPtr;
	switch(cmEnum)
	{
	case CUBEMAPSIDE::LEFT:
		{
			cmPtr = render->cmap.negX;
			break;
		}
	case CUBEMAPSIDE::RIGHT:
		{
			cmPtr = render->cmap.posX;
			break;
		}
	case CUBEMAPSIDE::UP:
		{
			cmPtr = render->cmap.posY;
			break;
		}
	case CUBEMAPSIDE::DOWN:
		{
			cmPtr = render->cmap.negY;
			break;
		}
	case CUBEMAPSIDE::FRONT:
		{
			cmPtr = render->cmap.negZ;
			break;
		}
	case CUBEMAPSIDE::BACK:
		{
			cmPtr = render->cmap.posZ;
			break;
		}
	}

	if(u < 0.0f) u = 0.0f;
	if(v < 0.0f) v = 0.0f;
	if(u > 1.0f) u = 1.0f;
	if(v > 1.0f) v = 1.0f;

	int u0 = u * (render->cmap.xSize - 1);
	int v0 = v * (render->cmap.ySize - 1);

	int u1, v1;

	if(u0 == render->cmap.xSize - 1) u1 = u0 - 1;
	else u1 = u0 + 1;
	if(v0 == render->cmap.ySize - 1) v1 = v0 - 1;
	else v1 = v0 + 1;

	float s = (u * (render->cmap.xSize - 1)) - float(u0);
	float t = (v * (render->cmap.ySize - 1)) - float(v0);

	GzColor *topLeft, *topRight, *bottomRight, *bottomLeft;
	topLeft = &cmPtr[v0 * render->cmap.xSize + u0];
	topRight = &cmPtr[v0 * render->cmap.xSize + u1];
	bottomRight = &cmPtr[v1 * render->cmap.xSize + u1];
	bottomLeft = &cmPtr[v1 * render->cmap.xSize + u0];

	float topLeftRatio = (1 - s) * (1 - t);
	float topRightRatio = s * (1 - t);
	float bottomRightRatio = s * t;
	float bottomLeftRatio = (1 - s) * t;

	color[0] = bottomRightRatio * (*bottomRight)[0] + bottomLeftRatio * (*bottomLeft)[0] + topRightRatio * (*topRight)[0] + topLeftRatio * (*topLeft)[0];
	color[1] = bottomRightRatio * (*bottomRight)[1] + bottomLeftRatio * (*bottomLeft)[1] + topRightRatio * (*topRight)[1] + topLeftRatio * (*topLeft)[1];
	color[2] = bottomRightRatio * (*bottomRight)[2] + bottomLeftRatio * (*bottomLeft)[2] + topRightRatio * (*topRight)[2] + topLeftRatio * (*topLeft)[2];

	/*
	u = fabsf(u);
	v = fabsf(v);
	int umin = int(render->cmap.xSize * u);
	int vmin = int(render->cmap.ySize * v);
	int umax = int(render->cmap.xSize * u) + 1;
	int vmax = int(render->cmap.ySize * v) + 1;
	float ucoef = fabsf(render->cmap.xSize * u - umin);
	float vcoef = fabsf(render->cmap.ySize * v - vmin);

	umin = min(max(umin, 0), render->cmap.xSize - 1);
	umax = min(max(umax, 0), render->cmap.xSize - 1);
	vmin = min(max(vmin, 0), render->cmap.ySize - 1);
	vmax = min(max(vmax, 0), render->cmap.ySize - 1);

	// What follows is a bilinear interpolation
	// along two coordinates u and v.

	color output = (1.0f - vcoef) * ((1.0f - ucoef) * tab[umin  + sizeU * vmin] + ucoef * tab[umax + sizeU * vmin]) +   vcoef * ((1.0f - ucoef) * tab[umin  + sizeU * vmax] + ucoef * tab[umax + sizeU * vmax]);
	*/
}