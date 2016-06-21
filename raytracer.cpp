#include <iostream>
#include <cstdio>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"
#include "limits"
#include <math.h>       
#include <algorithm> //max


using namespace std;
using glm::vec3;
using glm::mat3;

//using SDL_image::IMG_Load;

/* ----------------------------------------------------------------------------*/
/* STRUCTS 																																		*/
struct Intersection
{
	vec3 position;
	float distance;
	int triangleIndex;
	float u;
	float v;
};

struct SurfacePoint
{
	vec3 position;
	vec3 color;
	int triangleIndex;
};

struct Surface
{
	vec3 v0;
	vec3 v1;
	vec3 v2;
	vec3 normal;
	vec3 color;
	float reflection;
	float refraction;
	bool bump;
};

struct Light
{
	vec3 position;
	vec3 color;
};

/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */

const int SCREEN_WIDTH = 350;
const int SCREEN_HEIGHT = 350;
const int MAX_DEPTH = 3;
SDL_Surface* screen;

int t;
float epsilon = 0.00001;
float focalLength = SCREEN_HEIGHT/1.1;

// Camera
vec3 cameraPos( 0, 0, -2.8 );
mat3 R;
float yaw;

vector<Triangle> triangles;
vector<Surface> surfaces;
vector<Light> lights;
SurfacePoint bounce[SCREEN_WIDTH][SCREEN_HEIGHT];
vec3 preFXAA[SCREEN_WIDTH][SCREEN_HEIGHT];
float dx[SCREEN_WIDTH][SCREEN_HEIGHT];
float dy[SCREEN_WIDTH][SCREEN_HEIGHT];
vec3 bg(0.2,0.2,0.9);
float nMapScale = 0.01;
//bitmap_image textureImage("texture.bmp"); // doesnt work
//vec3 texture[SCREEN_WIDTH][SCREEN_HEIGHT];


/* Illumination variables */
//vec3 lightPos( 0, -0.5, -0.7 );
//vec3 lightColor = 14.f * vec3( 1, 1, 1 );
vec3 DirectLight( const Intersection& i,  vec3 lightPos, vec3 lightColor  );
vec3 indirectLight = 0.05f*vec3( 1, 1, 1 );

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */
void Update();
void Draw();
bool ClosestIntersection(vec3 start, vec3 dir, const vector<Surface>& surfaces, Intersection& closestIntersection );
void prepareSurfaces(const vector<Triangle>& triangles, const vector<Surface>& surfaces);
vec3 findColor(vec3 start, vec3 dir, int depth);
//vec3 getReflection(Surface surface, vec3 position, vec3 direction);

void prepareSurfaces(const vector<Triangle>& triangles, vector<Surface>& surfaces) {
	//cout << triangles.size();
	// Copy triangles into surfaces
	for (unsigned int i = 0; i < triangles.size(); i++){
		surfaces[i].v0 = triangles[i].v0;
		surfaces[i].v1 = triangles[i].v1;
		surfaces[i].v2 = triangles[i].v2;
		surfaces[i].normal = triangles[i].normal;
		surfaces[i].color = triangles[i].color;
		surfaces[i].reflection = 0;
		surfaces[i].refraction = 0;
	}
	
	// Make some changes to certain triangles
	for (unsigned int i = 0; i < 2; i++){
		surfaces[i].reflection = 0.6;
	}
	for (unsigned int i = 2; i < 10; i++){
		surfaces[i].reflection = 0.3;
	}
	for (unsigned int i = 10; i < 30; i++){
		surfaces[i].reflection = 0.1;
		surfaces[i].refraction = 1.2;
	}
}

void prepareLights(vector<Light>& lights) {
	lights.resize(2);

	// set properties
	lights[0].position = vec3( 0.2, -0.5, -0.7 );
	lights[0].color = 9.f * vec3( 1, 1, 1 );

	lights[1].position = vec3( -0.2, -0.4, -0.6 );
	lights[1].color = 7.f * vec3( 0.8, 1, 0.7);
}

float getNormalMap(float u, float v){
	return (sin(u*M_PI/180)+cos(v*M_PI/180));
}

int main( int argc, char* argv[] )
{
	
	//SDL_Surface *texture = SDL_image::IMG_Load ( "texture.jpeg" );
	//CImg<unsigned char> image("texture.jpeg");

	LoadTestModel( triangles );
	screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
	t = SDL_GetTicks();	// Set start value for timer.


	surfaces.resize(triangles.size());
	prepareSurfaces(triangles, surfaces);
	prepareLights(lights);
	//R[2] = vec3(25,14,12);

	while( NoQuitMessageSDL() )
	{
		Update();
		Draw();
	}

	SDL_SaveBMP( screen, "screenshot.bmp" );
	return 0;
}

void Update()
{
	// Compute frame time:
	int t2 = SDL_GetTicks();
	float dt = float(t2-t);
	t = t2;
	cout << "Render time: " << dt << " ms." << endl;

	Uint8* keystate = SDL_GetKeyState( 0 );

	// Rotate camera
	mat3 rotate;
	
	yaw += 0.1;
	rotate[0][0] = cos(yaw);
	rotate[1][0] = 0;
	rotate[2][0] = sin(yaw);

	rotate[0][2] = -sin(yaw);
	rotate[1][2] = 0;
	rotate[2][2] = cos(yaw);
	//R = R*rotate;

	if( keystate[SDLK_UP] )
	{
	// Move camera forward
		cameraPos.y += 0.1;
	}
	if( keystate[SDLK_DOWN] )
	{
	// Move camera backward
		cameraPos.y -= 0.1;
	}
	if( keystate[SDLK_LEFT] )
	{
	// Move camera to the left
		cameraPos.x += 0.1;
	}
	if( keystate[SDLK_RIGHT] )
	{
	// Move camera to the right
		cameraPos.x -= 0.1;
	}
	/*
	if( keystate[SDLK_w] )
		lightPos.z += 0.1;

	if( keystate[SDLK_s] )
		lightPos.z -= 0.1;
	
	if( keystate[SDLK_a] )
		lightPos.x += 0.1;

	if( keystate[SDLK_d] )
		lightPos.x -= 0.1;

	if( keystate[SDLK_q] )
		lightPos.y += 0.1;

	if( keystate[SDLK_e] )
		lightPos.y -= 0.1;
	*/

}
/*
vec3 getRefractionRay(float nr, vec3 N, vec3 V) {
	//cout << dot(N,V) << "\n";
	//cout << dot(N,V)/(N.length()*V.length()) << "\n";
	float thetai = acos(-dot(N,V)/(N.length()*V.length()));
	float thetat = asin(nr*sin(thetai));	// Snell's law
	vec3 M = (N*cos(thetai)-V)*(1/sin(thetai));
	vec3 T = (sin(thetat)*M) - (cos(thetat)*N);

	//cout << "theta_i=" << thetai << ", theta_t=" << thetat << "\n";

	return T;

}*/

vec3 getRefractionRay(float nr, vec3 N, vec3 V) {
	//cout << dot(N,V) << "\n";
	//cout << dot(N,V)/(N.length()*V.length()) << "\n";

	
	float I = -dot( N, V );
	float n = 1/nr;
	float cosTheta = 1.0f - n * n * (1.0f - I * I);
	vec3 T = (n * V) + (n * I - sqrt(cosTheta)) * N;
	//cout << "T=" << T.x << ", " << T.y << "\n";
	//cout << "theta_i=" << thetai << ", theta_t=" << thetat << "\n";

	return T;

}

vec3 findColor(vec3 start, vec3 dir, int depth) {
	vec3 color = bg;
	Intersection intersection;
	bool intersects = ClosestIntersection(start, normalize(dir), surfaces, intersection );

	if (intersects) {
		int i = intersection.triangleIndex;
		color = vec3(0,0,0);
		vec3 N = surfaces[i].normal;

		// Lighting
		for (unsigned int n = 0; n < lights.size(); n++) {
			vec3 D = DirectLight(intersection, lights[n].position, lights[n].color) + indirectLight;
			color += D * surfaces[i].color;
		}
		if (depth < MAX_DEPTH) {
			// Reflection
			if (surfaces[i].reflection > 0) {
				float dotProduct = -2*dot(dir,N);
				vec3 reflectionDir = (dotProduct*N) + dir;
				color += surfaces[i].reflection*findColor(intersection.position, reflectionDir, depth+1);
			}

			// Refraction
			if (surfaces[i].reflection > 0) {
				vec3 refractionDir = getRefractionRay(1.1, normalize(N), normalize(dir));
				color += surfaces[i].refraction*findColor(intersection.position, refractionDir, depth+1);
			}
		}
	}
	
	return color;
}

void DrawPixel(int x, int y){
	//vec3 color = bg;
	Intersection intersection;
	vec3 dir(x-(SCREEN_WIDTH/2), y-(SCREEN_HEIGHT/2), focalLength);
	vec3 color = findColor(cameraPos,normalize(dir),0);
	preFXAA[x][y] = color;
	//PutPixelSDL( screen, x, y, color );
	//SDL_UpdateRect( screen, 0, 0, 0, 0 );	// Just so we can quickly see whats happening
}

/*
void DrawSSAA(int x, int y){
	// Anti aliasing
	vec3 color[4];
	Intersection intersection;

	for (int ax = 0; ax < 2; ax++) {
		for (int ay = 0; ay < 2; ay++) {
			color[ax+(2*ay)] = bg;
			vec3 dir(x+ax-(SCREEN_WIDTH/2), y+ay-(SCREEN_HEIGHT/2), focalLength);
						bool intersects = ClosestIntersection(cameraPos, dir, triangles, intersection );
			if (intersects) {
				int i = intersection.triangleIndex;
				vec3 D = DirectLight(intersection);
				color[ax+(2*ay)] = triangles[i].color * (D + indirectLight);
			}
		}
	}
	
	vec3 aColor(0,0,0);
	//vec3 aColor = color[4];
	for (int i = 0; i < 4; i++) {
		aColor += color[i];			
	}
	aColor /= 4;
	PutPixelSDL( screen, x, y, aColor);
}
*/

void Draw()
{
	if( SDL_MUSTLOCK(screen) )
		SDL_LockSurface(screen);

	//DrawGI(0,0,1);	
	// /*
	for( int y=0; y<SCREEN_HEIGHT; ++y )
	{
		for( int x=0; x<SCREEN_WIDTH; ++x )
		{
			DrawPixel(x,y);
					
		}
	}
	// */
	
	// FXAA
	// Find dx and dy images
	for( int y=0; y<SCREEN_HEIGHT-1; ++y )
	{
		for( int x=0; x<SCREEN_WIDTH-1; ++x )
		{
			// Storing this and next pixel's colors
			float cr = preFXAA[x][y].x;
			float cxr = preFXAA[x+1][y].x;
			float cyr = preFXAA[x][y+1].x;
			float cg = preFXAA[x][y].y;
			float cxg = preFXAA[x+1][y].y;
			float cyg = preFXAA[x][y+1].y;
			float cb = preFXAA[x][y].z;
			float cxb = preFXAA[x+1][y].z;
			float cyb = preFXAA[x][y+1].z;
			// Find average differences (for each color)
			dx[x][y] = ((cxr-cr)+(cxg-cg)+(cxb-cb))/3.f;
			dx[x+1][y+1] = dx[x][y];
			dy[x][y] = ((cyr-cr)+(cyg-cg)+(cyb-cb))/3.f;
			dy[x+1][y+1] = dy[x][y];
			//cout << dx[x][y] << "\n";

		}
	}

	
	
	// Threshold
	float threshold = 0.0005;
	for( int y=0; y<SCREEN_HEIGHT; ++y )
	{
		for( int x=0; x<SCREEN_WIDTH; ++x )
		{
			if (dx[x][y] > threshold || dx[x][y] < -threshold) {
				dx[x][y] = 1;
			}
			else dx[x][y] = 0;
			if (dy[x][y] > threshold || dy[x][y] < -threshold) {
				dy[x][y] = 1;
			}
			else dy[x][y] = 0;
			//dx[x][y] = 1;
			//dy[x][y] = 1;
		}
	}

	// Blur according to dx and dy then draw
	for( int y=0; y<SCREEN_HEIGHT-1; ++y )
	{
		for( int x=0; x<SCREEN_WIDTH-1; ++x )
		{
			vec3 color = preFXAA[x][y];
			
			if (dx[x][y] > 0 && dy[x][y] > 0) {
				color = (preFXAA[x][y]+preFXAA[x+1][y]+preFXAA[x][y+1])/3.f;
			}
			else if (dx[x][y] > 0) {
				color = (preFXAA[x][y]+preFXAA[x+1][y])/2.f;
			}
			else if (dy[x][y] > 0) {
				color = (preFXAA[x][y]+preFXAA[x][y+1])/2.f;
			}

			PutPixelSDL(screen, x, y, color);
			PutPixelSDL(screen, x+1, y+1, color);
		}
	}
	if( SDL_MUSTLOCK(screen) )
		SDL_UnlockSurface(screen);

	SDL_UpdateRect( screen, 0, 0, 0, 0 );

}


bool ClosestIntersection(vec3 start, vec3 dir, const vector<Surface>& surfaces, Intersection& closestIntersection ){

	bool intersects = false;
	float m = std::numeric_limits<float>::max();
	closestIntersection.distance = m;

	for (unsigned int i = 0; i < surfaces.size(); i++){

		Surface surface = surfaces[i];

		vec3 v0 = surface.v0;
		vec3 v1 = surface.v1;
		vec3 v2 = surface.v2;
		vec3 e1 = v1 - v0;
		vec3 e2 = v2 - v0;
		vec3 b = start - v0;
		// /*
		mat3 A( -dir, e1, e2 );
		vec3 x = glm::inverse( A ) * b;
		float t = x.x, u = x.y, v = x.z;
		// */
		/*
		// Closed form
		vec3 p = cross(dir, e2);
		vec3 q = cross(start, e1);
		float t = (q.x*e2.x)+(q.y*e2.y)+(q.z*e2.z);
		//t = -t;
		float u = (p.x*start.x)+(p.y*start.y)+(p.z*start.z);
		float v = (q.x*dir.x)+(q.y*dir.y)+(q.z*dir.z);
		cout << t << ", " << u << ", " << v << "| ";
		// */
		
		if ( t > epsilon && u >= 0 && v >= 0 && (u+v) <= 1) {
			intersects = true;
			if (t < closestIntersection.distance){
				closestIntersection.position = start + (t*dir);
				closestIntersection.distance =  t;
				closestIntersection.triangleIndex = i;
				closestIntersection.u = u;
				closestIntersection.v = v;
			}
		}
	//cout <<"closest intersection co-ord: (" << closestIntersection.position.x << "," << closestIntersection.position.y << ")"; 
	}
	return intersects;
}

vec3 DirectLight( const Intersection& i, vec3 lightPos, vec3 lightColor ) {

	// Work out distance (light -> point)
	vec3 rhat = lightPos - i.position;
	float r = pow(pow(rhat.x,2) + pow(rhat.y,2) + pow(rhat.z,2),0.5);

	// Test for shadow
	Intersection shadowIntersect;
	ClosestIntersection(i.position, normalize(lightPos - i.position), surfaces, shadowIntersect);
	//cout << "S: " << shadowIntersect.distance << ", ";
	if (shadowIntersect.distance < r) {// && shadowIntersect.distance > 0.0001) {
		return vec3(0,0,0);
	}

	// Some bump mapping using sin/cos - not really..
	vec3 nhat = surfaces[i.triangleIndex].normal;
	//nhat.x += 1+nMapScale*sin(i.u*4*M_PI);
	//nhat.y += 1+nMapScale*cos(i.v*4*M_PI);
	float dotProduct = dot(rhat,nhat);
	//cout << "R: " << r << ", ";
	float A = 4 * M_PI * pow(r,2);
	vec3 B = lightColor/A;
	vec3 D = B * max(dotProduct,float(0));

	return D;
}
