#pragma once

#include <stddef.h>

#include <NvFlex.h>
#include <vector>

#include <core/types.h>
#include <core/maths.h>


#include <core/platform.h>
#include <core/mesh.h>
#include <core/voxelize.h>
#include <core/sdf.h>

#include <core/perlin.h>


#include <NvFlexExt.h>
#include <NvFlexDevice.h>


using namespace std;

inline float sqr(float x) { return x*x; }

struct Emitter
{
	Emitter() : mSpeed(0.0f), mEnabled(false), mLeftOver(0.0f), mWidth(8)   {}

	Vec3 mPos;
	Vec3 mDir;
	Vec3 mRight;
	float mSpeed;
	bool mEnabled;
	float mLeftOver;
	int mWidth;
};

struct RectEmitter
{
	RectEmitter() : mSpeed(0.0f), mEnabled(false), mLeftOver(0.0f), mDisc(0), mNoise(0)  {}

	Point3 mPos;
	Point3 mSize;
	Vec3 mRot;

	float mSpeed;
	bool mEnabled;
	float mLeftOver;

	int mDisc;
	int mNoise;
	float mNoiseThreshold;
	float mNoiseFreq;
	float mNoiseOffset;

};

struct VolumeBox
{
	VolumeBox()  {}

	Point3 mPos;
	Point3 mSize;
	Vec3 mRot;


};

struct SimBuffers
{
	NvFlexVector<Vec4> positions;
	NvFlexVector<Vec4> restPositions;
	NvFlexVector<Vec3> velocities;
	NvFlexVector<int> phases;
	NvFlexVector<float> densities;
	NvFlexVector<Vec4> anisotropy1;
	NvFlexVector<Vec4> anisotropy2;
	NvFlexVector<Vec4> anisotropy3;
	NvFlexVector<Vec4> normals;
	NvFlexVector<Vec4> smoothPositions;
	NvFlexVector<Vec4> diffusePositions;
	NvFlexVector<Vec4> diffuseVelocities;
	NvFlexVector<int> diffuseIndices;
	NvFlexVector<int> activeIndices;

	// convexes
	NvFlexVector<NvFlexCollisionGeometry> shapeGeometry;
	NvFlexVector<Vec4> shapePositions;
	NvFlexVector<Quat> shapeRotations;
	NvFlexVector<Vec4> shapePrevPositions;
	NvFlexVector<Quat> shapePrevRotations;
	NvFlexVector<int> shapeFlags;

	// rigids
	NvFlexVector<int> rigidOffsets;
	NvFlexVector<int> rigidIndices;
	NvFlexVector<int> rigidMeshSize;
	NvFlexVector<float> rigidCoefficients;
	NvFlexVector<Quat> rigidRotations;
	NvFlexVector<Vec3> rigidTranslations;
	NvFlexVector<Vec3> rigidLocalPositions;
	NvFlexVector<Vec4> rigidLocalNormals;

	// inflatables
	NvFlexVector<int> inflatableTriOffsets;
	NvFlexVector<int> inflatableTriCounts;
	NvFlexVector<float> inflatableVolumes;
	NvFlexVector<float> inflatableCoefficients;
	NvFlexVector<float> inflatablePressures;

	// springs
	NvFlexVector<int> springIndices;
	NvFlexVector<float> springLengths;
	NvFlexVector<float> springStiffness;

	NvFlexVector<int> triangles;
	NvFlexVector<Vec3> triangleNormals;
	NvFlexVector<Vec3> uvs;


	SimBuffers(NvFlexLibrary* l);
	~SimBuffers();

	void MapBuffers();
	void UnmapBuffers();
	void InitBuffers();

};



class FlexSystem {

	public:

	FlexSystem();
	~FlexSystem();

	void initSystem();
	void initParams();

	void getSimTimers();

	void GetParticleBounds(Vec3& lower, Vec3& upper);
	void CreateParticleGrid(Vec3 lower, int dimx, int dimy, int dimz, float radius, Vec3 velocity, float invMass, bool rigid, float rigidStiffness, int phase, float jitter);
	void CreateCenteredParticleGrid(Point3 position, Vec3 rotation, Point3 size, float restDistance, Vec3 velocity, float invMass, bool rigid, int phase, float jitter=0.005f);
	void CreateParticleShape(const Mesh *srcMesh, Vec3 lower, Vec3 scale, float rotation, float spacing, Vec3 velocity, float invMass, bool rigid, float rigidStiffness, int phase, bool skin, float jitter = 0.005f, Vec3 skinOffset = 0.0f, float skinExpand = 0.0f, Vec4 color = Vec4(0.0f), float springStiffness = 0.0f);

			void ClearShapes();
	void setShapes();

	void AddSphere(float radius, Vec3 position, Quat rotation);

	void AddBox(Vec3 halfEdge = Vec3(2.0f), Vec3 center = Vec3(0.0f), Quat quat = Quat(), bool dynamic = false);

	NvFlexTriangleMeshId CreateTriangleMesh(Mesh* m);
	void UpdateTriangleMesh(Mesh* m, NvFlexTriangleMeshId flexMeshId);
	void AddTriangleMesh(NvFlexTriangleMeshId mesh, Vec3 translation, Quat rotation, Vec3 scale);


	void emission();
	void update();

	void initScene();
	void postInitScene();

	void CreateSpring(int i, int j, float stiffness, float give = 0.0f);
	float SampleSDF(const float *sdf, int dim, int x, int y, int z);
	Vec3 SampleSDFGrad(const float *sdf, int dim, int x, int y, int z);

	NvFlexSolver *g_flex;
	NvFlexLibrary* g_flexLib;

	NvFlexExtForceField g_forcefield;
	NvFlexExtForceFieldCallback *g_forcefieldCallback;

	NvFlexParams g_params;
	NvFlexParams g_defaultParams;

	SimBuffers* g_buffers;

	vector<int> g_inactiveIndices;

	bool g_profile;
	bool g_useParticleShape;

	int activeParticles;
	int maxParticles;
	int numDiffuse;


	vector<RectEmitter> g_rectEmitters;
	vector<VolumeBox> g_volumeBoxes;


	int g_numSubsteps;
	float g_dt;	// the time delta used for simulation

	Vec3 g_sceneLower;
	Vec3 g_sceneUpper;

	Mesh *g_mesh;


	int g_maxDiffuseParticles;
	unsigned char g_maxNeighborsPerParticle;


	int nEmitter;
	int nVolumeBoxes;

	double time1;
	double time2;
	double time3;
	double time4;

	NvFlexTimers g_timers;

	float simLatency;

	int cursor;

	const char *g_meshPath;
	Vec3 g_shapePos;
	double g_shapeRot;
	double g_shapeSpacing;
	double g_shapeScale;
	double g_stiffness;

	vector<int> g_meshSkinIndices;
	vector<float> g_meshSkinWeights;
};
