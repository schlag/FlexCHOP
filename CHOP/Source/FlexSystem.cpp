#include "FlexSystem.h"


void ErrorCallback(NvFlexErrorSeverity, const char* msg, const char* file, int line)
{
	printf("Flex: %s - %s:%d\n", msg, file, line);
	//g_Error = true;
	//assert(0); asserts are bad for TeamCity
}


SimBuffers::SimBuffers(NvFlexLibrary* l) :
	positions(l), restPositions(l), velocities(l), phases(l), densities(l),
	anisotropy1(l), anisotropy2(l), anisotropy3(l), normals(l), smoothPositions(l),
	diffusePositions(l), diffuseVelocities(l), diffuseIndices(l), activeIndices(l),
	shapeGeometry(l), shapePositions(l), shapeRotations(l), shapePrevPositions(l),
	shapePrevRotations(l), shapeFlags(l), rigidOffsets(l), rigidIndices(l), rigidMeshSize(l),
	rigidCoefficients(l), rigidRotations(l), rigidTranslations(l),
	rigidLocalPositions(l), rigidLocalNormals(l), inflatableTriOffsets(l),
	inflatableTriCounts(l), inflatableVolumes(l), inflatableCoefficients(l),
	inflatablePressures(l), springIndices(l), springLengths(l),
	springStiffness(l), triangles(l), triangleNormals(l), uvs(l){


}

SimBuffers::~SimBuffers() {

	positions.destroy();
	restPositions.destroy();
	velocities.destroy();
	phases.destroy();
	densities.destroy();
	anisotropy1.destroy();
	anisotropy2.destroy();
	anisotropy3.destroy();
	normals.destroy();
	diffusePositions.destroy();
	diffuseVelocities.destroy();
	diffuseIndices.destroy();
	smoothPositions.destroy();
	activeIndices.destroy();

	// convexes
	shapeGeometry.destroy();
	shapePositions.destroy();
	shapeRotations.destroy();
	shapePrevPositions.destroy();
	shapePrevRotations.destroy();
	shapeFlags.destroy();

	// rigids
	rigidOffsets.destroy();
	rigidIndices.destroy();
	rigidMeshSize.destroy();
	rigidCoefficients.destroy();
	rigidRotations.destroy();
	rigidTranslations.destroy();
	rigidLocalPositions.destroy();
	rigidLocalNormals.destroy();

	// springs
	springIndices.destroy();
	springLengths.destroy();
	springStiffness.destroy();

	// inflatables
	inflatableTriOffsets.destroy();
	inflatableTriCounts.destroy();
	inflatableVolumes.destroy();
	inflatableCoefficients.destroy();
	inflatablePressures.destroy();

	// triangles
	triangles.destroy();
	triangleNormals.destroy();
	uvs.destroy();


}

void SimBuffers::MapBuffers() {

	positions.map();
	restPositions.map();
	velocities.map();
	phases.map();
	densities.map();
	anisotropy1.map();
	anisotropy2.map();
	anisotropy3.map();
	normals.map();
	diffusePositions.map();
	diffuseVelocities.map();
	diffuseIndices.map();
	smoothPositions.map();
	activeIndices.map();

	// convexes
	shapeGeometry.map();
	shapePositions.map();
	shapeRotations.map();
	shapePrevPositions.map();
	shapePrevRotations.map();
	shapeFlags.map();

	rigidOffsets.map();
	rigidIndices.map();
	rigidMeshSize.map();
	rigidCoefficients.map();
	rigidRotations.map();
	rigidTranslations.map();
	rigidLocalPositions.map();
	rigidLocalNormals.map();

	springIndices.map();
	springLengths.map();
	springStiffness.map();

	// inflatables
	inflatableTriOffsets.map();
	inflatableTriCounts.map();
	inflatableVolumes.map();
	inflatableCoefficients.map();
	inflatablePressures.map();

	triangles.map();
	triangleNormals.map();
	uvs.map();


}

void SimBuffers::UnmapBuffers() {

	// particles
	positions.unmap();
	restPositions.unmap();
	velocities.unmap();
	phases.unmap();
	densities.unmap();
	anisotropy1.unmap();
	anisotropy2.unmap();
	anisotropy3.unmap();
	normals.unmap();
	diffusePositions.unmap();
	diffuseVelocities.unmap();
	diffuseIndices.unmap();
	smoothPositions.unmap();
	activeIndices.unmap();

	// convexes
	shapeGeometry.unmap();
	shapePositions.unmap();
	shapeRotations.unmap();
	shapePrevPositions.unmap();
	shapePrevRotations.unmap();
	shapeFlags.unmap();

	// rigids
	rigidOffsets.unmap();
	rigidIndices.unmap();
	rigidMeshSize.unmap();
	rigidCoefficients.unmap();
	rigidRotations.unmap();
	rigidTranslations.unmap();
	rigidLocalPositions.unmap();
	rigidLocalNormals.unmap();

	// springs
	springIndices.unmap();
	springLengths.unmap();
	springStiffness.unmap();

	// inflatables
	inflatableTriOffsets.unmap();
	inflatableTriCounts.unmap();
	inflatableVolumes.unmap();
	inflatableCoefficients.unmap();
	inflatablePressures.unmap();

	// triangles
	triangles.unmap();
	triangleNormals.unmap();
	uvs.unmap();

}

void SimBuffers::InitBuffers() {

	positions.resize(0);
	velocities.resize(0);
	phases.resize(0);

	/*rigidOffsets.resize(0);
	rigidIndices.resize(0);
	rigidMeshSize.resize(0);
	rigidRotations.resize(0);
	rigidTranslations.resize(0);
	rigidCoefficients.resize(0);
	rigidLocalPositions.resize(0);
	rigidLocalNormals.resize(0);

	springIndices.resize(0);
	springLengths.resize(0);
	springStiffness.resize(0);
	triangles.resize(0);
	triangleNormals.resize(0);
	uvs.resize(0);*/

	shapeGeometry.resize(0);
	shapePositions.resize(0);
	shapeRotations.resize(0);
	shapePrevPositions.resize(0);
	shapePrevRotations.resize(0);
	shapeFlags.resize(0);

}

FlexSystem::FlexSystem()
{	

	g_profile = false;
	
	nEmitter = 0;
	nVolumeBoxes = 0;

	maxParticles = 9;
	g_maxDiffuseParticles=0;
	numDiffuse = 0;

	g_flex = NULL;


	time1 = 0;
	time2 = 0;
	time3 = 0;
	time4 = 0;

	memset(&g_timers, 0, sizeof(g_timers));

	cursor = 0;

}

FlexSystem::~FlexSystem(){


	if (g_flex)
	{
		if(g_buffers)
			delete g_buffers;

		NvFlexDestroySolver(g_flex);
		NvFlexShutdown(g_flexLib);
	}

	

	NvFlexDeviceDestroyCudaContext();

}

void FlexSystem::getSimTimers() {

	if (g_profile) {
		memset(&g_timers, 0, sizeof(g_timers));
		NvFlexGetTimers(g_flex, &g_timers);
		simLatency = g_timers.total;
	}
	else
		simLatency = NvFlexGetDeviceLatency(g_flex);

}

void FlexSystem::initSystem() {

	int g_device = -1;
	g_device = NvFlexDeviceGetSuggestedOrdinal();

	// Create an optimized CUDA context for Flex and set it on the 
	// calling thread. This is an optional call, it is fine to use 
	// a regular CUDA context, although creating one through this API
	// is recommended for best performance.
	bool success = NvFlexDeviceCreateCudaContext(g_device);

	if (!success)
	{
		printf("Error creating CUDA context.\n");
	}

	NvFlexInitDesc desc;
	desc.deviceIndex = g_device;
	desc.enableExtensions = true;
	desc.renderDevice = 0;
	desc.renderContext = 0;
	desc.computeType = eNvFlexCUDA;

	g_flexLib = NvFlexInit(NV_FLEX_VERSION, ErrorCallback, &desc);
}


void FlexSystem::initParams() {

	// sim params
	g_params.gravity[0] = 0.0f;
	g_params.gravity[1] = -9.8f;
	g_params.gravity[2] = 0.0f;

	g_params.wind[0] = 0.0f;
	g_params.wind[1] = 0.0f;
	g_params.wind[2] = 0.0f;

	g_params.radius = 0.15f;
	g_params.viscosity = 0.0f;
	g_params.dynamicFriction = 0.0f;
	g_params.staticFriction = 0.0f;
	g_params.particleFriction = 0.0f; // scale friction between particles by default
	g_params.freeSurfaceDrag = 0.0f;
	g_params.drag = 0.0f;
	g_params.lift = 0.0f;
	g_params.numIterations = 3;
	g_params.fluidRestDistance = 0.0f;
	g_params.solidRestDistance = 0.0f;

	g_params.anisotropyScale = 1.0f;
	g_params.anisotropyMin = 0.1f;
	g_params.anisotropyMax = 2.0f;
	g_params.smoothing = 1.0f;

	g_params.dissipation = 0.0f;
	g_params.damping = 0.0f;
	g_params.particleCollisionMargin = 0.0f;
	g_params.shapeCollisionMargin = 0.0f;
	g_params.collisionDistance = 0.0f;
	g_params.plasticThreshold = 0.0f;
	g_params.plasticCreep = 0.0f;
	g_params.fluid = true;
	g_params.sleepThreshold = 0.0f;
	g_params.shockPropagation = 0.0f;
	g_params.restitution = 0.001f;

	g_params.maxSpeed = FLT_MAX;
	g_params.maxAcceleration = 100.0f;	// approximately 10x gravity

	g_params.relaxationMode = eNvFlexRelaxationLocal;
	g_params.relaxationFactor = 1.0f;
	g_params.solidPressure = 1.0f;
	g_params.adhesion = 0.0f;
	g_params.cohesion = 0.1f;
	g_params.surfaceTension = 0.0f;
	g_params.vorticityConfinement = 80.0f;
	g_params.buoyancy = 1.0f;
	g_params.diffuseThreshold = 100.0f;
	g_params.diffuseBuoyancy = 1.0f;
	g_params.diffuseDrag = 0.8f;
	g_params.diffuseBallistic = 16;
	g_params.diffuseSortAxis[0] = 0.0f;
	g_params.diffuseSortAxis[1] = 0.0f;
	g_params.diffuseSortAxis[2] = 0.0f;
	g_params.diffuseLifetime = 2.0f;

	g_params.numPlanes = 0;

	g_forcefield.mPosition[0] = 0.0f;
	g_forcefield.mPosition[1] = 0.0f;
	g_forcefield.mPosition[2] = 0.0f;
	g_forcefield.mRadius = 1.0f;
	g_forcefield.mStrength = -30.0f;
	g_forcefield.mLinearFalloff = true;
	g_forcefield.mMode = eNvFlexExtModeForce;

	g_shapeScale = 1.0;
	g_shapeSpacing = 0.5;
	g_stiffness = 0.0;
}

void FlexSystem::initScene(){


	RandInit();

	cursor = 0;


	if (g_flex)
	{
		
		if (g_buffers)
			delete g_buffers;

		NvFlexDestroySolver(g_flex);
		g_flex = NULL;

	}

	// alloc buffers
	g_buffers = new SimBuffers(g_flexLib);

	// map during initialization
	g_buffers->MapBuffers();
	g_buffers->InitBuffers();

	g_rectEmitters.resize(0);
	g_volumeBoxes.resize(0);

	initParams();

	g_numSubsteps = 2;

	g_sceneLower = FLT_MAX;
	g_sceneUpper = -FLT_MAX;

	ClearShapes();


	maxParticles = 64;

	g_maxDiffuseParticles = 0;
	g_maxNeighborsPerParticle = 96;

}

void FlexSystem::postInitScene(){


	if (!g_params.fluid) {
		g_params.particleCollisionMargin = g_params.radius*0.5f;
		g_params.shapeCollisionMargin = g_params.radius*0.5f;
	}

	

	if(g_params.fluid)
		g_params.fluidRestDistance = g_params.radius*0.65f;

	// by default solid particles use the maximum radius
	if (g_params.fluid && g_params.solidRestDistance == 0.0f)
		g_params.solidRestDistance = g_params.fluidRestDistance;
	else
		g_params.solidRestDistance = g_params.radius;

	// collision distance with shapes half the radius
	if (g_params.collisionDistance == 0.0f)
	{
		g_params.collisionDistance = g_params.radius*0.5f;

		if (g_params.fluid)
			g_params.collisionDistance = g_params.fluidRestDistance*0.5f;
	}

	// default particle friction to 10% of shape friction
	if (g_params.particleFriction == 0.0f)
		g_params.particleFriction = g_params.dynamicFriction*0.1f;

	// add a margin for detecting contacts between particles and shapes
	if (g_params.shapeCollisionMargin == 0.0f)
		g_params.shapeCollisionMargin = g_params.collisionDistance*0.5f;


	g_maxDiffuseParticles = 0;

	if (g_useParticleShape) {
		if (g_mesh) {
			delete g_mesh;
		}
		if (g_meshPath && strlen(g_meshPath) > 0) {
			g_mesh = ImportMesh(g_meshPath);
			CreateParticleShape(g_mesh, g_shapePos, g_shapeScale, 0.0f, g_params.radius * g_shapeSpacing, Vec3(0.0f, 0.0f, 0.0f), 1.0f, true, 1.f, NvFlexMakePhase(0, eNvFlexPhaseSelfCollide), false, 0.0f);
		}
		
	} else {
		if (g_params.fluid) {

			for (int i = 0; i < nVolumeBoxes; i++) {

				CreateCenteredParticleGrid(Point3(g_volumeBoxes[i].mPos.x, g_volumeBoxes[i].mPos.y, g_volumeBoxes[i].mPos.z), g_volumeBoxes[i].mRot, Point3(g_volumeBoxes[i].mSize.x, g_volumeBoxes[i].mSize.y, g_volumeBoxes[i].mSize.z), g_params.fluidRestDistance, Vec3(0.0f), 1, false, NvFlexMakePhase(0, eNvFlexPhaseSelfCollide | eNvFlexPhaseFluid), g_params.fluidRestDistance*0.01f);
			}
		}
		else {

			for (int i = 0; i < nVolumeBoxes; i++) {

				CreateCenteredParticleGrid(Point3(g_volumeBoxes[i].mPos.x, g_volumeBoxes[i].mPos.y, g_volumeBoxes[i].mPos.z), g_volumeBoxes[i].mRot, Point3(g_volumeBoxes[i].mSize.x, g_volumeBoxes[i].mSize.y, g_volumeBoxes[i].mSize.z), g_params.radius, Vec3(0.0f), 1, false, NvFlexMakePhase(0, eNvFlexPhaseSelfCollide), g_params.radius*0.01f);
			}

		}
	}

	g_params.anisotropyScale = 3.0f/g_params.radius;


	uint32_t numParticles = g_buffers->positions.size(); //non zero if init volume boxes

	//**** IF PARTICLE GRID

	if (g_buffers->positions.size()) {
		g_buffers->activeIndices.resize(numParticles);
		for (size_t i = 0; i < g_buffers->activeIndices.size(); ++i)
			g_buffers->activeIndices[i] = i;
	}
	g_inactiveIndices.resize(maxParticles - numParticles);

	for (size_t i = 0; i < g_inactiveIndices.size(); ++i)
		g_inactiveIndices[i] = i + numParticles;


	g_buffers->diffusePositions.resize(g_maxDiffuseParticles);
	g_buffers->diffuseVelocities.resize(g_maxDiffuseParticles);
	g_buffers->diffuseIndices.resize(g_maxDiffuseParticles);
	
	// for fluid rendering these are the Laplacian smoothed positions
	g_buffers->smoothPositions.resize(maxParticles);

	g_buffers->normals.resize(0);
	g_buffers->normals.resize(maxParticles);


	// resize particle buffers to fit
	g_buffers->positions.resize(maxParticles);
	g_buffers->velocities.resize(maxParticles);
	g_buffers->phases.resize(maxParticles);

	g_buffers->densities.resize(maxParticles);
	g_buffers->anisotropy1.resize(maxParticles);
	g_buffers->anisotropy2.resize(maxParticles);
	g_buffers->anisotropy3.resize(maxParticles);

	// save rest positions
	g_buffers->restPositions.resize(g_buffers->positions.size());
	for (int i = 0; i < g_buffers->positions.size(); ++i)
		g_buffers->restPositions[i] = g_buffers->positions[i];


	g_flex = NvFlexCreateSolver(g_flexLib, maxParticles, g_maxDiffuseParticles, g_maxNeighborsPerParticle);
	if (g_forcefieldCallback)
		NvFlexExtDestroyForceFieldCallback(g_forcefieldCallback);
	g_forcefieldCallback = NvFlexExtCreateForceFieldCallback(g_flex);
}


void FlexSystem::AddBox(Vec3 halfEdge, Vec3 center, Quat quat, bool dynamic)
{
	// transform
	g_buffers->shapePositions.push_back(Vec4(center.x, center.y, center.z, 0.0f));
	g_buffers->shapeRotations.push_back(quat);

	g_buffers->shapePrevPositions.push_back(g_buffers->shapePositions.back());
	g_buffers->shapePrevRotations.push_back(g_buffers->shapeRotations.back());

	NvFlexCollisionGeometry geo;
	geo.box.halfExtents[0] = halfEdge.x;
	geo.box.halfExtents[1] = halfEdge.y;
	geo.box.halfExtents[2] = halfEdge.z;

	g_buffers->shapeGeometry.push_back(geo);
	g_buffers->shapeFlags.push_back(NvFlexMakeShapeFlags(eNvFlexShapeBox, dynamic));
}

void FlexSystem::AddSphere(float radius, Vec3 position, Quat rotation)
{
	NvFlexCollisionGeometry geo;
	geo.sphere.radius = radius;
	g_buffers->shapeGeometry.push_back(geo);

	g_buffers->shapePositions.push_back(Vec4(position, 0.0f));
	g_buffers->shapeRotations.push_back(rotation);

	g_buffers->shapePrevPositions.push_back(g_buffers->shapePositions.back());
	g_buffers->shapePrevRotations.push_back(g_buffers->shapeRotations.back());

	int flags = NvFlexMakeShapeFlags(eNvFlexShapeSphere, false);
	g_buffers->shapeFlags.push_back(flags);
}

void FlexSystem::GetParticleBounds(Vec3& lower, Vec3& upper)
{
	lower = Vec3(FLT_MAX);
	upper = Vec3(-FLT_MAX);

	for (size_t i=0; i < g_buffers->positions.size(); ++i)
	{
		lower = Min(Vec3(g_buffers->positions[i]), lower);
		upper = Max(Vec3(g_buffers->positions[i]), upper);
	}
}

void FlexSystem::CreateParticleGrid(Vec3 lower, int dimx, int dimy, int dimz, float radius, Vec3 velocity, float invMass, bool rigid, float rigidStiffness, int phase, float jitter=0.005f)
{


	for (int x=0; x < dimx; ++x)
	{
		for (int y=0; y < dimy; ++y)
		{
			for (int z=0; z < dimz; ++z)
			{

				Vec3 position = lower + Vec3(float(x), float(y), float(z))*radius + RandomUnitVector()*jitter;

				g_buffers->positions.push_back(Vec4(position.x, position.y, position.z, invMass));
				g_buffers->velocities.push_back(velocity);
				g_buffers->phases.push_back(phase);
			}
		}
	}


}

void FlexSystem::CreateCenteredParticleGrid(Point3 center, Vec3 rotation, Point3 size, float restDistance, Vec3 velocity, float invMass, bool rigid, int phase, float jitter)
{

	int dx = int(ceilf(size.x / restDistance));
	int dy = int(ceilf(size.y / restDistance));
	int dz = int(ceilf(size.z / restDistance));

	for (int x=0; x < dx; ++x)
	{
		for (int y=0; y < dy; ++y)
		{
			for (int z=0; z < dz; ++z)
			{
				Point3 position = restDistance*Point3(float(x) - 0.5*(dx-1), float(y) - 0.5*(dy-1), float(z) - 0.5*(dz-1)) + RandomUnitVector()*jitter;
				position = TranslationMatrix(center) * TransformMatrix(Rotation(rotation.y, rotation.z, rotation.x), Point3(0.f))*position;

				if(cursor<maxParticles-1) {

					g_buffers->positions.push_back(Vec4(position.x, position.y, position.z, invMass));
					g_buffers->velocities.push_back(velocity);
					g_buffers->phases.push_back(phase);

					cursor++;
				}

				
				
			}
		}
	}
}

void FlexSystem::CreateParticleShape(const Mesh *srcMesh, Vec3 lower, Vec3 scale, float rotation, float spacing, Vec3 velocity, float invMass, bool rigid, float rigidStiffness, int phase, bool skin, float jitter, Vec3 skinOffset, float skinExpand, Vec4 color, float springStiffness)
{
	if (rigid && g_buffers->rigidIndices.empty())
		g_buffers->rigidOffsets.push_back(0);

	if (!srcMesh)
		return;

	// duplicate mesh
	Mesh mesh;
	mesh.AddMesh(*srcMesh);

	int startIndex = int(g_buffers->positions.size());

	{
		mesh.Transform(RotationMatrix(rotation, Vec3(0.0f, 1.0f, 0.0f)));

		Vec3 meshLower, meshUpper;
		mesh.GetBounds(meshLower, meshUpper);

		Vec3 edges = meshUpper - meshLower;
		float maxEdge = max(max(edges.x, edges.y), edges.z);

		// put mesh at the origin and scale to specified size
		Matrix44 xform = ScaleMatrix(scale / maxEdge) * TranslationMatrix(Point3(-meshLower));

		mesh.Transform(xform);
		mesh.GetBounds(meshLower, meshUpper);

		// recompute expanded edges
		edges = meshUpper - meshLower;
		maxEdge = max(max(edges.x, edges.y), edges.z);

		// tweak spacing to avoid edge cases for particles laying on the boundary
		// just covers the case where an edge is a whole multiple of the spacing.
		float spacingEps = spacing * (1.0f - 1e-4f);

		// make sure to have at least one particle in each dimension
		int dx, dy, dz;
		dx = spacing > edges.x ? 1 : int(edges.x / spacingEps);
		dy = spacing > edges.y ? 1 : int(edges.y / spacingEps);
		dz = spacing > edges.z ? 1 : int(edges.z / spacingEps);

		int maxDim = max(max(dx, dy), dz);

		// expand border by two voxels to ensure adequate sampling at edges
		meshLower -= 2.0f * Vec3(spacing);
		meshUpper += 2.0f * Vec3(spacing);
		maxDim += 4;

		vector<uint32_t> voxels(maxDim * maxDim * maxDim);

		// we shift the voxelization bounds so that the voxel centers
		// lie symmetrically to the center of the object. this reduces the
		// chance of missing features, and also better aligns the particles
		// with the mesh
		Vec3 meshOffset;
		meshOffset.x = 0.5f * (spacing - (edges.x - (dx - 1) * spacing));
		meshOffset.y = 0.5f * (spacing - (edges.y - (dy - 1) * spacing));
		meshOffset.z = 0.5f * (spacing - (edges.z - (dz - 1) * spacing));
		meshLower -= meshOffset;

		//Voxelize(*mesh, dx, dy, dz, &voxels[0], meshLower - Vec3(spacing*0.05f) , meshLower + Vec3(maxDim*spacing) + Vec3(spacing*0.05f));
		Voxelize((const float *)&mesh.m_positions[0], mesh.m_positions.size(), (const int *)&mesh.m_indices[0], mesh.m_indices.size(), maxDim, maxDim, maxDim, &voxels[0], meshLower, meshLower + Vec3(maxDim * spacing));

		vector<int> indices(maxDim * maxDim * maxDim);
		vector<float> sdf(maxDim * maxDim * maxDim);
		MakeSDF(&voxels[0], maxDim, maxDim, maxDim, &sdf[0]);

		for (int x = 0; x < maxDim; ++x)
		{
			for (int y = 0; y < maxDim; ++y)
			{
				for (int z = 0; z < maxDim; ++z)
				{
					const int index = z * maxDim * maxDim + y * maxDim + x;

					// if voxel is marked as occupied the add a particle
					if (voxels[index])
					{
						if (rigid)
							g_buffers->rigidIndices.push_back(int(g_buffers->positions.size()));

						Vec3 position = lower + meshLower + spacing * Vec3(float(x) + 0.5f, float(y) + 0.5f, float(z) + 0.5f) + RandomUnitVector() * jitter;

						// normalize the sdf value and transform to world scale
						Vec3 n = SafeNormalize(SampleSDFGrad(&sdf[0], maxDim, x, y, z));
						float d = sdf[index] * maxEdge;

						if (rigid)
							g_buffers->rigidLocalNormals.push_back(Vec4(n, d));

						// track which particles are in which cells
						indices[index] = g_buffers->positions.size();

						g_buffers->positions.push_back(Vec4(position.x, position.y, position.z, invMass));
						g_buffers->velocities.push_back(velocity);
						g_buffers->phases.push_back(phase);
					}
				}
			}
		}
		mesh.Transform(ScaleMatrix(1.0f + skinExpand) * TranslationMatrix(Point3(-0.5f * (meshUpper + meshLower))));
		mesh.Transform(TranslationMatrix(Point3(lower + 0.5f * (meshUpper + meshLower))));

		if (springStiffness > 0.0f)
		{
			// construct cross link springs to occupied cells
			for (int x = 0; x < maxDim; ++x)
			{
				for (int y = 0; y < maxDim; ++y)
				{
					for (int z = 0; z < maxDim; ++z)
					{
						const int centerCell = z * maxDim * maxDim + y * maxDim + x;

						// if voxel is marked as occupied the add a particle
						if (voxels[centerCell])
						{
							const int width = 1;

							// create springs to all the neighbors within the width
							for (int i = x - width; i <= x + width; ++i)
							{
								for (int j = y - width; j <= y + width; ++j)
								{
									for (int k = z - width; k <= z + width; ++k)
									{
										const int neighborCell = k * maxDim * maxDim + j * maxDim + i;

										if (neighborCell > 0 && neighborCell < int(voxels.size()) && voxels[neighborCell] && neighborCell != centerCell)
										{
											CreateSpring(indices[neighborCell], indices[centerCell], springStiffness);
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	if (skin)
	{
		g_buffers->rigidMeshSize.push_back(mesh.GetNumVertices());

		int startVertex = 0;

		if (!g_mesh)
			g_mesh = new Mesh();

		// append to mesh
		startVertex = g_mesh->GetNumVertices();

		g_mesh->Transform(TranslationMatrix(Point3(skinOffset)));
		g_mesh->AddMesh(mesh);

		const Colour colors[7] =
				{
						Colour(0.0f, 0.5f, 1.0f),
						Colour(0.797f, 0.354f, 0.000f),
						Colour(0.000f, 0.349f, 0.173f),
						Colour(0.875f, 0.782f, 0.051f),
						Colour(0.01f, 0.170f, 0.453f),
						Colour(0.673f, 0.111f, 0.000f),
						Colour(0.612f, 0.194f, 0.394f)};

		for (uint32_t i = startVertex; i < g_mesh->GetNumVertices(); ++i)
		{
			int indices[4] = {-1, -1, -1, -1};
			float distances[4] = {FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX};

			if (LengthSq(color) == 0.0f)
				g_mesh->m_colours[i] = 1.25f * colors[phase % 7];
			else
				g_mesh->m_colours[i] = Colour(color);

			// find closest n particles
			for (int j = startIndex; j < g_buffers->positions.size(); ++j)
			{
				float dSq = LengthSq(Vec3(g_mesh->m_positions[i]) - Vec3(g_buffers->positions[j]));

				// insertion sort
				int w = 0;
				for (; w < 4; ++w)
					if (dSq < distances[w])
						break;

				if (w < 4)
				{
					// shuffle down
					for (int s = 3; s > w; --s)
					{
						indices[s] = indices[s - 1];
						distances[s] = distances[s - 1];
					}

					distances[w] = dSq;
					indices[w] = int(j);
				}
			}

			// weight particles according to distance
			float wSum = 0.0f;

			for (int w = 0; w < 4; ++w)
			{
				// convert to inverse distance
				distances[w] = 1.0f / (0.1f + powf(distances[w], .125f));

				wSum += distances[w];
			}

			float weights[4];
			for (int w = 0; w < 4; ++w)
				weights[w] = distances[w] / wSum;

			for (int j = 0; j < 4; ++j)
			{
				g_meshSkinIndices.push_back(indices[j]);
				g_meshSkinWeights.push_back(weights[j]);
			}
		}
	}

	if (rigid)
	{
		g_buffers->rigidCoefficients.push_back(rigidStiffness);
		g_buffers->rigidOffsets.push_back(int(g_buffers->rigidIndices.size()));
	}
}

void FlexSystem::ClearShapes()
{
	g_buffers->shapeGeometry.resize(0);
	g_buffers->shapePositions.resize(0);
	g_buffers->shapeRotations.resize(0);
	g_buffers->shapePrevPositions.resize(0);
	g_buffers->shapePrevRotations.resize(0);
	g_buffers->shapeFlags.resize(0);
}

void FlexSystem::setShapes(){

	if(g_buffers->shapeFlags.size()){
	NvFlexSetShapes(
				g_flex,
				g_buffers->shapeGeometry.buffer,
				g_buffers->shapePositions.buffer,
				g_buffers->shapeRotations.buffer,
				g_buffers->shapePrevPositions.buffer,
				g_buffers->shapePrevRotations.buffer,
				g_buffers->shapeFlags.buffer,
				g_buffers->shapeFlags.size());
		}
		

}


void FlexSystem::emission(){

	size_t e=0;


	for (; e < nEmitter; ++e)
	{
		if (!g_rectEmitters[e].mEnabled)
			continue;

		Vec3 emitterRot = g_rectEmitters[e].mRot;
		Point3 emitterSize = g_rectEmitters[e].mSize;
		Point3 emitterPos = g_rectEmitters[e].mPos;

		float r;
		int phase;

		if (g_params.fluid)
		{
			r = g_params.fluidRestDistance;
			phase = NvFlexMakePhase(0, eNvFlexPhaseSelfCollide | eNvFlexPhaseFluid);
		}
		else
		{
			r = g_params.solidRestDistance;
			phase = NvFlexMakePhase(0, eNvFlexPhaseSelfCollide);
		}

		float numSlices = (g_rectEmitters[e].mSpeed / r)*g_dt;

		// whole number to emit
		int n = int(numSlices + g_rectEmitters[e].mLeftOver);
				
		if (n)
			g_rectEmitters[e].mLeftOver = (numSlices + g_rectEmitters[e].mLeftOver)-n;
		else
			g_rectEmitters[e].mLeftOver += numSlices;

		//int circle = 1;

		int disc = g_rectEmitters[e].mDisc;

		int dy = 0;

				
		int dx = int(ceilf(emitterSize.x / g_params.fluidRestDistance));
		if(disc)
			dy = int(ceilf(emitterSize.x / g_params.fluidRestDistance));
		else
			dy = int(ceilf(emitterSize.y / g_params.fluidRestDistance));
		Mat44	tMat = TransformMatrix(Rotation(emitterRot.y, emitterRot.z, emitterRot.x), emitterPos);


		for (int z=0; z < n; ++z)
				{
			for (int x=0; x < dx; ++x)
			{
				for (int y=0; y < dy; ++y)
				{
						
					Point3 position = g_params.fluidRestDistance*Point3(float(x) - 0.5*(dx-1), float(y) - 0.5*(dy-1), float(z)) + RandomUnitVector()*g_params.fluidRestDistance*0.01f;

					int keep = 1;

					if(disc){
						if(position.x*position.x + position.y*position.y>0.25*emitterSize.x*emitterSize.x)
							keep=0;
					}

					if(g_rectEmitters[e].mNoise){
						Point3 scaledP = position*g_rectEmitters[e].mNoiseFreq + Point3(0,0,g_rectEmitters[e].mNoiseOffset);
						const float kNoise = Perlin3D(scaledP.x, scaledP.y, scaledP.z, 1, 0.25f);

						if(kNoise<g_rectEmitters[e].mNoiseThreshold)
							keep=0;

					}

					if(keep) {

						position = tMat*position;

						Vec3 vel = Vec3(0,0,g_rectEmitters[e].mSpeed);
						vel = tMat*vel;

						g_buffers->positions[cursor] = Vec4(Vec3(position), 1.0f);
						g_buffers->velocities[cursor] = vel;
						g_buffers->phases[cursor] = phase;

						if(g_buffers->activeIndices.size()<maxParticles)
							g_buffers->activeIndices.push_back(cursor);

						if(cursor<maxParticles-1)
							cursor++;
						else
							cursor = 0;

						
					}//end dist

				}
			}
		}

	}

	

}


void FlexSystem::update(){


		activeParticles = NvFlexGetActiveCount(g_flex);


		time1 = GetSeconds();

		emission();

		time2 = GetSeconds();


		g_buffers->UnmapBuffers();


		if (g_buffers->activeIndices.size()) {
			NvFlexSetActive(g_flex, g_buffers->activeIndices.buffer, g_buffers->activeIndices.size());
		}

		if (g_buffers->positions.size()) {


			NvFlexSetParticles(g_flex, g_buffers->positions.buffer, g_buffers->positions.size());
			NvFlexSetVelocities(g_flex, g_buffers->velocities.buffer, g_buffers->velocities.size());
			NvFlexSetPhases(g_flex, g_buffers->phases.buffer, g_buffers->phases.size());

		}

		setShapes();

		NvFlexSetParams(g_flex, &g_params);
		NvFlexExtSetForceFields(g_forcefieldCallback, &g_forcefield, 1);

		NvFlexUpdateSolver(g_flex, g_dt, g_numSubsteps, g_profile);


		if(g_buffers->positions.size()){


			NvFlexGetParticles(g_flex, g_buffers->positions.buffer, g_buffers->positions.size());
			NvFlexGetVelocities(g_flex, g_buffers->velocities.buffer, g_buffers->velocities.size());

		}

		activeParticles = NvFlexGetActiveCount(g_flex);

}

void FlexSystem::CreateSpring(int i, int j, float stiffness, float give)
{
	g_buffers->springIndices.push_back(i);
	g_buffers->springIndices.push_back(j);
	g_buffers->springLengths.push_back((1.0f + give) * Length(Vec3(g_buffers->positions[i]) - Vec3(g_buffers->positions[j])));
	g_buffers->springStiffness.push_back(stiffness);
}

float FlexSystem::SampleSDF(const float *sdf, int dim, int x, int y, int z)
{
	assert(x < dim && x >= 0);
	assert(y < dim && y >= 0);
	assert(z < dim && z >= 0);

	return sdf[z * dim * dim + y * dim + x];
}

// return normal of signed distance field
Vec3 FlexSystem::SampleSDFGrad(const float *sdf, int dim, int x, int y, int z)
{
	int x0 = max(x - 1, 0);
	int x1 = min(x + 1, dim - 1);

	int y0 = max(y - 1, 0);
	int y1 = min(y + 1, dim - 1);

	int z0 = max(z - 1, 0);
	int z1 = min(z + 1, dim - 1);

	float dx = (SampleSDF(sdf, dim, x1, y, z) - SampleSDF(sdf, dim, x0, y, z)) * (dim * 0.5f);
	float dy = (SampleSDF(sdf, dim, x, y1, z) - SampleSDF(sdf, dim, x, y0, z)) * (dim * 0.5f);
	float dz = (SampleSDF(sdf, dim, x, y, z1) - SampleSDF(sdf, dim, x, y, z0)) * (dim * 0.5f);

	return Vec3(dx, dy, dz);
}