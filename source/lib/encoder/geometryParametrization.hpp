#pragma once
#include <chrono>
#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <memory>
#include <cmath>
#include "geometryDecimate.hpp"

#include "misc.hpp"
#include "vector.hpp"
#include "mesh.hpp"
#include "kdtree.hpp"

//============================================================================

struct SubdivFitParameters {
	uint32_t subdivFitSubdivIterCount = 3;
	uint32_t samplingSubdivIterCount = 3;
	uint32_t fittingIterCount = 16;
	uint32_t adaptFittingIterMode = 1;
	bool     removeDuplicateVert = true;
	bool     fitSubdivSurface = true;
	double   smoothingCoeff = 0.25;
	double   smoothingCoeffDecayRatio = 0.75;
	double   missedVertSmoothingCoeff = 0.1;
	uint32_t missedVertSmoothingIterCount = 10;
	double   initialDeformNormalDeviationThreshold = 0.1;
	uint32_t initialDeformNNCount = 1;
	bool     initialDeformForceNormalDisp = false;
	bool     applySmoothingDeform = true;
	bool     smoothDeformUpdateNormals = true;
	bool     smoothingDeformUseInitialGeometry = true;
	bool     smoothingDeformSmoothMotion = true;
	double   smoothDeformTriangleNormalFlipThreshold = -0.5;
};

//============================================================================
class SubdivFit {
public:
	SubdivFit() = default;
	~SubdivFit() = default;

	template<typename T>
	bool generate(TriMesh<T>&                 refMesh,     // target
				  const TriMesh<T>&           deciMesh,    // source
				  TriMesh<T>&                 mapMesh,     // mapped
				  TriMesh<T>&                 baseMesh,    // base
				  TriMesh<T>&                 subdivMesh,  // deformed
				  const SubdivFitParameters&  params);

	template<typename T>
	bool subDivide(TriMesh<T>& mesh, const SubdivFitParameters& params);

	template<typename T>
	bool
		one2OneMapGen(const TriMesh<T>& target, TriMesh<T>& mapped, TriMesh<T>& deformed, const SubdivFitParameters& params);

	template<typename T>
	void one2OneMapGen(const TriMesh<T>& target, const KdTree<T>& kdtree, TriMesh<T>& output);

	template<typename T, bool isForceNormalDisp>
	void initMapping(const TriMesh<T>& target, TriMesh<T>& output, const SubdivFitParameters& params);

	template<typename T, bool isForceNormalDisp>
	void initMapping(const TriMesh<T>&          target,
					 const TriMesh<T>&          mapped,
					 TriMesh<T>&                output,
					 const SubdivFitParameters& params);

};



