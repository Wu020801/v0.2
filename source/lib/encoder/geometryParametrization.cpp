#include <chrono>
#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <memory>
#include <cmath>

#include "geometryParametrization.hpp"
#include "geometryDecimate.hpp"

//============================================================================

template<typename T>
bool SubdivFit::generate(TriMesh<T>&                 refMesh,     // reference
						  const TriMesh<T>&           deciMesh,    // decimate
						  TriMesh<T>&                 mapMesh,     // mapped
						  TriMesh<T>&                 baseMesh,    // base
						  TriMesh<T>&                 subdivFitMesh,  // submesh
						  const SubdivFitParameters&  params) {

	subdivFitMesh = deciMesh;
	subdivFitMesh.Normal.clear();
	subdivFitMesh.FaceNormal.clear();
	subDivide(subdivFitMesh, params);
	// fit subdivided mesh
	one2OneMapGen(refMesh, mapMesh, subdivFitMesh, params);

	baseMesh = deciMesh;
	return true;
}

//============================================================================
template bool SubdivFit::generate<float>(TriMesh<float>&            refMesh,
                                         TriMesh<float> const&      deciMesh,
                                         TriMesh<float>&            mapMesh,
                                         TriMesh<float>&            baseMesh,
                                         TriMesh<float>&            subdivMesh,
										 const SubdivFitParameters& params);

//============================================================================
template bool SubdivFit::generate<double>(TriMesh<double>&            refMesh,
										 TriMesh<double> const&      deciMesh,
										 TriMesh<double>&            mapMesh,
										 TriMesh<double>&            baseMesh,
										 TriMesh<double>&            subdivMesh,
										 const SubdivFitParameters& params);

//============================================================================

template<typename T>
bool SubdivFit::subDivide(TriMesh<T>&          mesh,
						  const SubdivFitParameters& params) {

	auto start = std::chrono::steady_clock::now();
	mesh.computeNormals();
	if (params.subdivFitSubdivIterCount != 0) {
		std::vector<SubdivLoD> levelOfDetails;
		std::vector<uint64_t>  subdivEdges;
		std::vector<uint64_t>  subdivUvEdges;
		mesh.subdivMidpoint(params.subdivFitSubdivIterCount, &levelOfDetails, &subdivEdges, &subdivUvEdges);
		mesh.Normal.resize(mesh.XYZ.size());
		subdivInterpolate<T, true>(mesh.Normal, levelOfDetails, subdivEdges);
	}

	auto end = std::chrono::steady_clock::now();
	auto deltams = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

#if DEBUG_ENC_LOG
	std::cout << " [done] " << deltams.count() << " ms\n";
	std::cout << "\t Subdivided mesh: " << mesh.XYZ.size() << "V " << mesh.FaceXYZ.size() << "T\n";
#endif

	return true;
}

//============================================================================

template<typename T>
int32_t
findNearestNeighbor(const KdTree<T>&       kdtTarget,
					const TriMesh<T>&      target,
					const Vector3D<T>&     point,
					const Vector3D<T>&     normal,
					std::vector<uint32_t>& indexes,
					std::vector<T>         dists,
					const int32_t          numDeform,
					double                 deformNormalThreshold) {
	kdtTarget.query(point.data(), numDeform, indexes.data(), dists.data());
	// int32_t nIdx = -1;
	for (size_t i = 0; i < numDeform; ++i) {
		const auto& idx = indexes[i];
		const auto& n = target.Normal[idx];
		if (normal * n >= deformNormalThreshold) { return idx; }
	}
	return -1;
}

//============================================================================

template<typename T>
void
findClosestPoint(const TriMesh<T>&  mapped,
				const TriMesh<T>&  target,
				const AdjInfo&     vert2FaceTarget,
				const Vector3D<T>& point,
				Vector3D<T>&       nearestXYZ,
				Vector3D<T>&       nearestDispl,
				int32_t            nIdx) {
	const auto  neighboursTarget = vert2FaceTarget.Neighbors;
	T           minDist = std::numeric_limits<T>::max();
	Vector3D<T> bcoords = Vector3D<T>(0, 0, 0);
	const auto  faceIdx = vert2FaceTarget.getNeighbors(nIdx);
	for (const auto& fidx : faceIdx) {
		const auto& face = mapped.FaceXYZ[fidx];
		const auto  closestPoint =
			closestPointInTriangle(point, mapped.XYZ[face[0]], mapped.XYZ[face[1]], mapped.XYZ[face[2]], &bcoords);
		const auto dist = (closestPoint - point).L2normSq();
		if (dist < minDist) {
			minDist = dist;
			nearestXYZ = closestPoint;
			nearestDispl = 0;
			for (size_t i = 0; i < 3; ++i) { nearestDispl += bcoords[i] * (target.XYZ[face[i]] - mapped.XYZ[face[i]]); }
		}
	}
}

//============================================================================

template<typename T, bool isForceNormalDisp>
void
SubdivFit::initMapping(const TriMesh<T>& target, TriMesh<T>& output, const SubdivFitParameters& params) {
	// initialization
	AdjInfo vert2FaceTarget;
	computeVertex2triangle(target.XYZ.size(), target.FaceXYZ, vert2FaceTarget);
	KdTree<T> kdtTarget(3, target.XYZ, 10);

	const auto            numDeform = params.initialDeformNNCount;
	std::vector<uint32_t> indexes(numDeform);
	std::vector<T>        dists(numDeform);
	std::vector<int32_t>  missedVertIdx;

	for (size_t vidx = 0, numXYZ = output.XYZ.size(); vidx < numXYZ; ++vidx) {
		const Vector3D<T>& point = output.XYZ[vidx];
		const Vector3D<T>& normal = output.Normal[vidx];

		// 1. find the nearest neighbor
		int32_t nIdx = findNearestNeighbor(
			kdtTarget, target, point, normal, indexes, dists, numDeform, params.initialDeformNormalDeviationThreshold);
		if (nIdx == -1) {
			missedVertIdx.push_back(vidx);
			continue;
		}

		// 2. find the closest point on the triangle
		Vector3D<T> nearestXYZ = target.XYZ[nIdx];
		Vector3D<T> nearestDispl = Vector3D<T>(0, 0, 0);
		findClosestPoint(target, target, vert2FaceTarget, point, nearestXYZ, nearestDispl, nIdx);

		// 3. update point
		if constexpr (isForceNormalDisp) {
			Vector3D<T> delta = nearestXYZ + nearestDispl - point;
			delta = delta * normal;
			output.XYZ[vidx] = point + delta * normal;
		}
		else {
			output.XYZ[vidx] = nearestXYZ + nearestDispl;
		}
	}
}

//============================================================================

template<typename T, bool isForceNormalDisp>
void
SubdivFit::initMapping(const TriMesh<T>&          target,
					   const TriMesh<T>&          mapped,
					   TriMesh<T>&                output,
					   const SubdivFitParameters& params) {
	// initialization
	AdjInfo vert2FaceMapped;
	computeVertex2triangle(mapped.XYZ.size(), mapped.FaceXYZ, vert2FaceMapped);
	KdTree<T>             kdtMapped(3, mapped.XYZ, 10);
	const auto            numDeform = params.initialDeformNNCount;
	std::vector<uint32_t> indexes(numDeform);
	std::vector<T>        dists(numDeform);

	for (size_t vidx = 0, numXYZ = output.XYZ.size(); vidx < numXYZ; ++vidx) {
		const Vector3D<T>& point = output.XYZ[vidx];
		const Vector3D<T>& normal = output.Normal[vidx];

		// 1. find the nearest neighbor
		int32_t nIdx = findNearestNeighbor(
			kdtMapped, target, point, normal, indexes, dists, numDeform, params.initialDeformNormalDeviationThreshold);
		if (nIdx == -1) { continue; }

		// 2. find the closest point on the triangle
		Vector3D<T> nearestXYZ = mapped.XYZ[nIdx];
		Vector3D<T> nearestDispl = target.XYZ[nIdx] - mapped.XYZ[nIdx];
		T           minDist = std::numeric_limits<T>::max();
		const auto  faceIdx = vert2FaceMapped.getNeighbors(nIdx);
		int         count = 0;
		for (const auto& fidx : faceIdx) {
			const auto& face = mapped.FaceXYZ[fidx];
			Vector3D<T> bcoords = Vector3D<T>(0, 0, 0);
			const auto& pt0 = mapped.XYZ[face[0]];
			const auto& pt1 = mapped.XYZ[face[1]];
			const auto& pt2 = mapped.XYZ[face[2]];
			const auto  closestPoint = closestPointInTriangle(point, pt0, pt1, pt2, &bcoords);
			const auto  dist = (closestPoint - point).L2normSq();

			if (dist < minDist) {
				minDist = dist;
				nearestXYZ = closestPoint;
				nearestDispl = bcoords[0] * (target.XYZ[face[0]] - pt0) + bcoords[1] * (target.XYZ[face[1]] - pt1)
					+ bcoords[2] * (target.XYZ[face[2]] - pt2);
			}
		}

		// 3. update point
		if constexpr (isForceNormalDisp) {
			Vector3D<T> delta = nearestXYZ + nearestDispl - point;
			delta = delta * normal;
			output.XYZ[vidx] = point + delta * normal;
		}
		else {
			output.XYZ[vidx] = nearestXYZ + nearestDispl;
		}
	}
}

//============================================================================


template<typename T>
bool
needInitMap(const TriMesh<T> target, const TriMesh<T> mapped) {
	return mapped.XYZ.size() != target.XYZ.size();
}

//============================================================================

template<typename T>
void
SubdivFit::one2OneMapGen(const TriMesh<T>& target, const KdTree<T>& kdtree, TriMesh<T>& output) {
	for (size_t vidx = 0, numXYZ = output.XYZ.size(); vidx < numXYZ; ++vidx) {
		uint32_t           closestIdx = 0;
		T                  minDist = std::numeric_limits<T>::max();
		const Vector3D<T>& point = output.XYZ[vidx];
		const Vector3D<T>& normal = output.Normal[vidx];
		kdtree.query(point.data(), 1, &closestIdx, &minDist);
		const auto delta = (target.XYZ[closestIdx] - point) * normal;
		output.XYZ[vidx] = point + delta * normal;
	}
}

//============================================================================

template<typename T>
bool
SubdivFit::one2OneMapGen(const TriMesh<T>&          target,    // refMesh
						 TriMesh<T>&                mapped,    // mapMesh
						 TriMesh<T>&                deformed,  // subdivMesh
						 const SubdivFitParameters& params) {

	auto start = std::chrono::steady_clock::now();

	// 1. Pre-processing - smoothing
	std::vector<Vector3D<T>> initXYZ;
	std::vector<Vector3D<T>> initFaceNormals;

	if (params.smoothingDeformUseInitialGeometry) {
		initXYZ = deformed.XYZ;
		deformed.computeTriangleNormals(initFaceNormals);
	}

	// 2. Initial one-to-one mapping it it doesn't exist
	if (needInitMap(target, mapped)) {
		mapped = target;
		if (params.initialDeformForceNormalDisp) {
			initMapping<T, true>(deformed, mapped, params);
		}
		else {
			initMapping<T, false>(deformed, mapped, params);
		}
	}

	// 3. Initial one-to-one mapping with subdivision
	auto subdivTarget = target;
	auto subdivMapped = mapped;

	if (params.initialDeformForceNormalDisp) {
		initMapping<T, true>(subdivTarget, subdivMapped, deformed, params);
	}
	else {
		initMapping<T, false>(subdivTarget, subdivMapped, deformed, params);
	}

	auto end = std::chrono::steady_clock::now();
	auto deltams = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

	return true;
}

//============================================================================
