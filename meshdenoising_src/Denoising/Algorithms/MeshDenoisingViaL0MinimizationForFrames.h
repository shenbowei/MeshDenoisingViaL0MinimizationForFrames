#ifndef MESHDENOISINGVIAL0MINIMIZATIONFORFRAMES_H
#define MESHDENOISINGVIAL0MINIMIZATIONFORFRAMES_H


#include "MeshDenoisingBase.h"
#include <vector>
#include <map>
#include <set>
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/SparseQR"
#include "Eigen/Eigen"

#include <ANN/ANN.h>

class MeshDenoisingViaL0MinimizationForFrames : public MeshDenoisingBase
{
public:
    MeshDenoisingViaL0MinimizationForFrames(DataManager *_data_manager, ParameterSet *_parameter_set);
    ~MeshDenoisingViaL0MinimizationForFrames() {}

private:
    void denoise();
    void initParameters();

    double getAverageDihedralAngle(TriMesh &mesh);
    void calculateAreaBasedEdgeOperator(TriMesh &mesh,
                                        std::vector<TriMesh::Point> &area_based_edge_operator,
                                        std::vector<std::vector<TriMesh::VertexHandle> > &edge_vertex_handle,
                                        std::vector< std::vector<double> > &coef);
    void solveDelta(std::vector<TriMesh::Point> &area_based_edge_operator, double lambda, double beta,
                    std::vector<TriMesh::Point> &delta);
    void getInitialVerticesMatrix(TriMesh &noisy_mesh, Eigen::MatrixXd &initial_vertices_matrix);
    void solveVertices(TriMesh &mesh, Eigen::MatrixXd &initial_vertices_matrix,
                       std::vector< std::vector<TriMesh::VertexHandle> > &edge_vertex_handle,
                       std::vector< std::vector<double> > &coef, std::vector<TriMesh::Point> &delta,
                       double alpha, double beta);
    //邻帧的比重系数
    double ratioPreFrame;
    double ratioNextFrame;
    double threshold;
    double annDistanceTotal[2];
    //前后帧的网格
    TriMesh meshPreFrame;
    TriMesh meshNextFrame;
    //前后帧对应的Ann搜索点集
    ANNpointArray annCloudPreFrame;
    ANNpointArray annCloudNextFrame;
    //annSearch中需要使用的变量
    ANNidxArray idxArr;
    ANNdistArray distArr;
    ANNkd_tree*	kdTreePreFrame;
    ANNkd_tree*	kdTreeNextFrame;
    int k;

    int annSearch(ANNkd_tree *kdTree, ANNpoint& queryPt);
    void initAnn();
    double getDistance(ANNpoint p1, ANNpoint p2);
    bool readMesh(TriMesh &mesh, const QString &fileName);

};

#endif // MESHDENOISINGVIAL0MINIMIZATIONFORFRAMES_H
