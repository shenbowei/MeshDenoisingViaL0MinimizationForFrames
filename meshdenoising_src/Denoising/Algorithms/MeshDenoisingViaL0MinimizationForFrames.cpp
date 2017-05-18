#include "MeshDenoisingViaL0MinimizationForFrames.h"
#include <QDebug>
#include <iostream>
#include <QException>

#define WITHFRAME

MeshDenoisingViaL0MinimizationForFrames::MeshDenoisingViaL0MinimizationForFrames(DataManager *_data_manager, ParameterSet *_parameter_set)
    : MeshDenoisingBase(_data_manager, _parameter_set)
{
    initParameters();
}

void MeshDenoisingViaL0MinimizationForFrames::denoise()
{
    // get data
    TriMesh mesh = data_manager_->getNoisyMesh();

    if(mesh.n_vertices() == 0)
        return;

    // get parameters
    double mu_beta, beta, beta_max;
    if(!parameter_set_->getValue(QString("mu_beta"), mu_beta))
        return;
    if(!parameter_set_->getValue(QString("beta"), beta))
        return;
    if(!parameter_set_->getValue(QString("beta_max"), beta_max))
        return;
    double mu_alpha, alpha;
    if(!parameter_set_->getValue(QString("mu_alpha"), mu_alpha))
        return;
    if(!parameter_set_->getValue(QString("alpha"), alpha))
        return;
    double lambda;
    if(!parameter_set_->getValue(QString("lambda"), lambda))
        return;
    QString preMeshFile;
    if(!parameter_set_->getValue(QString("pre_frame_mesh"), preMeshFile))
        return;
    QString nextMeshFile;
    if(!parameter_set_->getValue(QString("next_frame_mesh"), nextMeshFile))
        return;
    if(!parameter_set_->getValue(QString("pre_frame_ratio"), ratioPreFrame))
        return;
    if(!parameter_set_->getValue(QString("next_frame_ratio"), ratioNextFrame))
        return;
    if(!parameter_set_->getValue(QString("threshold"), threshold))
        return;

    //test
    //parameter_set_->print();
    //return;
    /*
    qDebug()<<"preMeshFile:"<<preMeshFile;
    qDebug()<<"nextMeshFile:"<<nextMeshFile;
    qDebug()<<"mu_beta:"<<mu_beta;
    qDebug()<<"beta:"<<beta;
    qDebug()<<"beta_max:"<<beta_max;
    qDebug()<<"mu_alpha:"<<mu_alpha;
    qDebug()<<"alpha:"<<alpha;
    qDebug()<<"lambda:"<<lambda;
    return;
    */

    //连续帧去噪的准备工作
    if(!readMesh(meshPreFrame, preMeshFile))
        return;
    if(!readMesh(meshNextFrame, nextMeshFile))
        return;
    initAnn();

    // do optimization
    Eigen::MatrixXd initial_vertices_matrix;
    getInitialVerticesMatrix(mesh, initial_vertices_matrix);

    std::vector<TriMesh::Point> area_based_edge_operator;//D(e)in(5)
    std::vector< std::vector<TriMesh::VertexHandle> > edge_vertex_handle;//D（e）式中右侧四个点的序号
    std::vector< std::vector<double> > coef;//D(e)式中左侧的向量
    std::vector<TriMesh::Point> delta;

    while(beta < beta_max){
        calculateAreaBasedEdgeOperator(mesh, area_based_edge_operator, edge_vertex_handle, coef);
        // update delta
        solveDelta(area_based_edge_operator, lambda, beta, delta);
        // update vertices
        solveVertices(mesh, initial_vertices_matrix, edge_vertex_handle, coef, delta, alpha, beta);
        beta *= mu_beta;
        alpha *= mu_alpha;
    }

    qDebug() << "L0 denoising finished!!!" ;
    qDebug() << "annDistanceTotal[0] = " << this->annDistanceTotal[0];
    qDebug() << "annDistanceTotal[1] = " << this->annDistanceTotal[1];
    //保存结果
    //OpenMesh::IO::Options opt = 0x0020;
    //OpenMesh::IO::write_mesh(mesh, "ply/sbw/00643_l0denoised_0.001.ply", opt);
    //OpenMesh::IO::write_mesh(mesh, "ply/sbw/00031_addNoise_zhui_l0denoised_0.001.ply", opt);
    // update data
    data_manager_->setMesh(mesh);
    data_manager_->setDenoisedMesh(mesh);
    emit(statusMessage("Applying algorithm --MeshDenoisingViaL0MinimizationForFrames-- done."));
}

void MeshDenoisingViaL0MinimizationForFrames::initParameters()
{
    parameter_set_->removeAllParameter();

    parameter_set_->addParameter(QString("mu_beta"), 1.414, QString("mu_beta"), QString("Update ratio for beta."),
                                 true, 1.001, 10000.0);
    parameter_set_->addParameter(QString("beta"), 0.001, QString("beta"), QString("Initial beta value in optimization."),
                                 true, 1.0e-9, 10000.0);
    parameter_set_->addParameter(QString("beta_max"), 1000.0, QString("beta_max"), QString("Max beta value in optimization."),
                                 true, 10.0, 1.0e9);
    parameter_set_->addParameter(QString("mu_alpha"), 0.5, QString("mu_alpha"), QString("Update ratio for alpha."),
                                 true, 1.0e-9, 0.9999);

    TriMesh mesh = data_manager_->getNoisyMesh();
    if(mesh.n_vertices() == 0)
    {
        parameter_set_->addParameter(QString("alpha"), 0.0001, QString("alpha"), QString("Initial alpha value in optimization."),
                                     true, 0.0, 1.0e9);
        parameter_set_->addParameter(QString("lambda"), 0.0001, QString("lambda"), QString("Lambda value in optimization."),
                                     true, 1.0e-9, 1.0e9);
    }
    else
    {
        double mean_edge_length = getAverageEdgeLength(mesh);
        double mean_dihedral_angle = getAverageDihedralAngle(mesh);
        parameter_set_->addParameter(QString("alpha"), 0.1 * mean_dihedral_angle, QString("alpha"), QString("Initial alpha value in optimization."),
                                     true, 0.0, 1.0e9);
        parameter_set_->addParameter(QString("lambda"), 0.2 * mean_edge_length * mean_edge_length * mean_dihedral_angle, QString("lambda"), QString("Lambda value in optimization."),
                                     true, 1.0e-9, 1.0e9);
    }

    //For Frames Parameters
    parameter_set_->addParameter(QString("pre_frame_mesh"), QString(""), QString("pre_frame_mesh"), QString("Pre-frame mesh."));
    parameter_set_->addParameter(QString("next_frame_mesh"), QString(""), QString("next_frame_mesh"), QString("Next-frame mesh."));

    parameter_set_->addParameter(QString("threshold"), double(0.2), QString("threshold"), QString("Threshold value between pre/next frames."),
                                 true, 1.0e-9, 10000.0);
    parameter_set_->addParameter(QString("pre_frame_ratio"), double(0.5), QString("pre_frame_ratio"), QString("Pre-frame mesh ratio."),
                                 true, 1.0e-9, 10000.0);
    parameter_set_->addParameter(QString("next_frame_ratio"), double(0.5), QString("next_frame_ratio"), QString("Next-frame mesh ratio."),
                                 true, 1.0e-9, 10000.0);

    parameter_set_->setName(QString("Mesh Denoising via L0 Minimization"));
    parameter_set_->setLabel(QString("Mesh Denoising via L0 Minimization"));
    parameter_set_->setIntroduction(QString("Mesh Denoising via L0 Minimization -- Parameters"));
}

double MeshDenoisingViaL0MinimizationForFrames::getAverageDihedralAngle(TriMesh &mesh)
{
    mesh.request_face_normals();
    mesh.update_face_normals();

    double mean_dihedral_angle = 0.0;
    double num = 0.0;
    for(TriMesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); e_it++)
    {
        if(!mesh.is_boundary(*e_it))
            mean_dihedral_angle += mesh.calc_dihedral_angle(*e_it);
        num++;
    }

    return mean_dihedral_angle / num;
}

void MeshDenoisingViaL0MinimizationForFrames::calculateAreaBasedEdgeOperator(TriMesh &mesh,
                                                                    std::vector<TriMesh::Point> &area_based_edge_operator,
                                                                    std::vector< std::vector<TriMesh::VertexHandle> > &edge_vertex_handle,
                                                                    std::vector< std::vector<double> > &coef)
{
    std::vector<double> face_area;
    getFaceArea(mesh, face_area);

    area_based_edge_operator.resize((int)mesh.n_edges(), TriMesh::Point(0.0, 0.0, 0.0));
    std::vector<double> temp_coef(4, 0.0);
    coef.resize(mesh.n_edges(), temp_coef);
    std::vector<TriMesh::VertexHandle> vertex_handle(4);
    edge_vertex_handle.resize(mesh.n_edges(), vertex_handle);
    for(TriMesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); e_it++)
    {
        if(!mesh.is_boundary(*e_it))
        {
            int index = e_it->idx();
            double edge_length = mesh.calc_edge_length(*e_it);

            // get four vertices correspond to edge *e_it
            TriMesh::HalfedgeHandle he = mesh.halfedge_handle(*e_it, 0);
            TriMesh::VertexHandle v1 = mesh.from_vertex_handle(he);
            TriMesh::VertexHandle v3 = mesh.to_vertex_handle(he);
            TriMesh::HalfedgeHandle he_next = mesh.next_halfedge_handle(he);
            TriMesh::VertexHandle v4 = mesh.to_vertex_handle(he_next);

            TriMesh::HalfedgeHandle he_oppo = mesh.opposite_halfedge_handle(he);
            TriMesh::HalfedgeHandle he_oppo_next = mesh.next_halfedge_handle(he_oppo);
            TriMesh::VertexHandle v2 = mesh.to_vertex_handle(he_oppo_next);

            // two faces
            TriMesh::FaceHandle f1 = mesh.face_handle(he);
            TriMesh::FaceHandle f2 = mesh.face_handle(he_oppo);

            // the area of two faces correspond to edge *e_it
            double area134 = face_area[f1.idx()];
            double area123 = face_area[f2.idx()];
            double totalArea = area123 + area134;

            TriMesh::Point p1 = mesh.point(v1);
            TriMesh::Point p2 = mesh.point(v2);
            TriMesh::Point p3 = mesh.point(v3);
            TriMesh::Point p4 = mesh.point(v4);

            TriMesh::Point p12 = p1 - p2;
            TriMesh::Point p13 = p1 - p3;
            TriMesh::Point p14 = p1 - p4;
            TriMesh::Point p23 = p2 - p3;
            TriMesh::Point p34 = p3 - p4;

            // calc coefficient
            temp_coef[0] = (area123 * (p34 | p13) - area134 * (p13 | p23)) / (edge_length * edge_length * totalArea);
            temp_coef[1] = area134 / totalArea;
            temp_coef[2] = (-area123 * (p13 | p14) - area134 * (p12 | p13)) / (edge_length * edge_length * totalArea);
            temp_coef[3] = area123 / totalArea;
            if(_isnan(temp_coef[0]))
            {
                qDebug()<<"temp_coef[0] nan:"<<edge_length<<" "<<totalArea;
                temp_coef[0] = 0.0;
            }
            if(_isnan(temp_coef[1]))
            {
                qDebug()<<"temp_coef[1] nan:"<<totalArea;
                temp_coef[1] = 0.0;
            }
            if(_isnan(temp_coef[2]))
            {
                qDebug()<<"temp_coef[2] nan:"<<edge_length<<" "<<totalArea;
                temp_coef[2] = 0.0;
            }
            if(_isnan(temp_coef[3]))
            {
                qDebug()<<"temp_coef[3] nan:"<<totalArea;
                temp_coef[3] = 0.0;
            }
            coef[index] = temp_coef;

            vertex_handle[0] = v1;
            vertex_handle[1] = v2;
            vertex_handle[2] = v3;
            vertex_handle[3] = v4;
            edge_vertex_handle[index] = vertex_handle;

            // calc area-based edge operator
            TriMesh::Point pt = p1 * temp_coef[0] + p2 * temp_coef[1] + p3 * temp_coef[2] + p4 * temp_coef[3];
            area_based_edge_operator[index] = pt;
        }
    }
}

void MeshDenoisingViaL0MinimizationForFrames::solveDelta(std::vector<TriMesh::Point> &area_based_edge_operator, double lambda, double beta,
                                                std::vector<TriMesh::Point> &delta)
{
    delta.resize((int)area_based_edge_operator.size(), TriMesh::Point(0.0, 0.0, 0.0));

    for(int i = 0; i < (int)area_based_edge_operator.size(); i++)
    {
        TriMesh::Point pt = area_based_edge_operator[i];
        if(pt.length() * pt.length() >= lambda/beta)
            delta[i] = pt;
    }
}

void MeshDenoisingViaL0MinimizationForFrames::getInitialVerticesMatrix(TriMesh &noisy_mesh, Eigen::MatrixXd &initial_vertices_matrix)
{
    initial_vertices_matrix.resize((int)noisy_mesh.n_vertices(), 3);

    for(TriMesh::VertexIter v_it = noisy_mesh.vertices_begin(); v_it != noisy_mesh.vertices_end(); v_it++)
    {
        int index = v_it->idx();
        TriMesh::Point p = noisy_mesh.point(*v_it);
        initial_vertices_matrix(index, 0) = p[0];
        initial_vertices_matrix(index, 1) = p[1];
        initial_vertices_matrix(index, 2) = p[2];
    }
}

void MeshDenoisingViaL0MinimizationForFrames::solveVertices(TriMesh &mesh, Eigen::MatrixXd &initial_vertices_matrix,
                                                   std::vector< std::vector<TriMesh::VertexHandle> > &edge_vertex_handle,
                                                   std::vector< std::vector<double> > &coef, std::vector<TriMesh::Point> &delta,
                                                   double alpha, double beta)
{
    Eigen::MatrixXd right_term = initial_vertices_matrix;
    Eigen::SparseMatrix<double> coef_matrix((int)mesh.n_vertices(), (int)mesh.n_vertices());

    std::vector< Eigen::Triplet<double> > triple; triple.clear();

    std::map<TriMesh::VertexHandle, double> vertex_coef;

    std::set<TriMesh::EdgeHandle> edge_handle;;
    ANNpoint queryPt = annAllocPt(3);
    for(TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
    {
        edge_handle.clear();
        vertex_coef.clear();

        //这里加上前后帧的比例系数
        vertex_coef[*v_it] = 1.0 + this->ratioNextFrame + this->ratioPreFrame;

        for(TriMesh::VertexFaceIter vf_it = mesh.vf_iter(*v_it); vf_it.is_valid(); vf_it++)
        {
            for(TriMesh::FaceEdgeIter fe_it = mesh.fe_iter(*vf_it); fe_it.is_valid(); fe_it++)
            {
                edge_handle.insert(*fe_it);
            }
        }

        TriMesh::Point right(0.0, 0.0, 0.0);
        for(std::set<TriMesh::EdgeHandle>::iterator s_it = edge_handle.begin(); s_it != edge_handle.end(); s_it++)
        {
            if(!mesh.is_boundary(*s_it))//判断*s_it是否为边缘
            {
                int index = (*s_it).idx();
                TriMesh::VertexHandle v1 = edge_vertex_handle[index][0],
                        v2 = edge_vertex_handle[index][1],
                        v3 = edge_vertex_handle[index][2],
                        v4 = edge_vertex_handle[index][3];
                double coe1 = coef[index][0],
                        coe2 = coef[index][1],
                        coe3 = coef[index][2],
                        coe4 = coef[index][3];
                TriMesh::Point temp_delta = delta[index];
                if(v1 == *v_it)
                {
                    vertex_coef[v1] = vertex_coef[v1] + alpha + beta * coe1 * coe1;
                    if(_isnan(vertex_coef[v1]))
                    {
                        qDebug()<<"if1 v1:"<<vertex_coef[v1]<<" "<<coe1<<" "<<coe1;
                    }
                    vertex_coef[v2] = vertex_coef[v2] - alpha + beta * coe1 * coe2;
                    if(_isnan(vertex_coef[v2]))
                    {
                        qDebug()<<"if1 v2:"<<vertex_coef[v2]<<" "<<coe1<<" "<<coe2;
                    }
                    vertex_coef[v3] = vertex_coef[v3] + alpha + beta * coe1 * coe3;
                    if(_isnan(vertex_coef[v3]))
                    {
                        qDebug()<<"if1 v3:"<<vertex_coef[v3]<<" "<<coe1<<" "<<coe3;
                    }
                    vertex_coef[v4] = vertex_coef[v4] - alpha + beta * coe1 * coe4;
                    if(_isnan(vertex_coef[v4]))
                    {
                        qDebug()<<"if1 v4:"<<vertex_coef[v4]<<" "<<coe1<<" "<<coe4;
                    }
                    right += temp_delta * beta * coe1;
                }
                else if(v2 == *v_it)
                {
                    vertex_coef[v1] = vertex_coef[v1] - alpha + beta * coe2 * coe1;
                    if(_isnan(vertex_coef[v1]))
                    {
                        qDebug()<<"if2 v1:"<<vertex_coef[v1]<<" "<<coe2<<" "<<coe1;
                    }
                    vertex_coef[v2] = vertex_coef[v2] + alpha + beta * coe2 * coe2;
                    if(_isnan(vertex_coef[v2]))
                    {
                        qDebug()<<"if2 v2:"<<vertex_coef[v2]<<" "<<coe2<<" "<<coe2;
                    }
                    vertex_coef[v3] = vertex_coef[v3] - alpha + beta * coe2 * coe3;
                    if(_isnan(vertex_coef[v3]))
                    {
                        qDebug()<<"if2 v3:"<<vertex_coef[v3]<<" "<<coe2<<" "<<coe3;
                    }
                    vertex_coef[v4] = vertex_coef[v4] + alpha + beta * coe2 * coe4;
                    if(_isnan(vertex_coef[v4]))
                    {
                        qDebug()<<"if2 v4:"<<vertex_coef[v4]<<" "<<coe2<<" "<<coe4;
                    }
                    right += temp_delta * beta * coe2;
                }
                else if(v3 == *v_it)
                {
                    vertex_coef[v1] = vertex_coef[v1] + alpha + beta * coe3 * coe1;
                    if(_isnan(vertex_coef[v1]))
                    {
                        qDebug()<<"if3 v1:"<<vertex_coef[v1]<<" "<<coe3<<" "<<coe1;
                    }
                    vertex_coef[v2] = vertex_coef[v2] - alpha + beta * coe3 * coe2;
                    if(_isnan(vertex_coef[v2]))
                    {
                        qDebug()<<"if3 v2:"<<vertex_coef[v2]<<" "<<coe3<<" "<<coe2;
                    }
                    vertex_coef[v3] = vertex_coef[v3] + alpha + beta * coe3 * coe3;
                    if(_isnan(vertex_coef[v3]))
                    {
                        qDebug()<<"if3 v3:"<<vertex_coef[v3]<<" "<<coe3<<" "<<coe3;
                    }
                    vertex_coef[v4] = vertex_coef[v4] - alpha + beta * coe3 * coe4;
                    if(_isnan(vertex_coef[v4]))
                    {
                        qDebug()<<"if3 v4:"<<vertex_coef[v4]<<" "<<coe3<<" "<<coe4;
                    }
                    right += temp_delta * beta * coe3;
                }
                else if(v4 == *v_it)
                {
                    vertex_coef[v1] = vertex_coef[v1] - alpha + beta * coe4 * coe1;
                    if(_isnan(vertex_coef[v1]))
                    {
                        qDebug()<<"if4 v1:"<<vertex_coef[v1]<<" "<<coe4<<" "<<coe1;
                    }
                    vertex_coef[v2] = vertex_coef[v2] + alpha + beta * coe4 * coe2;
                    if(_isnan(vertex_coef[v2]))
                    {
                        qDebug()<<"if4 v2:"<<vertex_coef[v2]<<" "<<coe4<<" "<<coe2;
                    }
                    vertex_coef[v3] = vertex_coef[v3] - alpha + beta * coe4 * coe3;
                    if(_isnan(vertex_coef[v3]))
                    {
                        qDebug()<<"if4 v3:"<<vertex_coef[v3]<<" "<<coe4<<" "<<coe3;
                    }
                    vertex_coef[v4] = vertex_coef[v4] + alpha + beta * coe4 * coe4;
                    if(_isnan(vertex_coef[v4]))
                    {
                        qDebug()<<"if4 v4:"<<vertex_coef[v4]<<" "<<coe4<<" "<<coe4;
                    }
                    right += temp_delta * beta * coe4;
                }
            }
            else
            {
                //qDebug()<<"*s_it is boundary:"<<(*s_it).idx();
            }
        }

        //在这里加上前后帧的对应位置！
        TriMesh::Point currentPt = mesh.point(*v_it);
        queryPt[0] = currentPt[0];
        queryPt[1] = currentPt[1];
        queryPt[2] = currentPt[2];
        int indexPreFrame = annSearch(kdTreePreFrame, queryPt);
        int indexNextFrame = annSearch(kdTreeNextFrame, queryPt);

        this->annDistanceTotal[0] += this->getDistance(queryPt, annCloudPreFrame[indexPreFrame]);
        this->annDistanceTotal[1] += this->getDistance(queryPt, annCloudNextFrame[indexNextFrame]);

        right_term(v_it->idx(), 0) += right[0];
        right_term(v_it->idx(), 1) += right[1];
        right_term(v_it->idx(), 2) += right[2];

        if (this->getDistance(queryPt, annCloudPreFrame[indexPreFrame]) <= this->threshold) {
            right_term(v_it->idx(), 0) += ratioPreFrame*annCloudPreFrame[indexPreFrame][0];
            right_term(v_it->idx(), 1) += ratioPreFrame*annCloudPreFrame[indexPreFrame][1];
            right_term(v_it->idx(), 2) += ratioPreFrame*annCloudPreFrame[indexPreFrame][2];
        } else {
            right_term(v_it->idx(), 0) += ratioPreFrame*queryPt[0];
            right_term(v_it->idx(), 1) += ratioPreFrame*queryPt[1];
            right_term(v_it->idx(), 2) += ratioPreFrame*queryPt[2];
        }

        if (this->getDistance(queryPt, annCloudNextFrame[indexNextFrame]) <= this->threshold) {
            right_term(v_it->idx(), 0) += ratioNextFrame*annCloudNextFrame[indexNextFrame][0];
            right_term(v_it->idx(), 1) += ratioNextFrame*annCloudNextFrame[indexNextFrame][1];
            right_term(v_it->idx(), 2) += ratioNextFrame*annCloudNextFrame[indexNextFrame][2];
        } else {
            right_term(v_it->idx(), 0) += ratioNextFrame*queryPt[0];
            right_term(v_it->idx(), 1) += ratioNextFrame*queryPt[1];
            right_term(v_it->idx(), 2) += ratioNextFrame*queryPt[2];
        }

        //right_term(v_it->idx(), 0) += right[0] + ratioPreFrame*annCloudPreFrame[indexPreFrame][0] + ratioNextFrame*annCloudNextFrame[indexNextFrame][0];
        //right_term(v_it->idx(), 1) += right[1] + ratioPreFrame*annCloudPreFrame[indexPreFrame][1] + ratioNextFrame*annCloudNextFrame[indexNextFrame][1];
        //right_term(v_it->idx(), 2) += right[2] + ratioPreFrame*annCloudPreFrame[indexPreFrame][2] + ratioNextFrame*annCloudNextFrame[indexNextFrame][2];

        for(std::map<TriMesh::VertexHandle, double>::iterator m_it = vertex_coef.begin(); m_it != vertex_coef.end(); m_it++)
        {
            /*
            if(v_it->idx()==63448||m_it->first.idx()==63448)
            {
                qDebug()<<v_it->idx()<<","<<m_it->first.idx()<<":"<<m_it->second<<endl;
            }
            */
            if(m_it->second==0.0)
            {
                qDebug()<<"zero at:"<<v_it->idx()<<","<<m_it->first.idx()<<endl;
            }
            if(_isnan(m_it->second))
            {
                qDebug()<<"nan at:"<<v_it->idx()<<","<<m_it->first.idx()<<endl;
            }
            triple.push_back(Eigen::Triplet<double>(v_it->idx(), m_it->first.idx(), m_it->second));
        }
    }
    //Fill the matrix *this with the list of triplets defined by the iterator range begin - end.
    coef_matrix.setFromTriplets(triple.begin(), triple.end());


    //print matrix size
    qDebug()<<"begin solver";
    qDebug()<<coef_matrix.rows();
    qDebug()<<coef_matrix.cols();

    Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > solver;
    solver.compute(coef_matrix);
    Eigen::MatrixXd vertices_term = solver.solve(right_term);
    if(solver.info() != Eigen::Success) {
         qDebug() << "Solve linear system failed.";
    }
    else
    {
        qDebug() << "Solve linear system success.";
    }
    for(TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
    {
        int index = v_it->idx();
        TriMesh::Point pt = TriMesh::Point(vertices_term(index,0), vertices_term(index,1), vertices_term(index,2));
        mesh.set_point(*v_it, pt);
    }
}

//把前后帧trimesh模型中的点加入annCloud
void MeshDenoisingViaL0MinimizationForFrames::initAnn()
{
    annDistanceTotal[0] = 0.0;
    annDistanceTotal[1] = 0.0;
    annCloudPreFrame = annAllocPts(meshPreFrame.n_vertices(), 3);
    qDebug() << "annCloudPreFrame size:" <<meshPreFrame.n_vertices();
    for(TriMesh::VertexIter v_it = meshPreFrame.vertices_begin(); v_it != meshPreFrame.vertices_end(); v_it++)
    {
        int index = v_it->idx();
        TriMesh::Point p = meshPreFrame.point(*v_it);
        annCloudPreFrame[index][0] = p[0];
        annCloudPreFrame[index][1] = p[1];
        annCloudPreFrame[index][2] = p[2];
    }

    annCloudNextFrame = annAllocPts(meshNextFrame.n_vertices(), 3);
    qDebug() << "annCloudNextFrame size:" <<meshNextFrame.n_vertices();
    for(TriMesh::VertexIter v_it = meshNextFrame.vertices_begin(); v_it != meshNextFrame.vertices_end(); v_it++)
    {
        int index = v_it->idx();
        TriMesh::Point p = meshNextFrame.point(*v_it);
        annCloudNextFrame[index][0] = p[0];
        annCloudNextFrame[index][1] = p[1];
        annCloudNextFrame[index][2] = p[2];
    }
    k = 1;
    idxArr = new ANNidx[k];
    distArr = new ANNdist[k];
    kdTreePreFrame = new ANNkd_tree(annCloudPreFrame, meshPreFrame.n_vertices(), 3);
    kdTreeNextFrame = new ANNkd_tree(annCloudNextFrame, meshNextFrame.n_vertices(), 3);
}

int MeshDenoisingViaL0MinimizationForFrames::annSearch(ANNkd_tree *kdTree, ANNpoint &queryPt)
{
    double eps = 0.0;

    kdTree->annkSearch(queryPt, k, idxArr, distArr, eps);

    //qDebug() << "\tNN:\tIndex\tDistance\n";
    for (int i = 0; i < k; ++i) {
        //qDebug() << "\t" << i << "\t" << idxArr[i] << "\t" << distArr[i] << "\n";
        //返回最近的节点下标
        return idxArr[i];
    }
}

double MeshDenoisingViaL0MinimizationForFrames::getDistance(ANNpoint p1, ANNpoint p2)
{
    double distance_2 = (p1[0] - p2[0])*(p1[0] - p2[0]) + (p1[1] - p2[1])*(p1[1] - p2[1]) + (p1[2] - p2[2])*(p1[2] - p2[2]);
    return sqrt(distance_2);
}

bool MeshDenoisingViaL0MinimizationForFrames::readMesh(TriMesh &mesh, const QString &fileName)
{
    try{
        if(!OpenMesh::IO::read_mesh(mesh, fileName.toStdString())) {
            qDebug()<< "read failed:"<<fileName;
            emit(statusMessage("read failed:" + fileName));
            return false;
        } else {
            qDebug()<< "read success:"<<fileName;
            return true;
        }
    }catch(QException &e){
        qDebug()<<e.what();
        emit(statusMessage("read failed:" + fileName));
        return false;
    }
}
