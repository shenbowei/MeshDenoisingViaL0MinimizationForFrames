#include "datamanager.h"

DataManager::DataManager()
{

}

bool DataManager::ImportMeshFromFile(std::string filename)
{
    //增加options 读取颜色
    OpenMesh::IO::Options opt = 0x0020;
    if(!OpenMesh::IO::read_mesh(mesh_, filename, opt))
        return false;
    else
    {
        original_mesh_ = mesh_;
        noisy_mesh_ = mesh_;
        denoised_mesh_ = mesh_;
        return true;
    }
}

bool DataManager::ExportMeshToFile(std::string filename)
{
    //增加options 导出颜色
    OpenMesh::IO::Options opt = 0x0020;
    return OpenMesh::IO::write_mesh(mesh_, filename, opt);
}
