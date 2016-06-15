#ifndef LIBRARY_ACCESS_H
#define LIBRARY_ACCESS_H

#include <string>
#include <vector>

class LibraryAccess{

  public:
    void LoadLibraryFromFile(std::string libraryfile, bool reflected, bool reflT0);
    const float* GetReflT0(size_t Voxel, int no_pmt);
    const float* GetReflCounts(size_t Voxel, int no_pmt);
    const float* GetCounts(size_t Voxel, int no_pmt);
    const float* GetLibraryEntries(int VoxID, bool wantReflected, int no_pmt);
    std::vector<int> GetVoxelCoords(int id, double position[3]);
    std::vector<double> PhotonLibraryAnalyzer(double _energy, const int _scint_yield, const double _quantum_efficiency, int _pmt_number, int _rand_voxel);

    LibraryAccess();

  private:
    std::vector<std::vector<float> > table_;
    std::vector<std::vector<float> > reflected_table_;
    std::vector<std::vector<float> > reflT_table_;

    const double gLowerCorner[3] = {2.5, -200, 0};
    const double gUpperCorner[3] = {202.5, 200, 500};
    const int gxSteps = 40;
    const int gySteps = 80;
    const int gzSteps = 100;



};



#endif
