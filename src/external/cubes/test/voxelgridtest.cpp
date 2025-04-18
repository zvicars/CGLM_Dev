#include "../VoxelGrid.hpp"
//Testing script for the voxelgrid object type to ensure that it is quantitatively accurate. Two things must remain true, the total density must add up to N/density and the 
//center of mass of identical mass particles must be conserved so long as the input parameters are reasonable
//to this end, a single atom will be added to the box in an arbitrary location and these quantities will be computed
//one caveat is that the com might be split across pbcs, so the added particle will be far away from the boundaries, this isn't ideal as an error in pbc treatment would go unnoticed
//however, backresolving the right periodic image from the grid would be ugly, so it's better to test the pbc code separately
//note that com will not be conserved perfectly if any density crosses pbc's or 2*sigma is less than any dimension
int main(){
    double testing_threshold = 1e-4;
    std::array<int,3> size = {50, 40, 45};
    std::array<double, 3> box_size = {20.0, 15.0, 10.0};
    double isovalue = 0.5;
    VoxelGrid v1(size, isovalue);
    std::array<std::size_t, 3> x_in = {10, 7, 5};
    v1.add_discrete(x_in, 1.0);
    double tot_mass = v1.getTot();
    auto tot_com = v1.getWeightedCOM();

    double diff = fabs(1 - v1.getTot());
    double diff2 = 0.0;
    
    for(int i = 0; i < 3; i++){
        diff2 += fabs(x_in[i]+0.5 - tot_com[i]); 
    }

    if(diff > testing_threshold) {
        std::cout << "Total mass is not being conserved." << std::endl;
        std::cout << "System mass = " << tot_mass << std::endl;
        std::cout << "Expected mass = 1.0" << std::endl;
        return 1;
    }
    if(diff2 > testing_threshold){
        std::cout << "Center of mass is not being conserved." << std::endl;
        for(int i = 0; i < 3; i++){
            std::cout << "Observed COM for dim " << i << " = " << tot_com[i] << std::endl;
            std::cout << "Expected COM for dim " << i << " = " << x_in[i] << std::endl;
        }
        return 1;
    }
    return 0;
}