#include "hamiltonian_CGLM.hpp"
#include "../../tools/stlmath.hpp"
#include "../lattice/lattice_pbc.hpp"
#include "../lattice/lattice.hpp"

real SubvolumePair::compute_chi(const Lattice& lattice){
	real chi2 = 0; 
  ProbeVolume* sv1 = pvs[0];
  ProbeVolume* sv2 = pvs[1];
  Vec<std::size_t> bounds1, bounds2;
  sv1->bounds(bounds1);
  sv2->bounds(bounds2);
	for(std::size_t k = bounds1[2]; k <= bounds1[5]; k++)
		for(std::size_t j = bounds1[1]; j <= bounds1[4]; j++)
			for(std::size_t i = bounds1[0]; i <= bounds1[3]; i++)
			{
        Vec3<int> ref_idx = {i,j,k};
        auto v1 = sv1->isInside({i,j,k});
				if(!v1) continue;
        std::size_t counter = 0;
				for(auto& offset : offsets_) //loop over all sites that have a non-zero chi_ij term
				{
          auto new_site = ref_idx + offset;
          auto new_idx = wrap3(new_site, lattice.size());
					bool v2 = sv2->isInside(new_idx);
					if(!v2) continue;
          bool state1 = lattice.getState({i,j,k});
					bool state2 = lattice.getState(new_idx);
          if(state1 && state2) chi2+=values_[counter];
          counter++;
				}
			}
	return chi2;
}

//this one deletes any chi values that are unused for this subvolume pair
real SubvolumePair::compute_chi_init(const Lattice& lattice){
  Vec<int> chi_counts(offsets_.size(), 0);
	real chi2 = 0; 
  ProbeVolume* sv1 = pvs[0];
  ProbeVolume* sv2 = pvs[1];
  Vec<std::size_t> bounds1, bounds2;
  sv1->bounds(bounds1);
  sv2->bounds(bounds2);
	for(std::size_t k = bounds1[2]; k <= bounds1[5]; k++)
		for(std::size_t j = bounds1[1]; j <= bounds1[4]; j++)
			for(std::size_t i = bounds1[0]; i <= bounds1[3]; i++)
			{
        Vec3<int> ref_idx = {i,j,k};
        auto v1 = sv1->isInside({i,j,k});
				if(!v1) continue;
        std::size_t counter = 0;
				for(auto& offset : offsets_) //loop over all sites that have a non-zero chi_ij term
				{
          auto new_site = ref_idx + offset;
          auto new_idx = wrap3(new_site, lattice.size());
					bool v2 = sv2->isInside(new_idx);
					if(!v2) continue;
          chi_counts[counter]++;
          bool state1 = lattice.getState({i,j,k});
					bool state2 = lattice.getState(new_idx);
          if(state1 && state2) chi2+=values_[counter];
          counter++;
				}
			}
  for(int i = 0; i < chi_counts.size(); i++){
    if(chi_counts[i] == 0) chi_counts.erase(chi_counts.begin() + i);
  }
	return chi2;
}

real SubvolumePair::compute_dchi(const Lattice& lattice){
	dchi = 0.0;
  auto active_idx = lattice.getActiveIndex();
	bool in_11, in_12, in_21, in_22;
	in_11 = pvs[0]->isInside(active_idx);
	in_12 = pvs[1]->isInside(active_idx);
	if(!(in_11 || in_12)) return 0.0; //if the changing site is not in either volume in the pair, exit and leave dchi as zero
	Vec3<int> active_site = {active_idx[0], active_idx[1], active_idx[2]};
  if(!self) // behavior for on-diagonal terms is different since it's only counted once
	{
		if(in_11)
		{
			for(unsigned int it = 1; it < offsets_.size(); it++) //iterator starts at one since chi_ij = 0, 0, 0 needs special treatment to account for the fact that the state of 1,1 and 2,1 are both being flipped
			{
        Vec3<int>& offset = offsets_[it];
        auto new_site = offset + active_site;
        auto new_idx = wrap3(new_site, lattice.size());
        bool in22 = pvs[1]->isInside(new_idx);
        if(!in22) continue;
				if(!lattice.getState(new_idx)) continue;
				dchi+=values_[it];
			}

		}
		if(in_12)
		{
			for(unsigned int it = 1; it < offsets_.size(); it++) //same motivation as above, 0, 0, 0 will be treated differently
			{
        Vec3<int>& offset = offsets_[it];
        auto new_site = offset + active_site;
        auto new_idx = wrap3(new_site, lattice.size());
        bool in21 = pvs[0]->isInside(new_idx);
				if(!in21) continue;
        if(!lattice.getState(new_idx))
				dchi += values_[it];
			}
		}
		if(in_11 && in_12) dchi += values_[0];
	}
	else
	{
		for(unsigned int it = 1; it < offsets_.size(); it++) //iterator starts at one since chi_ij = 0, 0, 0 needs special treatment to account for the fact that the state of 1,1 and 2,1 are both being flipped
		{
      Vec3<int>& offset = offsets_[it];
      auto new_site = offset + active_site;
      auto new_idx = wrap3(new_site, lattice.size());
      auto new_state = lattice.getState(new_idx);
			in_21 = pvs[0]->isInside(new_idx);	
			if(lattice.getState(new_idx)) dchi+=2.0*values_[it];
		}
		dchi += values_[0];
	}
	
	if(lattice.getActiveState()) dchi = -dchi;
	return dchi;
}

