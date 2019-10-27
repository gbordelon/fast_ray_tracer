#ifndef PHOTON_TRACER
#define PHOTON_TRACER

void trace_photons(const World w, const size_t num_photons, bool, bool);

PhotonMap *array_of_photon_maps(size_t num);

#endif
