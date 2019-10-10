#ifndef PM_H
#define PM_H

/* This is the photon.  The power is not compressed so the size is 28 bytes */
typedef struct Photon {
  double pos[3];                 // photon position
  short plane;                  // splitting plane for kd-tree
  unsigned char theta, phi;     // incoming direction
  double power[3];               // photon power (uncompressed)
} Photon;


/* This structure is used only to locate the nearest photons */
typedef struct NearestPhotons {
  long max;
  long found;
  int got_heap;
  double pos[3];
  double *dist2;
  Photon **index;
} NearestPhotons;

typedef struct {

  Photon *photons;

  long stored_photons;
  long half_stored_photons;
  long max_photons;
  long prev_scale;

  double costheta[256];
  double sintheta[256];
  double cosphi[256];
  double sinphi[256];

  double bbox_min[3];            // use bbox_min;
  double bbox_max[3];            // use bbox_max;
} PhotonMap;

void init_Photon_map( long max_phot, PhotonMap *pm );

void delete_Photon_map(PhotonMap *pm);

  void pm_store(
          PhotonMap *pm,
          double power[3],          // photon power
          double pos[3],            // photon position
          double dir[3] );          // photon direction

  void pm_scale_photon_power(
          PhotonMap *pm,
          double scale );           // 1/(number of emitted photons)

  void pm_balance(PhotonMap *pm);           // balance the kd-tree (before use!)

  void pm_irradiance_estimate(
    PhotonMap *pm,
    double irrad[3],                // returned irradiance
    double pos[3],                  // surface position
    double normal[3],               // surface normal at pos
    double max_dist,                // max distance to look for photons
    int nphotons,
    double cone_filter_k );                // number of photons to use

  void pm_locate_photons(
    PhotonMap *pm,
    NearestPhotons *np,            // np is used to locate the photons
    long index );                   // call with index = 1

  void pm_photon_dir(
    PhotonMap *pm,
    double *dir,                    // direction of photon (returned)
    Photon *p );                   // the photon

#endif
