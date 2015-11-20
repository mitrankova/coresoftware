#ifndef __PHG4TRUTHINFOCONTAINER_H__
#define __PHG4TRUTHINFOCONTAINER_H__

#include <phool/PHObject.h>
#include <map>
#include <set>

class PHG4Particle;
class PHG4VtxPoint;

class PHG4TruthInfoContainer: public PHObject {
  
public:

  typedef std::map<int,PHG4Particle *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  typedef std::map<int,PHG4VtxPoint *> VtxMap;
  typedef VtxMap::iterator VtxIterator;
  typedef VtxMap::const_iterator ConstVtxIterator;
  typedef std::pair<VtxIterator, VtxIterator> VtxRange;
  typedef std::pair<ConstVtxIterator, ConstVtxIterator> ConstVtxRange;

  PHG4TruthInfoContainer();
  virtual ~PHG4TruthInfoContainer();

  void Reset();
  void identify(std::ostream& os = std::cout) const;

  // --- particle storage ------------------------------------------------------
 
  //! Add a particle that the user has created
  ConstIterator AddParticle(const int particleid, PHG4Particle* newparticle);

  PHG4Particle* GetParticle(const int particleid);
  PHG4Particle* GetPrimaryParticle(const int particleid);

  ConstIterator AddPrimaryParticle(PHG4Particle *newparticle);
  
  //! Get a range of iterators covering the entire container
  Range GetParticleRange();
  ConstRange GetParticleRange() const;

  //! particle size
  unsigned int size( void ) const {return particlemap.size();}
  int GetNumPrimaryVertexParticles() {return primary_particle_map.size();}
  
  //! Get the map itself
  const Map& GetMap() const {return particlemap;}
  const Map& GetPrimaryMap() const {return primary_particle_map;}

  int maxtrkindex() const;
  int mintrkindex() const;

  int maxprimarytrkindex() const;
  int minprimarytrkindex() const;

  void delete_particle(Iterator piter);

  int GetLastParticleIndex() {return (primary_particle_map.rbegin())->first;}

  std::pair< std::set<int>::const_iterator, std::set<int>::const_iterator > GetEmbeddedTrkIds() const
    {return std::make_pair(embedded_trkid.begin(), embedded_trkid.end());}
  void AddEmbededTrkId(const int i) {embedded_trkid.insert(i);}

  int isEmbeded(const int trackid) const;
  
  // --- vertex storage --------------------------------------------------------
  
  //! Add a vertex and return an iterator to the user
  VtxIterator AddVertex(const int vtxid);
  //! Add a vertex and return an iterator to the user
  ConstVtxIterator AddVertex(const int vtxid, PHG4VtxPoint* vertex);

  //! Add a primary vertex and return index to the user
  int AddPrimaryVertex(PHG4VtxPoint *);
  
  PHG4VtxPoint* GetVtx(const int vtxid);
  PHG4VtxPoint* GetPrimaryVtx(const int vtxid);
  
  //! Get a range of iterators covering the entire vertex container
  VtxRange GetVtxRange();
  ConstVtxRange GetVtxRange() const;

  //! Get the number of vertices stored
  unsigned int GetNumVertices() const {return vtxmap.size();}

  const VtxMap& GetVtxMap() const {return vtxmap;}

  int maxvtxindex() const;
  int minvtxindex() const;
 
  void delete_vtx(VtxIterator viter);

  // returns the first primary vertex that was processed by Geant4
  int GetPrimaryVertexIndex() {return (vtxmap.lower_bound(1))->first;}

  std::set<int> GetSubEventIds() const;
  Range         GetSubEventPrimaryParticleRange(int subevent);
  Range         GetSubEventSecondaryParticleRange(int subevent);
  VtxRange      GetSubEventPrimaryVertexRange(int subevent);
  VtxRange      GetSubEventSecondaryVertexRange(int subevent); 
  void          MarkSubEventBoundary();
  
  // deprecated interface, confusingly named as we store particles not hits ----
  // do not call these functions in new code, i'm leaving these for now for
  // build compatibility outside of coresoftware
  
  ConstIterator AddHit(const int particleid, PHG4Particle *newparticle) {return AddParticle(particleid,newparticle);}
  PHG4Particle* GetHit(const int particleid) {return GetParticle(particleid);}
  PHG4Particle* GetPrimaryHit(const int particleid) {return GetPrimaryParticle(particleid);}
  Range GetHitRange() {return GetParticleRange();}
  ConstRange GetHitRange() const {return GetParticleRange();}
  void delete_hit(Iterator piter) {delete_particle(piter);}

  // end deprecated interface, confusingly named -------------------------------
  
 private:

  // map format description:
  // primary particles are appended in the positive direction
  // secondary particles are appended in the negative direction
  // subevent boundaries between geant runs are stored in the subevent markers

  // +N   primary particle id => particle*
  // +N-1 
  // ...
  // +J+1 subevent #2
  // +J   subevent #1 boundary (last particle flagged in particle_subevents)
  // +J-1 subevent #1
  // ..
  // +1   primary particle id => particle*
  // 0    no entry
  // -1   secondary particle id => particle*
  // ...
  // -K+1 subevent #1
  // -K   subevent #1 boundary (last particle flagged in particle_subevents)
  // -K-1 subevent #2
  // ..
  // -M+1
  // -M   secondary particle id => particle*
  
  Map particlemap;
  VtxMap vtxmap;

  std::map< int, std::pair<int,int> > particle_subevents; // subevent index => lower key and upper key
  std::map< int, std::pair<int,int> > vertex_subevents;   // subevent index => lower key and upper key
  std::map< int, int> particle_embed_flags;               // trackid => embed flag
  std::map< int, int> vertex_embed_flags;                 // trackid => embed flag

  // --- deprecated storage ----------
  Map primary_particle_map;
  VtxMap primary_vtxmap;
  // track ids of embedded particles
  std::set<int> embedded_trkid;

  ClassDef(PHG4TruthInfoContainer,1)
};

#endif
