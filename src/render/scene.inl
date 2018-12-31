// ========================================================================== //
// Copyright (c) 2017-2018 The University of Texas at Austin.                 //
// All rights reserved.                                                       //
//                                                                            //
// Licensed under the Apache License, Version 2.0 (the "License");            //
// you may not use this file except in compliance with the License.           //
// A copy of the License is included with this software in the file LICENSE.  //
// If your copy does not contain the License, you may obtain a copy of the    //
// License at:                                                                //
//                                                                            //
//     https://www.apache.org/licenses/LICENSE-2.0                            //
//                                                                            //
// Unless required by applicable law or agreed to in writing, software        //
// distributed under the License is distributed on an "AS IS" BASIS, WITHOUT  //
// WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.           //
// See the License for the specific language governing permissions and        //
// limitations under the License.                                             //
//                                                                            //
// ========================================================================== //

#if !defined(SPRAY_SCENE_INL_)
#error An implementation of Scene
#endif

#define DEBUG_SCENE
#undef DEBUG_SCENE

namespace spray {

template <typename CacheT, typename SurfaceBufT>
void Scene<CacheT, SurfaceBufT>::init(const std::string& desc_filename,
                                      const std::string& ply_path,
                                      const std::string& storage_basepath,
                                      int cache_size, int view_mode,
                                      bool insitu_mode, int num_partitions,
                                      int num_light_samples) {
  // load .domain file
  SceneLoader loader;
  loader.load(desc_filename, ply_path, num_light_samples, &domains_, &lights_);

  // merge domain bounds and find the scene bounds
  std::size_t max_num_vertices, max_num_faces;
  mergeDomainBounds(&max_num_vertices, &max_num_faces);

  // TODO: support ooc (ooc not supported at this moment)
  // see notes in TriMeshBuffer::loadShapes
  for (std::size_t i = 0; i < domains_.size(); ++i) {
    if (!domains_[i].shapes.empty()) {
      CHECK_LT(cache_size, 0);
    }
  }

  insitu_ = false;

  if (view_mode == VIEW_MODE_PARTITION) {
    partition_num_ = 0;
    num_partitions_ = num_partitions;

    // partition data if in-situ mode
    partition_.partition(getNumDomains(), getDomains(), getBound(),
                         num_partitions);

  } else if (insitu_mode) {
    insitu_ = true;

    // partition data if in-situ mode
    partition_.partition(getNumDomains(), getDomains(), getBound(),
                         mpi::size());
#ifdef SPRAY_GLOG_CHECK
    LOG(INFO) << "<<<<<INSITU MODE>>>>>>";
    for (auto& r : partition_.getDomains(mpi::rank())) {
      LOG(INFO) << "rank " << mpi::rank() << " domain " << r;
    }
#endif
    // NOTE: override cache size
    cache_size = partition_.getNumDomains(mpi::rank());
  }

  // copy to local disk
  if (!storage_basepath.empty()) {
    storage_basepath_ = storage_basepath;

    copyAllDomainsToLocalDisk(storage_basepath, insitu_mode);
  }

  // initialize cache
  if (!(view_mode == VIEW_MODE_DOMAIN || view_mode == VIEW_MODE_PARTITION)) {
    cache_.init(domains_.size(), cache_size, insitu_mode);

    // initialize mesh buffer
    surface_buf_.init(cache_.getCacheSize(), max_num_vertices, max_num_faces,
                      true /* compute_normals */);

    // warm up cache
    if (view_mode == VIEW_MODE_FILM || view_mode == VIEW_MODE_GLFW) {
      if (insitu_mode) {
        const std::list<int>& domains = partition_.getDomains(mpi::rank());
        for (int id : domains) load(id);
      } else if (cache_size < 0) {
        for (std::size_t id = 0; id < domains_.size(); ++id) load(id);
      }
    }
  }

  wbvh_.init(getBound(), getDomains());

  // for domain vis.
  glfw_domain_idx_ = 0;
}

template <typename CacheT, typename SurfaceBufT>
void Scene<CacheT, SurfaceBufT>::buildWbvh() {
#if defined(SPRAY_ISECT_PACKET1)
  wbvh_.build(WbvhEmbree::NORMAL1);

#elif defined(SPRAY_ISECT_PACKET8)

#if !defined(SPRAY_AVX2)
#warning Use AVX2 for a packet of 8 rays.
#endif

  wbvh_.build(WbvhEmbree::NORMAL8);

#elif defined(SPRAY_ISECT_PACKET16)

#if !defined(SPRAY_AVX512)
#error unsupported
#warning Use AVX512 for a packet of 16 rays.
#endif

  wbvh_.build(WbvhEmbree::NORMAL16);

#elif defined(SPRAY_ISECT_STREAM_1M)

  wbvh_.build(WbvhEmbree::STREAM_1M);

#else
#error unsupported

#endif
}

template <typename CacheT, typename SurfaceBufT>
void Scene<CacheT, SurfaceBufT>::load(int id) {
  int cache_block;
  if (cache_.load(id, &cache_block)) {
#ifdef DEBUG_SCENE
    LOG(INFO) << "loading cached domain " << id << " cache block "
              << cache_block << " $size " << cache_.getSize() << " $capacity "
              << cache_.getCacheSize();
#endif
    scene_ = surface_buf_.get(cache_block);
  } else {
#ifdef DEBUG_SCENE
    LOG(INFO) << "loading uncached domain " << id << " cache block "
              << cache_block << " $size " << cache_.getSize() << " $capacity "
              << cache_.getCacheSize();
#endif
    const glm::mat4& x = domains_[id].transform;
    bool apply_transform = (x != glm::mat4(1.0));

    scene_ = surface_buf_.load(domains_[id].filename, cache_block, x,
                               apply_transform, domains_[id].shapes);

    // cache_.setLoaded(cache_block);
  }
  cache_block_ = cache_block;
}

template <typename CacheT, typename SurfaceBufT>
void Scene<CacheT, SurfaceBufT>::load(int id, SceneInfo* sinfo) {
  int cache_block;
  if (cache_.load(id, &cache_block)) {
#ifdef DEBUG_SCENE
    LOG(INFO) << "loading cached domain " << id << " cache block "
              << cache_block << " $size " << cache_.getSize() << " $capacity "
              << cache_.getCacheSize();
#endif
    scene_ = surface_buf_.get(cache_block);
  } else {
#ifdef DEBUG_SCENE
    LOG(INFO) << "loading uncached domain " << id << " cache block "
              << cache_block << " $size " << cache_.getSize() << " $capacity "
              << cache_.getCacheSize();
#endif
    const glm::mat4& x = domains_[id].transform;
    bool apply_transform = (x != glm::mat4(1.0));

    scene_ = surface_buf_.load(domains_[id].filename, cache_block, x,
                               apply_transform, domains_[id].shapes);

    // cache_.setLoaded(cache_block);
  }
  sinfo->rtc_scene = scene_;
  sinfo->cache_block = cache_block;
}

template <typename CacheT, typename SurfaceBufT>
bool Scene<CacheT, SurfaceBufT>::intersect(RTCScene rtc_scene, int cache_block,
                                           RTCRayIntersection* isect) {
  rtcIntersect(rtc_scene, (RTCRay&)(*isect));

  if (isect->geomID != RTC_INVALID_GEOMETRY_ID) {
    surface_buf_.updateIntersection(cache_block, isect);
    return true;
  }
  return false;
}

template <typename CacheT, typename SurfaceBufT>
bool Scene<CacheT, SurfaceBufT>::occluded(RTCScene rtc_scene, RTCRay* ray) {
  rtcOccluded(rtc_scene, *ray);
  if (ray->geomID != RTC_INVALID_GEOMETRY_ID) {  // occluded
    return true;
  }
  return false;  // unoccluded
}

// copy only those domains mapped to this process.
// we don't have to copy everything in some cases.
template <typename CacheT, typename SurfaceBufT>
void Scene<CacheT, SurfaceBufT>::copyAllDomainsToLocalDisk(
    const std::string& dest_path, bool insitu_mode) {
  // clean up existing folder

  // std::string stuff_to_remove = dest_path + "/*";
  // std::string cmd_cleanup = "rm -rf " + stuff_to_remove;
  // std::cout << "cleaning up exiting files in " + dest_path << "\n";
  // std::system(cmd_cleanup.c_str());

  std::vector<int> ids;

  std::size_t domain_size;
  if (insitu_mode) {
    const std::list<int>& domains = partition_.getDomains(mpi::rank());
    ids.reserve(domains.size());
    for (int id : domains) {
      ids.push_back(id);
    }
  } else {
    ids.reserve(domains_.size());
    for (std::size_t i = 0; i < domains_.size(); ++i) {
      ids.push_back(i);
    }
  }

  int res;
  // for (std::size_t i = 0; i < domains_.size(); ++i) {
  for (std::size_t i = 0; i < ids.size(); ++i) {
    int id = ids[i];
    Domain& domain = domains_[id];
    CHECK_EQ(id, domain.id);

    // create a new folder
    std::string new_dir = dest_path + "/proc" + std::to_string(mpi::rank()) +
                          "_domain" + std::to_string(id);
    std::string cmd_make_dir = "mkdir " + new_dir;

    res = std::system(cmd_make_dir.c_str());

#ifdef SPRAY_GLOG_CHECK
    LOG(INFO) << cmd_make_dir << " " << res;
#endif
    std::cout << cmd_make_dir << " " << res << std::endl;
    CHECK_EQ(res, 0);

    // extract basename
    std::string bname = std::string(util::getFilename(domain.filename.c_str()));
    std::string destination_file = new_dir + "/" + bname;

    // cp model file to local disk
    std::string cmd_copy_domain =
        "cp " + domain.filename + " " + destination_file;

    int res = std::system(cmd_copy_domain.c_str());

#ifdef SPRAY_GLOG_CHECK
    LOG(INFO) << cmd_copy_domain << " " << res;
#endif
    std::cout << cmd_copy_domain << " " << res << std::endl;
    CHECK_EQ(res, 0);

    // update the descriptor so it now points to the copied model file
    domain.filename = destination_file;
  }
}

template <typename CacheT, typename SurfaceBufT>
void Scene<CacheT, SurfaceBufT>::deleteAllDomainsFromLocalDisk() {
  CHECK_EQ(storage_basepath_.empty(), false);
  int res;
  for (std::size_t i = 0; i < domains_.size(); ++i) {
    // create a new folder
    std::string new_dir = storage_basepath_ + "/proc" +
                          std::to_string(mpi::rank()) + "_domain" +
                          std::to_string(i);
    std::string cmd_rm_dir = "rm -rf " + new_dir;

    res = std::system(cmd_rm_dir.c_str());
  }
}

template <typename CacheT, typename SurfaceBufT>
void Scene<CacheT, SurfaceBufT>::mergeDomainBounds(
    std::size_t* max_num_vertices, std::size_t* max_num_faces) {
  Aabb world_space_bound;

  std::size_t num_domains = domains_.size();
  domains_.resize(num_domains);

#if defined(PRINT_DOMAIN_BOUNDS) && defined(SPRAY_GLOG_CHECK)
  std::size_t total_faces = 0;
#endif

  // evaluate bound and number of primitives
  std::size_t num_vertices = 0, num_faces = 0;

  for (std::size_t id = 0; id < num_domains; ++id) {
    Domain& d = domains_[id];

    // let's enforce that that domain world-space bound and
    // the number of faces are provided through preprocessing.
    CHECK_EQ(d.world_aabb.isValid(), true);
    if (d.shapes.empty()) {
      CHECK_GT(d.num_vertices, 0);
      CHECK_GT(d.num_faces, 0);
    }

    // maximum values
    if (d.num_vertices > num_vertices) num_vertices = d.num_vertices;
    if (d.num_faces > num_faces) num_faces = d.num_faces;

#if defined(PRINT_DOMAIN_BOUNDS) && defined(SPRAY_GLOG_CHECK)
    if (mpi::isRootProcess()) {
      total_faces += d.num_faces;
      LOG(INFO) << "[domain " << id << "] [bounds " << d.world_aabb
                << "] [faces " << d.num_faces << "]";
    }
#endif

    world_space_bound.merge(d.world_aabb);

#if defined(PRINT_DOMAIN_BOUNDS) && defined(SPRAY_GLOG_CHECK)
    LOG(INFO) << " [Scene::load()] [domain " << domains_[id].id
              << "] bound: " << domains_[id].world_aabb;
#endif
  }

  *max_num_vertices = num_vertices;
  *max_num_faces = num_faces;

  bound_ = world_space_bound;
  CHECK(bound_.isValid());

#if defined(PRINT_DOMAIN_BOUNDS) && defined(SPRAY_GLOG_CHECK)
  LOG_IF(INFO, mpi::isRootProcess()) << "total faces: " << total_faces;
  LOG_IF(INFO, mpi::isRootProcess()) << "world bound: " << bound_;
#endif
}

template <typename CacheT, typename SurfaceBufT>
void Scene<CacheT, SurfaceBufT>::drawDomains() {
  glPushAttrib(GL_ENABLE_BIT | GL_LINE_BIT);
  glDisable(GL_LIGHTING);
  glLineWidth(.1);
  glEnable(GL_DEPTH_TEST);

  glPolygonOffset(1.0, 1.0);

  glDepthMask(GL_FALSE);

  glm::vec4 color(0.3f, 0.3f, 0.3f, 0.5f);
  glm::vec4 select(0.0f, 1.0f, 1.0f, .5f);

  for (std::size_t i = 0; i < domains_.size(); ++i) {
    const Domain& d = domains_[i];
    if (i != glfw_domain_idx_) {
      d.world_aabb.draw(color);
    }
  }

  glLineWidth(2.0);
  for (std::size_t i = 0; i < domains_.size(); ++i) {
    const Domain& d = domains_[i];
    if (i == glfw_domain_idx_) {
      d.world_aabb.draw(select);
    }
  }

  glDepthMask(GL_TRUE);
  glPopAttrib();
}

template <typename CacheT, typename SurfaceBufT>
void Scene<CacheT, SurfaceBufT>::drawPartitions() {
  glPushAttrib(GL_ENABLE_BIT | GL_LINE_BIT);
  glDisable(GL_LIGHTING);
  glLineWidth(.1);
  glEnable(GL_DEPTH_TEST);

  glPolygonOffset(1.0, 1.0);

  glDepthMask(GL_FALSE);

  glm::vec4 color(0.3f, 0.3f, 0.3f, 0.5f);
  glm::vec4 select(0.0f, 1.0f, 1.0f, .5f);

  std::vector<int> domain_to_partition;
  domain_to_partition = partition_.computePartition(num_partitions_);

  for (std::size_t i = 0; i < domains_.size(); ++i) {
    const Domain& d = domains_[i];
    if (domain_to_partition[d.id] != partition_num_) {
      d.world_aabb.draw(color);
    }
  }
  glLineWidth(2.0);
  for (std::size_t i = 0; i < domains_.size(); ++i) {
    const Domain& d = domains_[i];
    if (domain_to_partition[d.id] == partition_num_) {
      d.world_aabb.draw(select);
    }
  }

  glDepthMask(GL_TRUE);
  glPopAttrib();
}

}  // namespace spray
