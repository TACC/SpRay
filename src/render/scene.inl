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

namespace spray {

template <typename CacheT, typename SurfaceBufT>
void Scene<CacheT, SurfaceBufT>::getAssignedGeometrySizes(
    std::size_t* max_num_vertices, std::size_t* max_num_faces) {
  const auto& ids = partition_.getDomains(mpi::rank());

  std::size_t nvertices = 0;
  std::size_t nfaces = 0;

  for (auto id : ids) {
    const Domain& domain = domains_[id];
    if (domain.getNumVertices() > nvertices) {
      nvertices = domain.getNumVertices();
    }
    if (domain.getNumFaces() > nfaces) {
      nfaces = domain.getNumFaces();
    }
  }

  *max_num_vertices = nvertices;
  *max_num_faces = nfaces;
}

template <typename CacheT, typename SurfaceBufT>
void Scene<CacheT, SurfaceBufT>::init(const Config& cfg) {
  const std::string& storage_basepath = cfg.local_disk_path;
  int cache_size = cfg.cache_size;
  int view_mode = cfg.view_mode;
  bool insitu_mode = (cfg.partition == spray::Config::INSITU);

  int num_partitions = cfg.num_partitions;
  int num_light_samples = cfg.ao_samples;

  // load scene file
  SceneLoader loader;
  loader.load(cfg.model_descriptor_filename, cfg.ply_path, num_light_samples,
              &domains_, &lights_);

  // populate domain info and merge scene bounds
  std::size_t max_num_vertices, max_num_faces;
  loadAndPopulateDomainInfo(&max_num_vertices, &max_num_faces);

  // TODO: support ooc (ooc not supported at this moment)
  // see notes in TriMeshBuffer::loadShapes
  for (std::size_t i = 0; i < domains_.size(); ++i) {
    if (domains_[i].hasShapes()) {
      CHECK_LT(cache_size, 0);
    }
  }

  insitu_ = false;

  if (view_mode == VIEW_MODE_PARTITION) {
    partition_num_ = 0;
    num_partitions_ = num_partitions;

    // partition data if in-situ mode
    partition_.partition(getNumDomains(), getDomains(), world_aabb_,
                         num_partitions);

  } else if (insitu_mode) {
    insitu_ = true;

    // partition data if in-situ mode
    partition_.partition(getNumDomains(), getDomains(), world_aabb_,
                         mpi::size());
#ifdef SPRAY_GLOG_CHECK
    LOG(INFO) << "<<<<<INSITU MODE>>>>>>";
    for (auto& r : partition_.getDomains(mpi::rank())) {
      LOG(INFO) << "rank " << mpi::rank() << " domain " << r;
    }
#endif
  }

  // TODO: enable this
  /*
    // copy to local disk
    if (!storage_basepath.empty()) {
      storage_basepath_ = storage_basepath;

      copyAllDomainsToLocalDisk(storage_basepath, insitu_mode);
    }
  */

  // initialize cache
  if (!(view_mode == VIEW_MODE_DOMAIN || view_mode == VIEW_MODE_PARTITION)) {
    std::size_t sf_cache_size;
    if (insitu_mode) {
      sf_cache_size = partition_.getNumDomains(mpi::rank());

      // override max_num_vertices and max_num_faces
      // max number of vertices/faces assigned to this task
      getAssignedGeometrySizes(&max_num_vertices, &max_num_faces);

    } else if (cache_size < 0 || cache_size > domains_.size()) {
      sf_cache_size = domains_.size();

    } else {
      sf_cache_size = cache_size;
    }

    cache_.init(domains_.size(),
                sf_cache_size /*unused for InsituCache and InfiniteCache*/);

    surface_buf_.init(cfg.use_spray_color, sf_cache_size, max_num_vertices,
                      max_num_faces);

    // warm up cache
    if (view_mode == VIEW_MODE_FILM || view_mode == VIEW_MODE_GLFW) {
      if (insitu_mode) {
        const std::list<int>& domains = partition_.getDomains(mpi::rank());
        for (int id : domains) {
          load(id);
        }
      } else if (cache_size < 0) {
        for (std::size_t id = 0; id < domains_.size(); ++id) {
          load(id);
        }
      }
    }
  }

  wbvh_.init(world_aabb_, getDomains());

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
    scene_ = surface_buf_.get(cache_block);
  } else {
    scene_ = surface_buf_.load(cache_block, domains_[id]);
  }
  cache_block_ = cache_block;
}

template <typename CacheT, typename SurfaceBufT>
void Scene<CacheT, SurfaceBufT>::load(int id, SceneInfo* sinfo) {
  int cache_block;
  if (cache_.load(id, &cache_block)) {
    scene_ = surface_buf_.get(cache_block);
  } else {
    scene_ = surface_buf_.load(cache_block, domains_[id]);
  }
  sinfo->rtc_scene = scene_;
  sinfo->cache_block = cache_block;
}

template <typename CacheT, typename SurfaceBufT>
bool Scene<CacheT, SurfaceBufT>::intersect(RTCScene rtc_scene, int cache_block,
                                           RTCRayIntersection* isect) const {
  rtcIntersect(rtc_scene, (RTCRay&)(*isect));

  if (isect->geomID != RTC_INVALID_GEOMETRY_ID) {
    surface_buf_.updateIntersection(cache_block, isect);
    return true;
  }
  return false;
}

template <typename CacheT, typename SurfaceBufT>
bool Scene<CacheT, SurfaceBufT>::occluded(RTCScene rtc_scene,
                                          RTCRay* ray) const {
  rtcOccluded(rtc_scene, *ray);
  if (ray->geomID != RTC_INVALID_GEOMETRY_ID) {  // occluded
    return true;
  }
  return false;  // unoccluded
}

/*
// TODO: enable this, move this to scene loader
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
    CHECK_EQ(id, domain.getId());

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

    for (auto& model : domain.getModels()) {
      // extract basename
      std::string bname =
          std::string(util::getFilename(model.getFilename().c_str()));
      std::string destination_file = new_dir + "/" + bname;

      // cp model file to local disk
      std::string cmd_copy_domain =
          "cp " + model.getFilename() + " " + destination_file;

      int res = std::system(cmd_copy_domain.c_str());

#ifdef SPRAY_GLOG_CHECK
      LOG(INFO) << cmd_copy_domain << " " << res;
#endif
      std::cout << cmd_copy_domain << " " << res << std::endl;
      CHECK_EQ(res, 0);

      // update the descriptor so it now points to the copied model file
      model.setFilename(destination_file);
    }
  }
}
*/

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
      d.getWorldAabb().draw(color);
    }
  }

  glLineWidth(2.0);
  for (std::size_t i = 0; i < domains_.size(); ++i) {
    const Domain& d = domains_[i];
    if (i == glfw_domain_idx_) {
      d.getWorldAabb().draw(select);
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
    if (domain_to_partition[d.getId()] != partition_num_) {
      d.getWorldAabb().draw(color);
    }
  }
  glLineWidth(2.0);
  for (std::size_t i = 0; i < domains_.size(); ++i) {
    const Domain& d = domains_[i];
    if (domain_to_partition[d.getId()] == partition_num_) {
      d.getWorldAabb().draw(select);
    }
  }

  glDepthMask(GL_TRUE);
  glPopAttrib();
}

template <typename CacheT, typename SurfaceBufT>
void Scene<CacheT, SurfaceBufT>::getModelInfo(
    std::size_t total_num_models, std::vector<ModelInfo>* model_info_recvbuf) {
  std::vector<ModelInfo> model_info_sendbuf;

  // assign models to mpi tasks
  std::size_t num_assigned_models = total_num_models / mpi::size();

  // if no models are assigned, just assume it has one model assigned for a
  // place hodler.
  if (num_assigned_models == 0) num_assigned_models = 1;

  model_info_sendbuf.resize(num_assigned_models * mpi::size());
  model_info_recvbuf->resize(num_assigned_models * mpi::size());

  std::size_t begin = mpi::rank() * num_assigned_models;

  // read model file if my assigned models are in the valid range
  if (begin < total_num_models) {
    std::size_t end = begin + num_assigned_models;
    if (end > total_num_models) end = total_num_models;

    PlyLoader::LongHeader long_header;
    PlyLoader::Header header;
    std::string extension;

    std::size_t model_id = 0;

    for (std::size_t id = 0; id < domains_.size(); ++id) {
      const Domain& domain = domains_[id];
      std::size_t num_models = domain.getNumModels();

      for (std::size_t m = 0; m < num_models; ++m) {
        // load assigned models only
        if (model_id >= begin && model_id < end) {
          const SurfaceModel& model = domain.getModel(m);
          const std::string& filename = model.getFilename();
          CHECK_EQ(filename.empty(), false);
          extension = spray::util::getFileExtension(filename);
          CHECK_EQ(extension, "ply");

          bool model_valid_vertex_count = model.isNumVerticesSet();
          bool model_valid_face_count = model.isNumFacesSet();
          bool model_valid_aabb = model.isObjectAabbSet();
          bool model_configured = model_valid_vertex_count &&
                                  model_valid_face_count && model_valid_aabb;

          if (!model_configured) {
            ModelInfo* info = &model_info_sendbuf[model_id];

            if (!model_valid_aabb) {
              PlyLoader::readLongHeader(filename, &long_header);

              info->num_vertices = long_header.num_vertices;
              info->num_faces = long_header.num_faces;

              for (int i = 0; i < 3; ++i) {
                info->obj_bounds_min[i] = long_header.bounds.getMin()[i];
              }
              for (int i = 0; i < 3; ++i) {
                info->obj_bounds_max[i] = long_header.bounds.getMax()[i];
              }
            } else {
              PlyLoader::quickHeaderRead(filename, &header);

              info->num_vertices = header.num_vertices;
              info->num_faces = header.num_faces;
            }
          }
        }
        ++model_id;
      }
    }
  }

  // gather model info
  std::size_t count = sizeof(ModelInfo) * num_assigned_models;

  MPI_Allgather(&model_info_sendbuf[begin], count, MPI_UNSIGNED_CHAR,
                &(*model_info_recvbuf)[0], count, MPI_UNSIGNED_CHAR,
                MPI_COMM_WORLD);
}

template <typename CacheT, typename SurfaceBufT>
void Scene<CacheT, SurfaceBufT>::updateDomains(
    bool all_domains_configured,
    const std::vector<ModelInfo>& model_info_recvbuf,
    std::size_t* max_num_vertices, std::size_t* max_num_faces) {
  std::size_t model_id = 0;
  std::size_t max_nvertices = 0, max_nfaces = 0;
  std::size_t nvertices = 0, nfaces = 0;
  std::size_t total_num_faces = 0;

  world_aabb_.reset();

  for (std::size_t id = 0; id < domains_.size(); ++id) {
    Domain* domain = &domains_[id];

    if (!all_domains_configured) {
      std::size_t num_models = domain->getNumModels();

      for (std::size_t m = 0; m < num_models; ++m) {
        const SurfaceModel& model = domain->getModel(m);
        const ModelInfo& info = model_info_recvbuf[model_id];

        if (!model.isNumVerticesSet()) {
          domain->setNumVertices(m, info.num_vertices);
        }
        if (!model.isNumFacesSet()) {
          domain->setNumFaces(m, info.num_faces);
        }
        if (!model.isObjectAabbSet()) {
          domain->setObjectAabb(m, info.obj_bounds_min, info.obj_bounds_max);
        }
        ++model_id;
      }
    }

    domain->updateDomainInfo();

    CHECK_EQ(domain->getWorldAabb().isValid(), true)
        << id << "," << domains_.size();

    if (!domain->hasShapes()) {
      CHECK_GT(domain->getNumVertices(), 0);
      CHECK_GT(domain->getNumFaces(), 0);
    }

    // maximum values

    if (domain->getNumVertices() > max_nvertices) {
      max_nvertices = domain->getNumVertices();
    }

    if (domain->getNumFaces() > max_nfaces) {
      max_nfaces = domain->getNumFaces();
    }

    if (mpi::rank() == 0) {
      total_num_faces += domain->getNumFaces();
      std::cout << "[INFO] [domain " << id << "] [bounds "
                << domain->getWorldAabb() << "] [face count "
                << domain->getNumFaces() << "]" << std::endl;
    }

    world_aabb_.merge(domain->getWorldAabb());
  }  // for (std::size_t id = 0; id < num_domains; ++id) {

  *max_num_vertices = max_nvertices;
  *max_num_faces = max_nfaces;

  CHECK_EQ(world_aabb_.isValid(), true);

  if (mpi::rank() == 0) {
    std::cout << "[INFO] total face count: " << total_num_faces << std::endl;
    std::cout << "[INFO] scene bounds: " << world_aabb_ << std::endl;
  }
}

template <typename CacheT, typename SurfaceBufT>
void Scene<CacheT, SurfaceBufT>::loadAndPopulateDomainInfo(
    std::size_t* max_num_vertices, std::size_t* max_num_faces) {
  std::size_t num_domains = domains_.size();
  std::size_t total_num_models = 0;

  // count the number of configured domains and the total number of models
  std::size_t num_configured_domains = 0;
  for (std::size_t id = 0; id < num_domains; ++id) {
    total_num_models += domains_[id].getNumModels();
    const Domain& domain = domains_[id];
    num_configured_domains += domain.isConfigured();
  }

  // if not all domains are configured and there are geometric model files (i.e.
  // ply files), configure those not configured

  bool all_domains_configured = (num_configured_domains == num_domains);
  std::vector<ModelInfo> model_info_recvbuf;

  if (!all_domains_configured && total_num_models) {
    getModelInfo(total_num_models, &model_info_recvbuf);
  }

  updateDomains(all_domains_configured, model_info_recvbuf, max_num_vertices,
                max_num_faces);
}

}  // namespace spray
