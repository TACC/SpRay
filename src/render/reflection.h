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

#pragma once

#include "embree/random_sampler.h"

#include "glm/glm.hpp"
#include "glog/logging.h"

#include "render/sampler.h"
#include "utils/math.h"

namespace spray {

enum SampleType {
  BSDF_REFLECTION = 1 << 0,
  BSDF_TRANSMISSION = 1 << 1,
  BSDF_ALL = BSDF_REFLECTION | BSDF_TRANSMISSION
};

inline bool hasReflection(int sample_type) {
  return ((sample_type & BSDF_REFLECTION) == BSDF_REFLECTION);
}

inline bool hasTransmission(int sample_type) {
  return ((sample_type & BSDF_TRANSMISSION) == BSDF_TRANSMISSION);
}

//! Evaluates reflection ray.
/*!
  \param I Incident ray direction.
  \param N Surface normal.
  \return Returns outgoing ray direction (reflection).
*/
glm::vec3 reflect(const glm::vec3& I, const glm::vec3& N);

//! Evaluates Snell's law.
/*!
  \param I Incident ray direction.
  \param N Surface normal.
  \param n1 Refractive index of the input medium.
  \param n2 Refractive index of the output medium.
  \param T Outgoing ray direction (refraction).
  \return Returns true for refraction and false for total internal reflection.
*/
bool refract(const glm::vec3& I, const glm::vec3& N, float n1, float n2,
             glm::vec3* T);

//! Evaluates Fresnel's law.
/*!
  \param I Incident ray direction.
  \param N Surface normal.
  \param n1 Refractive index of the input medium.
  \param n2 Refractive index of the output medium.
  \return Returns probability of reflection.
*/
float reflectanceFresnel(const glm::vec3& I, const glm::vec3& N, float n1,
                         float n2);

//! Evaluates Schlick's law.
/*!
  \param I Incident ray direction.
  \param N Surface normal.
  \param n1 Refractive index of the input medium.
  \param n2 Refractive index of the output medium.
  \return Returns probability of reflection.
*/
float reflectanceSchlick(const glm::vec3& I, const glm::vec3& N, float n1,
                         float n2);

// BSDF Inline Functions
inline float CosTheta(const glm::vec3& w) { return w.z; }
inline float Cos2Theta(const glm::vec3& w) { return w.z * w.z; }
inline float AbsCosTheta(const glm::vec3& w) { return glm::abs(w.z); }
inline float Sin2Theta(const glm::vec3& w) {
  return glm::max(0.0f, 1.0f - Cos2Theta(w));
}

inline float SinTheta(const glm::vec3& w) { return glm::sqrt(Sin2Theta(w)); }

inline float TanTheta(const glm::vec3& w) { return SinTheta(w) / CosTheta(w); }

inline float Tan2Theta(const glm::vec3& w) {
  return Sin2Theta(w) / Cos2Theta(w);
}

inline float CosPhi(const glm::vec3& w) {
  float sinTheta = SinTheta(w);
  return (sinTheta == 0.0f) ? 1.0f : glm::clamp(w.x / sinTheta, -1.0f, 1.0f);
}

inline float SinPhi(const glm::vec3& w) {
  float sinTheta = SinTheta(w);
  return (sinTheta == 0.0f) ? 0.0f : glm::clamp(w.y / sinTheta, -1.0f, 1.0f);
}

inline float Cos2Phi(const glm::vec3& w) { return CosPhi(w) * CosPhi(w); }

inline float Sin2Phi(const glm::vec3& w) { return SinPhi(w) * SinPhi(w); }

inline float CosDPhi(const glm::vec3& wa, const glm::vec3& wb) {
  return glm::clamp(
      (wa.x * wb.x + wa.y * wb.y) /
          glm::sqrt((wa.x * wa.x + wa.y * wa.y) * (wb.x * wb.x + wb.y * wb.y)),
      -1.0f, 1.0f);
}

inline glm::vec3 Reflect(const glm::vec3& wo, const glm::vec3& normal_ff) {
  return (-wo + ((2.0f * glm::dot(wo, normal_ff)) * normal_ff));
}

inline bool SameHemisphere(const glm::vec3& w, const glm::vec3& wp) {
  return w.z * wp.z > 0;
}

inline float FrDielectric(float cosThetaI, float etaI, float etaT,
                          const glm::vec3& wo, const glm::vec3& normal_ff,
                          glm::vec3* wt) {
  float sin2ThetaI = glm::max(0.0f, 1.0f - (cosThetaI * cosThetaI));
  float eta = etaI / etaT;
  float sin2ThetaT = eta * eta * sin2ThetaI;

  // Handle total internal reflection
  if (sin2ThetaT >= 1.0f) return 1.0f;

  // float cosThetaT = glm::sqrt(glm::max(0.0f, 1.0f - sinThetaT * sinThetaT));
  float cosThetaT = glm::sqrt(1.0f - sin2ThetaT);
  float Rparl = ((etaT * cosThetaI) - (etaI * cosThetaT)) /
                ((etaT * cosThetaI) + (etaI * cosThetaT));
  float Rperp = ((etaI * cosThetaI) - (etaT * cosThetaT)) /
                ((etaI * cosThetaI) + (etaT * cosThetaT));

  // probability of reflection
  float fr = (Rparl * Rparl + Rperp * Rperp) / 2.0f;

  // direction of transmitted ray
  *wt = (eta * -wo) + (normal_ff * (eta * cosThetaI - cosThetaT));

  return fr;
}

inline bool Refract(float cosThetaI, float etaI, float etaT,
                    const glm::vec3& wo, const glm::vec3& normal_ff,
                    glm::vec3* wt) {
  // Compute _cosThetaT_ using Snell's law
  float sin2ThetaI = glm::max(0.0f, 1.0f - (cosThetaI * cosThetaI));
  float eta = etaI / etaT;
  float sin2ThetaT = eta * eta * sin2ThetaI;

  // Handle total internal reflection
  if (sin2ThetaT >= 1.0f) return false;
  float cosThetaT = glm::sqrt(1.0f - sin2ThetaT);
  // direction of transmitted ray
  *wt = (eta * -wo) + (normal_ff * (eta * cosThetaI - cosThetaT));
  return true;
}

// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
inline glm::vec3 FrConductor(float cosThetaI, const glm::vec3& etai,
                             const glm::vec3& etat, const glm::vec3& k) {
  cosThetaI = glm::clamp(cosThetaI, -1.0f, 1.0f);
  glm::vec3 eta = etat / etai;
  glm::vec3 etak = k / etai;

  float cosThetaI2 = cosThetaI * cosThetaI;
  float sinThetaI2 = 1. - cosThetaI2;
  glm::vec3 eta2 = eta * eta;
  glm::vec3 etak2 = etak * etak;

  glm::vec3 t0 = eta2 - etak2 - sinThetaI2;
  glm::vec3 a2plusb2 = glm::sqrt(t0 * t0 + 4.0f * eta2 * etak2);
  glm::vec3 t1 = a2plusb2 + cosThetaI2;
  glm::vec3 a = glm::sqrt(0.5f * (a2plusb2 + t0));
  glm::vec3 t2 = (float)2 * cosThetaI * a;
  glm::vec3 Rs = (t1 - t2) / (t1 + t2);

  glm::vec3 t3 = cosThetaI2 * a2plusb2 + sinThetaI2 * sinThetaI2;
  glm::vec3 t4 = t2 * sinThetaI2;
  glm::vec3 Rp = Rs * (t3 - t4) / (t3 + t4);

  return 0.5f * (Rp + Rs);
}

inline glm::vec3 blinnPhong(const float costheta, const glm::vec3 kd,
                            const glm::vec3 ks, const float shininess,
                            const glm::vec3& li, const glm::vec3& wi,
                            const glm::vec3& n_hat, const glm::vec3 wo) {
  //
  glm::vec3 half_hat = glm::normalize(wi + wo);
  float n_dot_h = glm::clamp(glm::dot(n_hat, half_hat), 0.0f, 1.0f);

  glm::vec3 cs = ks * glm::pow(n_dot_h, shininess);
  glm::vec3 cd = kd * costheta;

  return (li * (cd + cs));
}

/**
 * Interface for BSDFs.
 */
class Bsdf {
 public:
  virtual ~Bsdf() {}
  virtual void sampleDelta(bool entering, float cos_theta_i,
                           const glm::vec3& wo, const glm::vec3& normal_ff,
                           uint32_t* sample_type, float* fr,
                           glm::vec3* wt) const = 0;

  virtual void sampleRandom(const glm::vec3& normal, RandomSampler* sampler,
                            glm::vec3* wi, float* pdf) const = 0;
  /**
   * Get the emission value of the surface material. For non-emitting surfaces
   * this would be a zero energy spectrum.
   * \return emission spectrum of the surface material
   */
  virtual glm::vec3 getEmission() const = 0;

  // virtual glm::vec3 getCoefficient() const = 0;

  virtual bool isDelta() const = 0;
};

class DiffuseBsdf : public Bsdf {
 public:
  DiffuseBsdf(const glm::vec3& albedo) : albedo_(albedo) {}

  void sampleRandom(const glm::vec3& normal, RandomSampler* sampler,
                    glm::vec3* wi, float* pdf) const override {
    glm::vec2 u = RandomSampler_get2D(*sampler);
    getCosineHemisphereSample(u.x, u.y, normal, wi, pdf);
  }

  void sampleDelta(bool entering, float cos_theta_i, const glm::vec3& wo,
                   const glm::vec3& normal_ff, uint32_t* sample_type, float* fr,
                   glm::vec3* wt) const override {
    LOG(FATAL) << "forbidden";
  }

  glm::vec3 getEmission() const override { return glm::vec3(0.0f); }
  bool isDelta() const override { return false; }

 private:
  glm::vec3 albedo_;
};

class MirrorBsdf : public Bsdf {
 public:
  MirrorBsdf(const glm::vec3& reflectance) : reflectance_(reflectance) {}

  glm::vec3 getEmission() const override { return glm::vec3(0.0f); }
  bool isDelta() const override { return true; }
  void sampleRandom(const glm::vec3& normal, RandomSampler* sampler,
                    glm::vec3* wi, float* pdf) const override {
    LOG(FATAL) << "forbidden";
  }
  void sampleDelta(bool entering, float cos_theta_i, const glm::vec3& wo,
                   const glm::vec3& normal_ff, uint32_t* sample_type, float* fr,
                   glm::vec3* wt) const override {
    *fr = 1.0f;
    *sample_type = BSDF_REFLECTION;
  }

 private:
  glm::vec3 reflectance_;
};

/**
 * Glass BSDF.
 */
class GlassBsdf : public Bsdf {
 public:
  GlassBsdf(float eta_exterior, float eta_interior)
      : etaI_(eta_exterior), etaT_(eta_interior) {}
  glm::vec3 getEmission() const override { return glm::vec3(0.0f); }
  bool isDelta() const override { return true; }

  void sampleRandom(const glm::vec3& normal, RandomSampler* sampler,
                    glm::vec3* wi, float* pdf) const override {
    LOG(FATAL) << "forbidden";
  }

  void sampleDelta(bool entering, float cos_theta_i, const glm::vec3& wo,
                   const glm::vec3& normal_ff, uint32_t* sample_type, float* fr,
                   glm::vec3* wt) const override {
    float etaI = entering ? etaI_ : etaT_;
    float etaT = entering ? etaT_ : etaI_;

    float prob_reflect =
        FrDielectric(cos_theta_i, etaI, etaT, wo, normal_ff, wt);
    *fr = prob_reflect;

    if (prob_reflect == 1.0f) {  // reflection only (total internal reflection)
      *sample_type = BSDF_REFLECTION;

    } else if (prob_reflect == 0.0f) {  // refraction only
      *sample_type = BSDF_TRANSMISSION;

    } else {  // both reflection and refraction
#ifdef SPRAY_GLOG_CHECK
      CHECK_GT(prob_reflect, 0.0f);
#endif
      *sample_type = BSDF_REFLECTION | BSDF_TRANSMISSION;
    }
  }

 private:
  float etaI_, etaT_;
};

/**
 * Refraction BSDF.
 */
class TransmissionBsdf : public Bsdf {
 public:
  TransmissionBsdf(float eta_exterior, float eta_interior)
      : etaI_(eta_exterior), etaT_(eta_interior) {}
  glm::vec3 getEmission() const override { return glm::vec3(0.0f); }
  bool isDelta() const override { return true; }

  void sampleRandom(const glm::vec3& normal, RandomSampler* sampler,
                    glm::vec3* wi, float* pdf) const override {
    LOG(FATAL) << "forbidden";
  }

  void sampleDelta(bool entering, float cos_theta_i, const glm::vec3& wo,
                   const glm::vec3& normal_ff, uint32_t* sample_type, float* fr,
                   glm::vec3* wt) const override {
    float etaI = entering ? etaI_ : etaT_;
    float etaT = entering ? etaT_ : etaI_;
    // float eta = etaI / etaT;

    bool transmission = Refract(cos_theta_i, etaI, etaT, wo, normal_ff, wt);

    if (transmission) {  // refraction only
      *sample_type = BSDF_TRANSMISSION;
      *fr = 0.0f;  // prob. of reflection

    } else {  // reflection only (total internal reflection)
      *sample_type = BSDF_REFLECTION;
      *fr = 1.0f;  // prob. of reflection
    }
  }

 private:
  float etaI_, etaT_;
};

}  // namespace spray

