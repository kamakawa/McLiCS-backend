#pragma once

#ifdef __CUDACC__
  #define HD __host__ __device__
  #define DINL __device__ __forceinline__
#else
  #define HD
  #define DINL inline
#endif

#include <cstdint>
#include <cmath>   // cosf, sinf, sqrtf, fabsf

// ============================
// Anchoring types (GPU-safe)
// ============================
enum class AnchoringType : int32_t {
  Strong = 0,
  Strong_GHRL = 1,
  RP = 2,
  RP_GHRL = 3,
  FG = 4,
  FG_GHRL = 5,
  Homeotropic = 6,
  Homeotropic_GHRL = 7
};

// ============================
// POD params uploaded to GPU
// ============================
// NOVO: layout previsível (ordem e alinhamento)
struct alignas(16) AnchoringDeviceParams {
  AnchoringType type;
  int32_t _pad0;   // NOVO: padding explícito p/ manter alinhamento e layout estável
  float W;
  float theta_s;   // graus
  float phi_s;     // graus
  // (se no futuro precisar de mais coisas, adiciona aqui)
};

// ============================
// Small math helpers
// ============================
DINL HD float deg2rad(float deg) {
  // NOVO: constante local (evita depender de M_PI em device)
  constexpr float PI = 3.14159265358979323846f;
  return deg * (PI / 180.0f);
}

DINL HD float dot3(const float a[3], const float b[3]) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

DINL HD float signf_nz(float x) { // NOVO: sinal como no CPU (nij > 0 ? 1 : -1)
  return (x > 0.0f) ? 1.0f : -1.0f;
}

// Normaliza vetor 3D (robusto)
DINL HD void normalize3(float v[3]) {
  float n2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
  if (n2 > 0.0f) {
    float inv = 1.0f / sqrtf(n2);
    v[0] *= inv; v[1] *= inv; v[2] *= inv;
  }
}

// ============================
// Surface anchoring potentials
// ============================
// IMPORTANT: no CPU você tem várias estratégias.
// Na GPU vamos roteá-las por enum.
//
// Assinatura mínima: (ni, s) -> energia.
// - ni: diretor local (unit)
// - s : normal local da superfície (unit) (ou direção “surface axis” usada no seu modelo)
// - p : params do ancoramento (W, theta_s, phi_s, tipo)
// - ghrl_*: constantes usadas nas variantes GHRL
//
// OBS: No CPU, alguns casos "Strong" usam (theta_s, phi_s) para definir n_s fixo.
// Já "Homeotropic" usa s (normal local).
// RP/FG aqui mantemos como no seu código atual (forma quadrática em dot).
//
// NOVO: você pode chamar isso tanto no host quanto no device (HD).
DINL HD float anchoring_surface_potential(
    const AnchoringDeviceParams& p,
    const float ni[3],
    const float s[3],
    float ghrl_lambda,
    float ghrl_mu,
    float ghrl_nu,
    float ghrl_rho,
    float ghrl_sigma
) {
  // Early exit
  if (p.W == 0.0f) return 0.0f;

  // Constrói n_s se precisar (Strong / RP usam theta/phi no seu CPU)
  // No seu CPU, Strong/RP usam:
  // n_s = (cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta))
  float n_s[3];
  {
    float th = deg2rad(p.theta_s);
    float ph = deg2rad(p.phi_s);
    float sTH = sinf(th);
    n_s[0] = cosf(ph) * sTH;
    n_s[1] = sinf(ph) * sTH;
    n_s[2] = cosf(th);
  }

  // dot products
  const float ni_dot_ns = dot3(ni, n_s);
  const float ni_dot_s  = dot3(ni, s);

  switch (p.type) {

    // =========================
    // Strong (CPU: -W*(ni·n_s)^2)
    // =========================
    case AnchoringType::Strong:
    case AnchoringType::Strong_GHRL: {
      return -p.W * ni_dot_ns * ni_dot_ns;
    }

    // =========================
    // Rapini-Papoular (no seu CPU atual está igual ao Strong)
    // -W*(ni·n_s)^2
    // =========================
    case AnchoringType::RP: {
      return -p.W * ni_dot_ns * ni_dot_ns;
    }

    // =========================
    // RP_GHRL (copiado do seu CPU)
    // nj = n_s, s = normal (argumento)
    // =========================
    case AnchoringType::RP_GHRL: {
      const float el = ghrl_lambda;
      const float em = ghrl_mu;
      const float en = ghrl_nu;
      const float er = ghrl_rho;
      const float es = ghrl_sigma;

      const float v15 = 1.5f;
      const float v05 = 0.5f;

      // ai = ni·s ; aj = nj·s ; nij = ni·nj
      const float ai  = ni_dot_s;
      const float aj  = dot3(n_s, s);
      const float nij = dot3(ni, n_s);

      const float pij = v15 * nij * nij - v05;

      // cross = (ni x nj) · s
      const float cross =
          (ni[2]*n_s[1] - ni[1]*n_s[2]) * s[0] +
          (ni[0]*n_s[2] - ni[2]*n_s[0]) * s[1] +
          (ni[1]*n_s[0] - ni[0]*n_s[1]) * s[2];

      const float E1 = (v15*ai*ai) + (v15*aj*aj) - 1.0f;

      // CPU: return W * ((E1*(er*pij+el) + em*(ai*aj*nij) - 1/9) + en*pij + es*sign(nij)*cross)
      return p.W * ((E1 * (er * pij + el) + em * (ai * aj * nij) - (1.0f / 9.0f)) +
                    en * pij + es * signf_nz(nij) * cross);
    }

    // =========================
    // Fournier-Galatola (CPU atual: +W*(ni·s)^2)
    // =========================
    case AnchoringType::FG: {
      return +p.W * ni_dot_s * ni_dot_s;
    }

    // =========================
    // FG_GHRL (copiado do seu CPU)
    // =========================
    case AnchoringType::FG_GHRL: {
      const float el = ghrl_lambda;
      const float em = ghrl_mu;
      const float en = ghrl_nu;
      const float er = ghrl_rho;
      const float es = ghrl_sigma;

      const float v15 = 1.5f;
      const float v05 = 0.5f;

      // CPU:
      // mod = sqrt(|ni·s|)
      // nj = ni - s*mod
      // normalize(nj)
      float mod = sqrtf(fabsf(ni_dot_s));
      float nj[3] = { ni[0] - s[0]*mod, ni[1] - s[1]*mod, ni[2] - s[2]*mod };
      normalize3(nj);

      const float ai  = ni_dot_s;
      const float aj  = dot3(nj, s);
      const float nij = dot3(ni, nj);

      const float pij = v15 * nij * nij - v05;

      // cross = (ni x nj) · s
      const float cross =
          (ni[2]*nj[1] - ni[1]*nj[2]) * s[0] +
          (ni[0]*nj[2] - ni[2]*nj[0]) * s[1] +
          (ni[1]*nj[0] - ni[0]*nj[1]) * s[2];

      const float E1 = (v15*ai*ai) + (v15*aj*aj) - 1.0f;

      return p.W * ((E1 * (er * pij + el) + em * (ai * aj * nij) - (1.0f / 9.0f)) +
                    en * pij + es * signf_nz(nij) * cross);
    }

    // =========================
    // Homeotropic (CPU: -W*(ni·s)^2)
    // =========================
    case AnchoringType::Homeotropic: {
      return -p.W * ni_dot_s * ni_dot_s;
    }

    // =========================
    // Homeotropic_GHRL (copiado do seu CPU)
    // =========================
    case AnchoringType::Homeotropic_GHRL: {
      const float el = ghrl_lambda;
      const float em = ghrl_mu;
      const float en = ghrl_nu;
      const float er = ghrl_rho;
      const float es = ghrl_sigma;

      const float v15 = 1.5f;
      const float v05 = 0.5f;

      // CPU:
      // ai = ni·s
      // aj = 1
      // nij = ni·s
      const float ai  = ni_dot_s;
      const float aj  = 1.0f;
      const float nij = ni_dot_s;

      const float pij = v15 * nij * nij - v05;

      // cross = (ni x s) · s  (no CPU está assim, mantém compatibilidade)
      const float cross =
          (ni[2]*s[1] - ni[1]*s[2]) * s[0] +
          (ni[0]*s[2] - ni[2]*s[0]) * s[1] +
          (ni[1]*s[0] - ni[0]*s[1]) * s[2];

      const float E1 = (v15*ai*ai) + (v15*aj*aj) - 1.0f;

      return p.W * ((E1 * (er * pij + el) + em * (ai * aj * nij) - (1.0f / 9.0f)) +
                    en * pij + es * signf_nz(nij) * cross);
    }

    default:
      return 0.0f;
  }
}