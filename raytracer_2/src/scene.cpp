#include "scene.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "hit.h"
#include "image.h"
#include "material.h"
#include "ray.h"

using namespace std;

pair<ObjectPtr, Hit> Scene::castRay(Ray const &ray) const {
  // Find hit object and distance
  Hit min_hit(numeric_limits<double>::infinity(), Vector());
  ObjectPtr obj = nullptr;
  for (unsigned idx = 0; idx != objects.size(); ++idx) {
    Hit hit(objects[idx]->intersect(ray));
    if (hit.t < min_hit.t) {
      min_hit = hit;
      obj = objects[idx];
    }
  }

  return pair<ObjectPtr, Hit>(obj, min_hit);
}

Color Scene::trace(Ray const &ray, unsigned depth) {
  pair<ObjectPtr, Hit> mainhit = castRay(ray);
  ObjectPtr obj = mainhit.first;
  Hit min_hit = mainhit.second;

  // No hit? Return background color.
  if (!obj) return Color(0.0, 0.0, 0.0);

  Material const &material = obj->material;
  Point hit = ray.at(min_hit.t);
  Vector V = -ray.D;

  // Pre-condition: For closed objects, N points outwards.
  Vector N = min_hit.N;

  // The shading normal always points in the direction of the view,
  // as required by the Phong illumination model.
  Vector shadingN;
  if (N.dot(V) >= 0.0)
    shadingN = N;
  else
    shadingN = -N;

  Color matColor = material.color;

  // Add ambient once, regardless of the number of lights.
  Color color = material.ka * matColor;

  // Add diffuse and specular components.
  for (auto const &light : lights) {
    Vector L = (light->position - hit).normalized();

    // check for shadows
    if(renderShadows) {
      // offset shadow ray using shading normal
      //Ray shadowRay(hit + epsilon * shadingN, L);
      Ray shadowRay(hit + epsilon * shadingN, L);
      pair<ObjectPtr, Hit> shadowHit = castRay(shadowRay);

      // if the shadow ray hits an object
      // and does not intersect with the object of hitpoint (self-shadowing)
      if (shadowHit.first && shadowHit.second.t < (light->position - hit).length()) {
        // then the point is in shadow
        // skip diffuse and specular light
        continue;
      }
    }
    // Add diffuse.
    double dotNormal = shadingN.dot(L);
    double diffuse = std::max(dotNormal, 0.0);
    color += diffuse * material.kd * light->color * matColor;

    // Add specular.
    if (dotNormal > 0) {
      Vector reflectDir = reflect(-L, shadingN);
      double specAngle = std::max(reflectDir.dot(V), 0.0);
      double specular = std::pow(specAngle, material.n);

      color += specular * material.ks * light->color;
    }
  }

  // Reflection and refraction (doesn't work yet)
  if (depth > 0 and material.isTransparent) {
    // The object is transparent, and thus refracts and reflects light.
    
    // // Compute direction of the reflection ray
    Vector specReflectDir = reflect(-V, shadingN);
    Ray reflectRay(hit + epsilon * shadingN, specReflectDir);

    // Compute direction of the refracted ray
    Vector refractDir = refract(-V, shadingN, material.nt);
    Ray refractRay(hit - epsilon * shadingN, refractDir);

    // // Compute the specular reflection coefficient
    // // using Schlick's approximation
    double ratio = schlick(1.0, material.nt, shadingN, -V);
    // // Recursively trace the reflected and refracted rays
    color += trace(reflectRay, depth - 1);
    color += trace(refractRay, depth - 1);

    // Reflection only
  } else if (depth > 0 and material.ks > 0.0) {
    // The object is not transparent, but opaque.
    // Compute direction of the reflection ray
    Vector specReflectDir = reflect(-V, shadingN);
    // Avoid intersection with current object using eps
    Ray reflectRay(hit + epsilon * shadingN, specReflectDir);
    // Recursively trace the reflected ray
    color += material.ks * trace(reflectRay, depth - 1);
  }

  return color;
}

void Scene::render(Image &img) {
  unsigned w = img.width();
  unsigned h = img.height();

  for (unsigned y = 0; y < h; ++y)
    for (unsigned x = 0; x < w; ++x) {
      Point pixel(x + 0.5, h - 1 - y + 0.5, 0);
      Ray ray(eye, (pixel - eye).normalized());
      Color col = trace(ray, recursionDepth);
      col.clamp();
      img(x, y) = col;
    }
}

// --- Misc functions ----------------------------------------------------------

// Defaults
Scene::Scene() : objects(), lights(), eye(), renderShadows(false), recursionDepth(0), supersamplingFactor(1) {}

void Scene::addObject(ObjectPtr obj) { objects.push_back(obj); }

void Scene::addLight(Light const &light) { lights.push_back(LightPtr(new Light(light))); }

void Scene::setEye(Triple const &position) { eye = position; }

unsigned Scene::getNumObject() { return objects.size(); }

unsigned Scene::getNumLights() { return lights.size(); }

void Scene::setRenderShadows(bool shadows) { renderShadows = shadows; }

void Scene::setRecursionDepth(unsigned depth) { recursionDepth = depth; }

void Scene::setSuperSample(unsigned factor) { supersamplingFactor = factor; }
