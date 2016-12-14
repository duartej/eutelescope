// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#if defined USE_GEAR

// EUTelescope includes:
#include "EUTelTripletGBLUtility.h"
//#include "EUTelTripletGBLDUTscatInstance.h"

#include "EUTELESCOPE.h"

using namespace std;
using namespace eutelescope;


EUTelTripletGBLUtility::EUTelTripletGBLUtility(){}

TMatrixD EUTelTripletGBLUtility::JacobianPointToPoint( double ds ) {
  /* for GBL:
     Jacobian for straight line track
     track = q/p, x', y', x, y
     0,   1,  2,  3, 4
     */
  TMatrixD jac( 5, 5 );
  jac.UnitMatrix();
  jac[3][1] = ds; // x = x0 + xp * ds
  jac[4][2] = ds; // y = y0 + yp * ds
  return jac;
}

EUTelTripletGBLUtility::track::track(triplet up, triplet down) : upstream(up), downstream(down) {}

double EUTelTripletGBLUtility::track::kink_x() {
  return (downstream.slope().x - upstream.slope().x);
}

double EUTelTripletGBLUtility::track::kink_y() {
  return (downstream.slope().y - upstream.slope().y);
}

EUTelTripletGBLUtility::hit EUTelTripletGBLUtility::track::intersect() {
  hit inter;
  // Re-check what this actually is...
  // and simplify using triplet class members...
  inter.x = ( upstream.base().x - upstream.slope().x * upstream.base().z - downstream.base().x + downstream.slope().x * downstream.base().z ) / kink_x();
  inter.y = ( upstream.base().y - upstream.slope().y * upstream.base().z - downstream.base().y + downstream.slope().y * downstream.base().z ) / kink_y();
  return inter;
}

EUTelTripletGBLUtility::triplet EUTelTripletGBLUtility::track::get_upstream() {
  return upstream;
}

EUTelTripletGBLUtility::triplet EUTelTripletGBLUtility::track::get_downstream() {
  return downstream;
}

EUTelTripletGBLUtility::hit EUTelTripletGBLUtility::track::gethit(int plane) {
  if(plane < 3) return upstream.gethit(plane);
  else return downstream.gethit(plane);
}



EUTelTripletGBLUtility::triplet::triplet() : linked_dut(false), hits() {
  // Empty default constructor
}

EUTelTripletGBLUtility::triplet::triplet(hit hit0, hit hit1, hit hit2) : linked_dut(false), hits() {
  triplet();
  filltriplet(hit0, hit1, hit2);
}

EUTelTripletGBLUtility::hit EUTelTripletGBLUtility::triplet::getpoint_at(double z) {
  hit impact;
  impact.z = z - base().z;
  impact.x = base().x + slope().x * impact.z;
  impact.y = base().y + slope().y * impact.z;
  return impact;
}

double EUTelTripletGBLUtility::triplet::getx_at(double z) {
  return this->base().x + this->slope().x * (z - this->base().z);
}

double EUTelTripletGBLUtility::triplet::getdx() {
  return this->hits.rbegin()->second.x - this->hits.begin()->second.x;
}

double EUTelTripletGBLUtility::triplet::getdx(int ipl) {
  return this->hits[ipl].x - this->base().x - this->slope().x * (this->hits[ipl].z - this->base().z);
}

double EUTelTripletGBLUtility::triplet::getdx(hit point) {
  return point.x - this->base().x - this->slope().x * (point.z - this->base().z);
}

double EUTelTripletGBLUtility::triplet::gety_at(double z) {
  return this->base().y + this->slope().y * (z - this->base().z);
}

double EUTelTripletGBLUtility::triplet::getdy() {
  return this->hits.rbegin()->second.y - this->hits.begin()->second.y;
}

double EUTelTripletGBLUtility::triplet::getdy(int ipl) {
  return this->hits[ipl].y - this->base().y - this->slope().y * (this->hits[ipl].z - this->base().z);
}

double EUTelTripletGBLUtility::triplet::getdy(hit point) {
  return point.y - this->base().y - this->slope().y * (point.z - this->base().z);
}

double EUTelTripletGBLUtility::triplet::getdz() {
  return this->hits.rbegin()->second.z - this->hits.begin()->second.z;
}

EUTelTripletGBLUtility::hit EUTelTripletGBLUtility::triplet::gethit(int plane) {
  return this->hits[plane];
}

EUTelTripletGBLUtility::hit EUTelTripletGBLUtility::triplet::base() {
  hit center;
  center.x = 0.5*( this->hits.begin()->second.x + this->hits.rbegin()->second.x );
  center.y = 0.5*( this->hits.begin()->second.y + this->hits.rbegin()->second.y );
  center.z = 0.5*( this->hits.begin()->second.z + this->hits.rbegin()->second.z );
  return center;
}

EUTelTripletGBLUtility::hit EUTelTripletGBLUtility::triplet::slope() {
  hit sl;
  double dz = (this->hits.rbegin()->second.z - this->hits.begin()->second.z);
  sl.x = (this->hits.rbegin()->second.x - this->hits.begin()->second.x) / dz;
  sl.y = (this->hits.rbegin()->second.y - this->hits.begin()->second.y) / dz;
  return sl;
}

#endif
