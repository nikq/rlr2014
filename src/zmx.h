/*
 *
 * ZEMAX lens parser & tracer for rlr2014,
 * Copyright(c) 2014 by Hajime UCHIMURA aka nikq.
 *
 * Please contact me nikutama@gmail.com before commercial use.
 *
 * This source code is provided "as-is".
 * I am not responsible at all.
 *
 */

#ifndef __ZMX_H
#define __ZMX_H

#include <stdio.h>
#include <vector>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "vectormath.h"
#include "mtseq.h"

namespace NUMA { // レンズは沼です.

  typedef VECTORMATH::Vector3<double> Vector;

  class Intersection{
  public:
    bool   hit_;
    Vector from_;
    Vector dir_;
    Vector point_;
    Vector normal_;
    Intersection() : hit_(false) {;}
    void set( bool hit, const Vector& from, const Vector& dir, const Vector& point, const Vector& norm ){
      hit_   = hit;
      from_  = from;
      dir_   = dir;
      point_ = point;
      normal_= norm.normal();
    }
  };
  
  class Surface{
  public:

    typedef enum {
      NONE,     // 空間
      APERTURE, // 絞り面
      STANDARD, // 球面レンズ
    } TYPE;

    TYPE  type_;
    double center_;
    double radius_;
    double diameter_;
    double ior_;
    double abbe_;

    Surface(){
      init();
    }

    void init(){
      type_ = NONE;
      center_  = 0.f;
      radius_  = 0.f;
      diameter_= 0.f;
      ior_     = 1.f;
      abbe_    = 0.f;
    }

    bool intersect( const Vector& origin, const Vector& dir, Intersection& hit ){
      Vector center( center_, 0.f, 0.f );
      Vector oc = origin - center;

      double a = dir.dot( dir );
      double b = 2.f * dir.dot( oc );
      double c = oc.dot( oc ) - radius_ * radius_;
      double d2 = b*b - 4.f * a * c;

      if( d2 < 0.f )
        return false;

      double d  = sqrtf( d2 );
      double t1 = (-b - d) / (2.f * a);
      double t2 = (-b + d) / (2.f * a);
      if( t1 <= 0.f && t2 <= 0.f )
        return false;

      double t;
      if( t2 >= 0.f && t1 >= 0.f ){
        // 両方正なので両方テストする必要がある.
        Vector p1 = origin + dir * t1;
        Vector n1 = p1 - center;
        Vector p2 = origin + dir * t2;
        Vector n2 = p2 - center;
        if( radius_ < 0.f ){
          // 右(x+)向きに凸レンズ.
          // n.xが正の方を採用.
          if( n1.x > 0.f ){
            double r = sqrtf( p1.y * p1.y + p1.z * p1.z );
            if( r < diameter_ ){
              hit.set( true, origin, dir, p1, n1 );
              return true;
            }
          }else if( n2.x > 0.f ){
            double r = sqrtf( p2.y * p2.y + p2.z * p2.z );
            if( r < diameter_ ){
              hit.set( true, origin, dir, p2, n2 );
              return true;
            }
          }
          return false;
        }else{
          // 左(x-)向きに凸レンズ.
          // n.xが負の方を採用.
          if( n1.x < 0.f ){
            double r = sqrtf( p1.y * p1.y + p1.z * p1.z );
            if( r < diameter_ ){
              hit.set( true, origin, dir, p1, n1 );
              return true;
            }
          }else if( n2.x < 0.f ){
            double r = sqrtf( p2.y * p2.y + p2.z * p2.z );
            if( r < diameter_ ){
              hit.set( true, origin, dir, p2, n2 );
              return true;
            }
          }
          return false;
        }
      }else{
        t = (t1 > t2) ? t1 : t2; // どちらかが負なら大きい方が距離.
        Vector p = origin + dir * t;
        Vector n = p - center;
        if( radius_ < 0.f ){
          // 右(x+)向きに凸レンズ.
          // n.xが正の方を採用.
          if( n.x > 0.f ){
            double r = sqrtf( p.y * p.y + p.z * p.z );
            if( r < diameter_ ){
              hit.set( true, origin, dir, p, n );
              return true;
            }
          }
          return false;
        }else{
          // 左(x-)向きに凸レンズ.
          // n.xが負の方を採用.
          if( n.x < 0.f ){
            double r = sqrtf( p.y * p.y + p.z * p.z );
            if( r < diameter_ ){
              hit.set( true, origin, dir, p, n );
              return true;
            }
          }
          return false;
        }
      }
      return false;
    }
  };

  typedef std::vector<Surface> SurfaceSet;

  class Lens {
  public:

    SurfaceSet surfaces_;
    float      imageSize_; // 像高
    
    Lens() {
      surfaces_.clear();
    }

    void dump( void ){
      for(unsigned i=0;i<surfaces_.size();i++){
        printf("%d : type %d\n",i,surfaces_[i].type_);
        printf(" center %f\n",surfaces_[i].center_);
        printf(" diam   %f\n",surfaces_[i].diameter_);
        printf(" radius %f\n",surfaces_[i].radius_);
        printf(" ior  %f\n",surfaces_[i].ior_);
      }
    }
    double getLastDistance( void ){
      int i = surfaces_.size() - 1;
      return surfaces_[i].center_ - surfaces_[i].radius_;
    }

    // create a ray sample from lens cylinders.
    bool createRaySample( RLRUTIL::MTSequence& mt,double fx, double fy, double fz , Vector& origin, Vector& dir ){
      Intersection hit;

      double ior_now = 1.f; // 空気.
      Vector p( fx, fy, fz );
      Vector lookat;
      Vector d( -1.f, 0.f, 0.f );

      double r1 = 0.f, r2 = 0.f;
      double iris = 1. / 1.4;
      do{
        r1 = mt.genrand_real1();
        r2 = mt.genrand_real1();
      }while( r1 * r1 + r2 * r2 > 1.f);
      
      //printf("p %f,%f,%f\nr %f %f\n",p.x,p.y,p.z,r1,r2);

      Surface& last( surfaces_.back() );
      lookat.set( last.center_ - last.radius_,
                  last.diameter_ * (r1 - 0.5) * iris,
                  last.diameter_ * (r2 - 0.5) * iris);
      //printf("l %f,%f,%f\n",lookat.x,lookat.y,lookat.z);
      d = lookat - p;
      d.normalize();
      // printf("d %f,%f,%f\n",d.x,d.y,d.z);

      int i,i0=0,i1=0,ii=0;
      if( d.x > 0.f ){
        i1 = surfaces_.size();ii=1;
      }else{
        i0 = surfaces_.size()-1;
        i1 = 0;
        ii = -1;
      }

      for( i = i0; i != i1; i += ii ){
        // printf("trace %d\n",i);
        bool result = surfaces_[i].intersect( p, d, hit );
        if( result ){
          // next ray.
          //Vector nextDir;
          float nextIor = (ii<0) ? (i > 0 ? surfaces_[i-1].ior_ : 1.f) : surfaces_[i].ior_;
          refract( d, hit.normal_, ior_now / nextIor, d );
          d.normalize() ;
          p = hit.point_;
          ior_now = surfaces_[i].ior_;
        }else{
          break;
        }
      }
      if( i == i1 ){
        // 最後のレンズまで到達した！
        origin = p;
        dir    = d;
        return true;
      }
      return false;
    }
  };

  namespace LOADER{
    namespace ZEMAX{

      char *tokenize( char *s, char *token ) {
        char *p = (char *)s;
        while( *p && (*p == ' ' || *p == '\t') )
          p++;
        if( *p == 0 )
          return NULL;
        while(*p != ' ' && *p != '\t' && *p != '\r' && *p != '\n' && *p) {
          //printf("%c\n",*p);
          *token = *p;
          token++; p++;
        }
        *token = 0;
        if( *p == '\r' || *p == '\n' || *p == 0 )
          return NULL;
        return p;
      }

      void load( const char *filename, Lens& lens ){

        FILE *fp = fopen( filename, "rb" );
        if(!fp)return;
        int surfaceIndex = -1;
        Surface surface;

        bool isAperture = false;
        double sumz = 0.f;
        while( !feof(fp) ){
          char line[1024],*p,token[256];

          fgets(line,1024,fp);
          p = line;
          p = tokenize( p, token );
          //printf("token:%s p:%s\n",token,p);

          if( strcmp( token, "SURF" ) == 0 ){
            if( surfaceIndex >= 0 ){
              if( !isAperture && surface.diameter_ > 0.f ){
                surface.center_ = surface.center_ + surface.radius_; // 中心位置を調整しておく.
                //surface.type_ = Surface::APERTURE; // STOPは絞り.
                lens.surfaces_.push_back( surface ); // レンズフラッシュ
              }
            }
            p = tokenize( p, token );
            surfaceIndex = atoi( token );
            surface.init();
            isAperture = false;
          }

          if( strcmp( token, "TYPE" ) == 0 ){
            p = tokenize( p, token );
            if( strcmp( token, "STANDARD" ) == 0 )
              surface.type_ = Surface::STANDARD;
            else{
              printf("unknown surface: %s\n",token);
            }
          }
          if( strcmp( token, "YFLD") == 0 ){
            p = tokenize( p, token );
            p = tokenize( p, token );
            lens.imageSize_ = atof( token );
          }
          if( strcmp( token, "STOP") == 0 )
            isAperture = true;

          if( strcmp( token, "CURV" ) == 0 ){
            p = tokenize( p, token );
            double curve = atof( token );
            if( curve != 0.f )
              surface.radius_ = 1.f / curve;
            else
              surface.radius_ = 0.f;
          }

          if( strcmp( token, "DISZ" ) == 0 ){
            p = tokenize( p, token );
            double disz = atof( token );
            if( strcmp( token , "INFINITY" ) == 0 )
              disz = 0.f;
            surface.center_ = sumz / 2.f;
            sumz += disz;
          }
          if( strcmp( token, "DIAM" ) == 0 ){
            p = tokenize( p, token );
            surface.diameter_ = atof( token );
          }

          if( strcmp( token, "GLAS" ) == 0 ){
            p = tokenize( p, token ); // name
            p = tokenize( p, token ); // nazo1
            p = tokenize( p, token ); // nazo2

            p = tokenize( p, token ); // ior
            surface.ior_ = atof( token );
            p = tokenize( p, token ); // abbe
            surface.abbe_ = atof( token );
          }
        }
        fclose(fp);
      }
    }
  }
}

#endif
