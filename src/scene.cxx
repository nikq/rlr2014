/*
 * scene.cxx
 *
 * implements SPPM integrator, main loop.
 *
 */

#include "vectormath.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "scene.h"
#include "modelholder.h"
#include "bvh.h"
#include "kdtree.h"
#include "mtseq.h"

#include "zmx.h"

#include <time.h>

__NS_MTSEQ::MTSequence global_rand;

class Camera{
private:
public:
  Camera(){; }
  virtual ~Camera(){;}
  
  VECTORMATH::Vector pos_, lookat_, dir_, dirU_, dirV_;
  NUMA::Lens lens_;
  
  double film_z_;     // フィルム距離
  double film_w_;
  double film_h_;

  void setup(
    VECTORMATH::Vector p,
    VECTORMATH::Vector l,
    NUMA::Lens lens,
    double z = 1.,
    double w = 0.,
    double h = 1. ){
    VECTORMATH::Vector v;
    lens_   = lens;
    pos_    = p;
    lookat_ = l;
    dir_    = lookat_ - pos_;
    
    printf("[CAM] position : "); pos_.dump();
    printf("[CAM] lookat   : "); lookat_.dump();
    //printf("[CAM] direction: "); dir_.dump();
    //dist = dir.dot( dir );

    dir_.normalize();
    printf("[CAM] direction: "); dir_.dump();

    v.set( 0.0, 1.0, 0.0 );
    dirU_ = dir_.cross( v );
    dirU_.normalize();
    dirV_ = dirU_.cross( dir_ );
    dirV_.normalize();
    film_z_ = z;
    film_w_ = w;
    film_h_ = h;
    
    printf("[CAM]  world x : "); dirU_.dump();
    printf("[CAM]  world y : "); dirV_.dump();
    printf("[CAM]  world z : "); dir_.dump();
  }
  bool ray( __NS_MTSEQ::MTSequence& mt, double x, double y, RLR::Ray& r ){
    VECTORMATH::Vector origin, dir;
    bool result = lens_.createRaySample( mt, film_z_, film_w_ * (x-0.5) * -1, film_h_ * (y-0.5) * -1., origin, dir );
    if( !result )
      return false;
    r.pos_ = pos_ + dir_ * -origin.x + dirU_ * origin.y + dirV_ * origin.z;
    r.dir_ = dir_ * -dir.x + dirU_ * dir.y + dirV_ * dir.z;
    //printf("%f %f > %f %f %f / %f %f %f\n",x,y,origin.x,origin.y,origin.z,dir.x,dir.y,dir.z);
    //printf(" %f %f %f\n",r.pos_.x,r.pos_.y,r.pos_.z);
    //printf(" %f %f %f\n",r.dir_.x,r.dir_.y,r.dir_.z);
    return true;
  }
};


template<typename ACCEL>
class SPPM{
public:

  typedef struct HitPoint_t {
    __NS_RLR::Color  adj;
    __NS_RLR::Color  flux;
    __NS_RLR::Vector pos;
    __NS_RLR::Vector norm;
    double       r2;
    unsigned int n;
    int          pix;
    int          index;
  public: __NS_RLR::Vector& position(){ return pos; }
  } HitPoint;
  
  typedef SPPM<ACCEL> self;
  typedef __NS_KDTREE::KdTree<HitPoint> HitPointKD;
  
  double max_r2_;
  double initial_radius_;
  double total_photons_;
  
  HitPointKD            hitpoint_kd_;

  __NS_RLR::Scene              &scene_;
  ACCEL                        &accel_;
  Camera             &camera_;
  __NS_RLR::FrameBuffer< self > film_;
  
  SPPM(__NS_RLR::Scene &scene, ACCEL& accel, Camera &camera, int width, int height, int samples )
       :  total_photons_(0.),
          scene_(scene),
          accel_(accel),
          camera_(camera){
            film_.setup( width, height, samples);
            __NS_RLR::AABB abb = scene_.getAABB();
            __NS_RLR::Vector v = abb.hi() - abb.lo();
            //initial_radius_   = ((v.x + v.y + v.z) / 3.0) / ((film_.width_ + film_.height_) / 2.0) * 10.;
            initial_radius_   = ((v.x + v.y + v.z) / 3.0) / ((film_.width_ + film_.height_) / 2.0) * 2.0;
            max_r2_ = initial_radius_ * initial_radius_;
          }

  inline double MAX(double a,double b){ return (a<b)?b:a; }

  bool trace( const __NS_RLR::Ray &ray, int depth, int mode, __NS_RLR::Color adj, int pix ) {

    depth++;
    __NS_RLR::Intersection isect;

    accel_.testRAY( ray, isect );

    if( !isect.hit() || depth >= 16 ){
      return false;
    }

    __NS_RLR::Vector hit  ( isect.getPoint()  );
    __NS_RLR::Vector norm ( isect.getNormal() );

    __NS_RLR::Vector nl   = (norm.dot(ray.direction())) < 0. ? norm : norm * -1.; // ray.dirと逆方向のnorm
    __NS_RLR::SceneObject   * object   = isect.getObject();
    __NS_RLR::SceneMaterial * material = object->getMaterial();
    __NS_RLR::Color           color( material->color_ );
    __NS_RLR::Color           emit ( material->color_ * material->emitter_ );


    float r = global_rand.genrand_real1();
    if( r < material->diffuse_ ){
      
      if( mode == 1 ){

        int index = hitpoint_kd_.find( pix );
        if( index >= 0 ){
          hitpoint_kd_.get(index).adj = color.mul(adj);
          hitpoint_kd_.get(index).pos = hit;
          hitpoint_kd_.get(index).norm= norm;
        }else{
          HitPoint hitpoint;
          hitpoint.n   = 0;
          hitpoint.pix = pix;
          hitpoint.adj = color.mul( adj );
          hitpoint.pos = hit;
          hitpoint.norm= norm;
          hitpoint.flux= __NS_RLR::Vector(0.,0.,0.);
          hitpoint.r2  = initial_radius_ * initial_radius_; // 適当
          hitpoint_kd_.add( hitpoint );
        }
        
        return true;
        
      }else{
        
        typename HitPointKD::LIST result;
        hitpoint_kd_.query( hit, max_r2_, result );
        
        bool stored = false;
        if( result.size() > 0 ){
          // PPM integration.
          for( typename HitPointKD::LIST::iterator it = result.begin(); it != result.end(); ++it ){
            HitPoint& hitpoint( hitpoint_kd_.get( *it ) );
            __NS_RLR::Vector dist = hitpoint.pos - hit;
            if ((hitpoint.norm.dot(norm) > 1e-3) && (dist.dot(dist) <= hitpoint.r2)) {
              const double ALPHA = 0.7;
              double g = (hitpoint.n * ALPHA + ALPHA) / (hitpoint.n * ALPHA + 1.0);
              hitpoint.r2 = hitpoint.r2 * g;
              hitpoint.n++;
              hitpoint.flux = (hitpoint.flux + hitpoint.adj.mul(adj))*g;
              stored = true;
            }
          }
        }
        double phi= 2. * M_PI * global_rand.genrand_real1();
        double r  = global_rand.genrand_real1();

        double x  = cos( phi ) * sqrt( 1.0 - r );
        double y  = sin( phi ) * sqrt( 1.0 - r );
        double z  = sqrt( r );

        __NS_RLR::Vector basis2 = nl;
        __NS_RLR::Vector basis1( 0.f, 0.f, 0.f );
        if( -0.6 < nl.x && nl.x < 0.6 )
          basis1.x = 1.;
        else if( -0.6 < nl.y && nl.y < 0.6 )
          basis1.y = 1.;
        else if( -0.6 < nl.z && nl.z < 0.6 )
          basis1.z = 1.;
        else
          basis1.x = 1.;

        __NS_RLR::Vector basis0 = basis1.cross( basis2 );
        basis0.normalize();
        basis1 = basis2.cross( basis0 );
        basis1.normalize();

        __NS_RLR::Vector dx = basis0 * x;
        __NS_RLR::Vector dy = basis1 * y;
        __NS_RLR::Vector dz = basis2 * z;
        __NS_RLR::Vector dir = dx + dy + dz;
        __NS_RLR::Vector hit2 = (dir * 0.00001) + hit;
        __NS_RLR::Ray ray2;
        ray2.pos_ = hit2;
        ray2.dir_ = dir;

        double p = stored ? MAX( color.x, MAX( color.y, color.z ) ) : 1.;

        if (global_rand.genrand_real1() < p)
          return trace(ray2,depth,mode,color.mul(adj)*(1./p),pix) | stored;
        else
          return stored ;
      }

    } else if( r < material->diffuse_ + material->reflect_ ) {

      // reflection
      __NS_RLR::Vector refl = nl * (nl.dot(ray.dir_)) * 2.;
      __NS_RLR::Vector dir = ray.dir_ - refl;
      __NS_RLR::Vector hit2 = (dir * 0.00001) + hit;
      __NS_RLR::Ray ray;
      ray.pos_ = hit2;
      ray.dir_ = dir;

      return trace( ray, depth, mode, color.mul( adj ), pix );

    } else if( r < material->diffuse_ + material->reflect_ + material->refract_ ) {

      __NS_RLR::Vector refl = nl * (nl.dot(ray.dir_)) * 2.;
      __NS_RLR::Vector dir = ray.dir_ - refl;
      __NS_RLR::Vector hit2 = (dir * 0.00001) + hit;
      __NS_RLR::Ray ray2;
      ray2.pos_ = hit2;
      ray2.dir_ = dir;

      bool into = norm.dot( nl ) > 0;                // Ray from outside going in?
      double nc = 1;
      double nt = material->ior_;
      double nnt= into ? nc/nt : nt/nc;
      double ddn= ray.dir_.dot(nl);
      double cos2t = 1. - nnt*nnt*(1. - ddn*ddn);
      
      if (cos2t < 0.)    // Total internal reflection
        return trace( ray2, depth, mode, color.mul( adj ), pix );

      __NS_RLR::Vector tdir = (ray.dir_*nnt - norm*((into?1.:-1.)*(ddn*nnt+sqrt(cos2t)))).normal();

      double a  = nt - nc;
      double b  = nt + nc;
      double R0 = a*a/(b*b);
      double c  = 1 - (into?-ddn:tdir.dot(norm));
      double Re = R0+(1-R0)*c*c*c*c*c;
      // double Tr = 1-Re;
      // double P  =.25+.5*Re;
      // double RP = Re/P;
      // double TP = Tr/(1-P);

      __NS_RLR::Ray ray3;
      __NS_RLR::Vector hit3 = (dir * -0.00001) + hit;
      ray3.pos_ = hit3;
      ray3.dir_ = tdir;

      float r2 = global_rand.genrand_real1();

      return (r2 < Re) ? trace( ray2, depth, mode, color.mul( adj ), pix ) : trace(ray3, depth, mode, color.mul( adj ), pix );
    }
    
    return false;
  }


  void eyepass( void ){
    int x, y;
    for(y=0;y<film_.height_;y++){
      int progress = y * 100 / film_.height_;
      printf("eyepass %3d%%\r",progress);
      fflush(stdout);
      for(x=0;x<film_.width_;x++){
        
        int pix = y * film_.width_ + x;
        
        double camera_x = (x / (double)film_.height_);
        double camera_y = (y / (double)film_.height_);
        __NS_RLR::Ray ray;
        if( camera_.ray( global_rand, camera_x, camera_y, ray ) )
          trace( ray, 0, 1, __NS_RLR::Color(1.,1.,1.), pix );
      }
    }
    printf("eyepass done\n");
    hitpoint_kd_.build();
  }
  
  void photonpass( int photons ){
    
    __NS_RLR::Scene::SceneObjectList  emitters( scene_.getEmitterList() );
    __NS_RLR::Scene::SceneObjectList::iterator it;
    if( emitters.size() == 0 ){
      printf("no emitter\n");
      return;
    }
    
    int p = 0;
    //for( int i = 0; i < photons; i ++ ) {
    for( p = 0; p < photons; p ++ ) {
      
      int e = global_rand.genrand_int31() % emitters.size();
      __NS_RLR::Vector light_power( emitters[e]->getMaterial()->color_ * emitters[e]->getMaterial()->emitter_ );
      __NS_RLR::Vector light_pos = emitters[e]->getRandomPoint( global_rand );
      __NS_RLR::Vector light_dir;
      if( p % 1023 == 0 ){
        int progress = p * 100 / photons;
        printf("light photon %3d%%\r",progress);
        fflush(stdout);
      }
      do{
        light_dir.x = global_rand.genrand_real1() * 2. - 1.;
        light_dir.y = global_rand.genrand_real1() * 2. - 1.;
        light_dir.z = global_rand.genrand_real1() * 2. - 1.;
      }while( light_dir.dot( light_dir ) > 1. );
      light_dir.normalize();
      
      __NS_RLR::Ray ray;
      ray.pos_ = light_pos + light_dir * 0.00001;
      ray.dir_ = light_dir;
      
      // trace( ray, 0, false, __NS_RLR::Color(1.,1.,1.), 0 );
      if( trace( ray, 0, false, light_power, 0 ) )
        p++;
      // storeされたphotonだけカウントする
    }
    total_photons_ += (double) p;
    printf("light photon done\n");
    eyepass();
  }
  
  void exposure( void ){

    printf("darkroom\n");
    // clean up.
    for( int i=0;i<film_.width_*film_.height_;i++)
      film_.film_[i] = __NS_RLR::Color(0.,0.,0.);
    
    max_r2_ = 0.;
    for( int i = 0; i < hitpoint_kd_.count();i++){
      HitPoint& hitpoint( hitpoint_kd_.get(i) );
      
      if( max_r2_ < hitpoint.r2 )
        max_r2_ = hitpoint.r2;
      int p = hitpoint.pix;
      film_.film_[p] = film_.film_[p] + hitpoint.flux * (1.0 / ( M_PI * hitpoint.r2 * total_photons_ ));
    }
    
    for( int i=0;i<film_.width_*film_.height_;i++)
      film_.film_[i] = film_.film_[i] * 10000.;
  }
};

template<typename ACCEL>
class RayTrace{
public:
  typedef RayTrace<ACCEL> self;
  __NS_RLR::Scene              &scene_;
  ACCEL                        &accel_;
  Camera             &camera_;
  __NS_RLR::FrameBuffer< self > film_;
  
  RayTrace(__NS_RLR::Scene &scene, ACCEL& accel, Camera& camera,int width, int height,int sample)
       :scene_(scene),accel_(accel),camera_(camera){
         film_.setup( width, height, sample );
       }
  
  inline double MAX(double a,double b){ return (a<b)?b:a; }

  __NS_RLR::Color directLight( __NS_RLR::Vector &pos, __NS_RLR::Vector &n ){
    __NS_RLR::Scene::SceneObjectList  emitters( scene_.getEmitterList() );
    __NS_RLR::Scene::SceneObjectList::iterator it;
    __NS_RLR::Color light(0.,0.,0.);
    
    for( it = emitters.begin(); it != emitters.end(); ++ it ) {
      __NS_RLR::Vector light_power( (*it)->getMaterial()->color_ * (*it)->getMaterial()->emitter_ );
      __NS_RLR::Vector light_pos  ( (*it)->getCenter() );
      __NS_RLR::Vector light_dir;
      __NS_RLR::Intersection isect;
      __NS_RLR::Ray          ray;
      
      ray.pos_ = pos + n * 0.0001;
      light_dir = light_pos - pos;
      double dist = light_dir.length() - 0.0001;
      ray.dir_ = light_dir.normal();
      double s    = ray.dir_.dot( n );
      if( s > 0 ){
        accel_.testRAY( ray, isect );
        if( !isect.hit() || isect.distance() >= dist )
          light = light + light_power * s;
      }
    }
    return light;
  }
  
  __NS_RLR::Color trace(const __NS_RLR::Ray &ray, int depth ) {
    depth++;
    __NS_RLR::Intersection isect;
    
    accel_.testRAY( ray, isect );
    
    if( !isect.hit() )
      return __NS_RLR::Color(0.,0.,0.);

    __NS_RLR::Vector hit  ( isect.getPoint()  );
    __NS_RLR::Vector norm ( isect.getNormal() );
    
    __NS_RLR::Vector nl   = (norm.dot(ray.direction())) < 0. ? norm : norm * -1.; // ray.dirと逆方向のnorm
    __NS_RLR::SceneObject   * object   = isect.getObject();
    __NS_RLR::SceneMaterial * material = object->getMaterial();
    __NS_RLR::Color           color( material->color_ );
    
    return color;
    
    double p = MAX( color.x, MAX( color.y, color.z ) );
    
    // russian roulette
    if( depth > 5 ){
      if( global_rand.genrand_real1() < p )
        color = color * (1./ p);
      else
        return __NS_RLR::Color( material->emitter_, material->emitter_, material->emitter_ );
    }

    float r = global_rand.genrand_real1();
    if( r < material->diffuse_ ){
      // diffuse.
      // double c = (norm.y + 1.)/2.;
      
      return color.mul( directLight( hit, norm ) )  + __NS_RLR::Color( material->emitter_, material->emitter_, material->emitter_ );
    } else if( r < material->diffuse_ + material->reflect_ ){
      // reflection
      __NS_RLR::Vector refl = nl * (nl.dot(ray.dir_)) * 2.;
      __NS_RLR::Vector dir = ray.dir_ - refl;
      __NS_RLR::Vector hit2 = (dir * 0.00001) + hit;
      __NS_RLR::Ray ray2;
      ray2.pos_ = hit2;
      ray2.dir_ = dir;
      return color.mul( trace( ray2, depth ) ) + __NS_RLR::Color( material->emitter_, material->emitter_, material->emitter_ );
    } else {
      __NS_RLR::Vector refl = nl * (nl.dot(ray.dir_)) * 2.;
      __NS_RLR::Vector dir = ray.dir_ - refl;
      __NS_RLR::Vector hit2 = (dir * 0.00001) + hit;
      __NS_RLR::Ray ray2;
      ray2.pos_ = hit2;
      ray2.dir_ = dir;
      
      bool into = norm.dot( nl ) > 0;                // Ray from outside going in?
      double nc = 1;
      double nt = 1.5;
      double nnt= into ? nc/nt : nt/nc;
      double ddn= ray.dir_.dot(nl);
      double cos2t = 1. - nnt*nnt*(1. - ddn*ddn);
      if (cos2t < 0.)    // Total internal reflection
        return color.mul( trace( ray2, depth ) ) + __NS_RLR::Color( material->emitter_, material->emitter_, material->emitter_ );
      
      __NS_RLR::Vector tdir = (ray.dir_*nnt - norm*((into?1.:-1.)*(ddn*nnt+sqrt(cos2t)))).normal();
      
      double a  = nt - nc;
      double b  = nt + nc;
      double R0 = a*a/(b*b);
      double c  = 1 - (into?-ddn:tdir.dot(norm));
      double Re = R0+(1-R0)*c*c*c*c*c;
      double Tr = 1-Re;
      double P  =.25+.5*Re;
      double RP = Re/P;
      double TP = Tr/(1-P);
      
      __NS_RLR::Ray ray3;
      __NS_RLR::Vector hit3 = (dir * -0.00001) + hit;
      ray3.pos_ = hit3;
      ray3.dir_ = tdir;
      
      float r2 = global_rand.genrand_real1();

      if( depth > 2 ){
        return __NS_RLR::Color( material->emitter_, material->emitter_, material->emitter_ )
          + color.mul( (r2 < P) ? trace(ray2,depth) * RP : trace(ray3,depth)*TP );
      }
      // else
      return __NS_RLR::Color( material->emitter_, material->emitter_, material->emitter_ )
        + color.mul( trace( ray2, depth ) * Re + trace( ray3, depth )*Tr );
    }
    return __NS_RLR::Color( material->emitter_, material->emitter_, material->emitter_ );
  }
  
  void eyepass(void) {
    ;
  }
  void photonpass( int photons ){(void)photons;}
  void exposure(){
    int x, y, s, progress, next; // ステージ座標
    __NS_RLR::Color rgb,sum;
    
    next = 0;
    for(y=0;y<film_.height_;y++){
      progress = y * 100 / film_.height_;
      printf("%3d%%\r",progress);
      fflush(stdout);
      for(x=0;x<film_.width_;x++){
        for(s=0;s<film_.sample_;s++){
          double camera_x = (x / (double)film_.height_);
          double camera_y = (y / (double)film_.height_);
          __NS_RLR::Ray ray;
          camera_.ray( global_rand, camera_x, camera_y, ray );
          rgb = trace( ray, 0 );
          film_.film_  [ y*film_.width_ + x ] = film_.film_[ y*film_.width_ + x] + rgb;
          film_.weight_[ y*film_.width_ + x ] += 1.;
        }
      }
    }
    printf("done\n");
    
  }
};

int main( int argc, char *argv[] )
{

  __NS_RLRUTIL::ModelHolder models;;
  __NS_RLR::Scene           scene;

  if( !models.load( argv[1] ) )
    return 0;
  
  models.registScene( scene );
  
  __NS_RLR::AABB aabb = scene.getAABB();
  aabb.dump();
  typedef __NS_BVH::BVH< __NS_RLR::Scene::SceneObjectList > ACCEL;
  typedef SPPM< ACCEL >                                     TRACE;
  //typedef PathTrace< ACCEL >                                TRACE;
  //typedef RayTrace< ACCEL >                               TRACE;
  //typedef DepthTrace< ACCEL >                             TRACE;
  //typedef NormTrace< ACCEL >                              TRACE;
  //typedef BVHDiag< ACCEL >                                TRACE;
  
  ACCEL bvh;
  bvh.build( scene.getObjectList() );
  
  Camera camera;
  NUMA::Lens lens;
  NUMA::LOADER::ZEMAX::load( argv[2], lens );
  lens.dump();
  camera.setup( VECTORMATH::Vector(100.f,100.f,-200.f),
                VECTORMATH::Vector(0.f,10.f,0.f),
                lens, lens.getLastDistance() + 146.f,
                lens.imageSize_ * 2.f,
                lens.imageSize_ * 2.f);
  TRACE trace( scene, bvh, camera, 1024, 1024, 1 );
  
  trace.eyepass();
  time_t t0,t1,t2;
  time(&t0);
  time(&t1);
  int pass=1,index = 1;
  while(1){
    printf("pass %d\n",pass);
    trace.photonpass( 5000000 );
    trace.exposure  ( );
    //trace.film_.saturate();
    trace.film_.normalize();
    time(&t2);
    double dt = difftime(t2,t1);
    printf("dt %f\n",dt);
    if( dt >= 60.f ) {
      char fn[256];
      sprintf(fn,"%s_%03d.png",argv[3],index);
      trace.film_.save_png(fn,1.,1.);
      // sprintf(fn,"%s_%03d.hdr",argv[3],index);
      // trace.film_.save_hdr(fn,1.,2.2);
      index ++;
      t1 = t2 - 60.f;
    }
    dt = difftime( t2, t0 );
    if( dt >= 60.f * 25.f ){
      char fn[256];
      sprintf(fn,"%s_result.png",argv[3]);
      trace.film_.save_png(fn,1.,1.);
      sprintf(fn,"%s_result.hdr",argv[3]);
      trace.film_.save_hdr(fn,1.,2.2);
      break;
    }
    pass++;
  }
  printf("end\n");
  return 0;
}
