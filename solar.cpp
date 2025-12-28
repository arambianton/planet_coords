#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <limits>
#include <cstdlib>
#include <algorithm>

namespace fs = std::filesystem;

long long days_from_civil(int y,unsigned m,unsigned d){
    y-=m<=2; int era=(y>=0?y:y-399)/400;
    unsigned yoe=unsigned(y-era*400);
    unsigned doy=(153*(m+(m>2?-3:9))+2)/5+d-1;
    unsigned doe=yoe*365+yoe/4-yoe/100+doy;
    return 1LL*era*146097 + doe - 719468;
}
long long seconds_utc(int y,int m,int d,int H,int M,int S){
    return days_from_civil(y,m,d)*86400LL + H*3600LL + M*60LL + S;
}

struct Vec3{
    double x,y,z;
    Vec3(double X=0, double Y=0, double Z=0):x(X),y(Y),z(Z){}
    Vec3& operator+=(const Vec3& o){ x+=o.x; y+=o.y; z+=o.z; return *this; }
};
inline Vec3 operator-(const Vec3&a, const Vec3&b){ return {a.x-b.x,a.y-b.y,a.z-b.z}; }
inline Vec3 operator*(double k, const Vec3&v){ return {k*v.x,k*v.y,k*v.z}; }

class SolarBody{
    double m,qx,qy,qz,px,py,pz; std::string name;
public:
    SolarBody(double M,double X,double Y,double Z,
              double PX,double PY,double PZ,const std::string& n):
        m(M),qx(X),qy(Y),qz(Z),px(PX),py(PY),pz(PZ),name(n){}
    SolarBody():m(1),qx(0),qy(0),qz(0),px(0),py(0),pz(0),name(""){}
    double get_mass()const{return m;}
    Vec3   get_q() const{return {qx,qy,qz};}
    Vec3   get_p() const{return {px,py,pz};}
    void   set_q(const Vec3&v){ qx=v.x; qy=v.y; qz=v.z; }
    void   set_p(const Vec3&v){ px=v.x; py=v.y; pz=v.z; }
    const std::string& get_name()const{return name;}
};

class SolarSystem{
    std::vector<SolarBody> b;
    static constexpr double G = 6.67430e-11;
public:
    void push(const SolarBody& s){ b.push_back(s); }
    std::size_t size()const{return b.size();}
    SolarBody& body(std::size_t i){ return b[i]; }
    const SolarBody& body(std::size_t i)const{ return b[i]; }

    std::vector<Vec3> forces()const{
        std::size_t n = b.size(); std::vector<Vec3> F(n);
        for(std::size_t i=0;i<n;++i){
            Vec3 qi=b[i].get_q(); double mi=b[i].get_mass();
            for(std::size_t j=i+1;j<n;++j){
                Vec3 qj = b[j].get_q(); 
                double mj = b[j].get_mass();
                Vec3 dr = qi - qj; 
                double r2 = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;
                double inv_r3 = 1.0 / (r2*std::sqrt(r2));
                Vec3 Fij = -G*mi*mj*inv_r3*dr;
                F[i] += Fij; F[j] += Vec3(-Fij.x,-Fij.y,-Fij.z);
            }
        }
        return F;
    }

    void step(double h){
        auto F=forces();
        for(std::size_t i=0;i<size();++i){
            Vec3 p=b[i].get_p(); p+=0.5*h*F[i]; b[i].set_p(p);
        }
        for(std::size_t i=0;i<size();++i){
            Vec3 q=b[i].get_q(); Vec3 p=b[i].get_p();
            q+=(h/b[i].get_mass())*p; b[i].set_q(q);
        }
        F=forces();
        for(std::size_t i=0;i<size();++i){
            Vec3 p=b[i].get_p(); p+=0.5*h*F[i]; b[i].set_p(p);
        }
    }

    double energy()const{
        double Ek=0,Ep=0; std::size_t n=b.size();
        for(std::size_t i=0;i<n;++i){
            Vec3 p=b[i].get_p(); double m=b[i].get_mass();
            Ek+=(p.x*p.x+p.y*p.y+p.z*p.z)/(2*m);
            for(std::size_t j=i+1;j<n;++j){
                Vec3 dr=b[i].get_q()-b[j].get_q();
                double r=std::sqrt(dr.x*dr.x+dr.y*dr.y+dr.z*dr.z);
                Ep-=G*m*b[j].get_mass()/r;
            }
        }
        return Ek+Ep;
    }

    void integrate_interval(double dt, double& h, double epsrel, double h_min=1.0){
        double sign = (dt>=0)? 1.0 : -1.0;
        double t_left = std::fabs(dt);
        while(t_left > 1e-9){
            if(h > t_left) h = t_left;                   
            bool ok=false;
            while(!ok){
                double Eold = energy();
                std::vector<SolarBody> backup = b;      
                step(sign*h);
                double Enew = energy();
                double rel = std::fabs(Enew-Eold)/std::fabs(Eold);
                if(rel <= epsrel || h<=h_min){
                    ok=true; 
                }else{
                    b = std::move(backup);
                    h *= 0.5;
                }
            }
            t_left -= h;
        }
    }
};

SolarBody parseBody(const fs::path& f){
    std::ifstream file(f); if(!file) throw std::runtime_error("no "+f.string());
    double M=1,X=0,Y=0,Z=0,VX=0,VY=0,VZ=0; std::string name="";
    std::string line;
    while(std::getline(file,line)){
        if(line.find("Mass")!=std::string::npos){
            auto pos=line.find_last_of(" ="); M=std::stod(line.substr(pos+1));
        }
        if(line.find("X =")!=std::string::npos){
            auto px=line.find("X ="),py=line.find("Y =",px),pz=line.find("Z =",py);
            X=1e3*std::stod(line.substr(px+3));
            Y=1e3*std::stod(line.substr(py+3));
            Z=1e3*std::stod(line.substr(pz+3));
        }
        if(line.find("VX=")!=std::string::npos){
            auto px=line.find("VX="),py=line.find("VY=",px),pz=line.find("VZ=",py);
            VX=1e3*std::stod(line.substr(px+3));
            VY=1e3*std::stod(line.substr(py+3));
            VZ=1e3*std::stod(line.substr(pz+3));
        }
        if(line.find("Revised")!=std::string::npos){
            std::istringstream ss(line); std::string tmp; ss>>tmp>>tmp>>tmp>>tmp>>name;
        }
    }
    return SolarBody(M,X,Y,Z,M*VX,M*VY,M*VZ,name);
}

/* ==============================   main   ============================== */
int main(int argc,char*argv[]){
    if(argc<4){
        std::cerr<<"usage: "<<argv[0]<<" <dir> <DD.MM.YYYY> <HH:MM> [h0] [EPSREL]\n";
        return 1;
    }
    fs::path dir = argv[1];
    std::string dateStr=argv[2], timeStr=argv[3];
    double h = (argc>=5)? std::stod(argv[4]) : 600.0;   // начальный шаг, с
    double EPSREL = (argc>=6)? std::stod(argv[5]) : 1e-8;

    /* --- читаем планеты --- */
    SolarSystem sys;
    for(auto& f: fs::recursive_directory_iterator(dir))
        if(f.is_regular_file()) sys.push(parseBody(f.path()));
    if(sys.size()==0){ std::cerr<<"no planet files\n"; return 1; }

    /* --- парсим дату --- */
    int d,m,y,H,M; char c1,c2,c3;
    if(std::sscanf(dateStr.c_str(),"%d%c%d%c%d",&d,&c1,&m,&c2,&y)!=5||c1!='.'||c2!='.'){
        std::cerr<<"bad date\n"; return 1;}
    if(std::sscanf(timeStr.c_str(),"%d%c%d",&H,&c3,&M)!=3||c3!=':'){
        std::cerr<<"bad time\n"; return 1;}

    const long long t0 = seconds_utc(2025,2,14,0,0,0);
    const long long tgt= seconds_utc(y,m,d,H,M,0);
    long long dtSec = tgt - t0;

    std::cout<<"EPSREL = "<<EPSREL<<"\n";

    /* --- интегрируем с адаптацией --- */
    sys.integrate_interval(static_cast<double>(dtSec), h, EPSREL);
    std::cout<<"Final h ≈ "<<h<<" s\n";

    /* --- вывод координат --- */
    const bool PRINT_AU=false; const double AU=1.495978707e+8;
    auto val=[&](double km){return PRINT_AU?km/AU:km;};
    const char* unit=PRINT_AU?"AU":"km";
    constexpr int PREC=8; std::cout<<std::scientific<<std::uppercase<<std::setprecision(PREC);
    const int W = 15; 
    const char SEP = ' ';

    std::cout<<"\nCoordinates for "<<dateStr<<' '<<timeStr<<"\n\n"
             <<std::left<<std::setw(10)<<"Body"
             <<SEP<<std::right<<std::setw(W)<<('X' + std::string{" ["}+unit+']')
             <<SEP<<std::setw(W)<<('Y' + std::string{" ["}+unit+']')
             <<SEP<<std::setw(W)<<('Z' + std::string{" ["}+unit+']')
             <<SEP<<std::setw(W)<<"VX [km/s]"
             <<SEP<<std::setw(W)<<"VY [km/s]"
             <<SEP<<std::setw(W)<<"VZ [km/s]\n"
             <<std::string(10+6*W,'-')<<'\n';

    for(std::size_t i=0;i<sys.size();++i){
        const auto& B=sys.body(i);
        Vec3 q=B.get_q(), p=B.get_p(); double m=B.get_mass();
        double X=q.x/1e3,Y=q.y/1e3,Z=q.z/1e3;
        double VX=(p.x/m)/1e3,VY=(p.y/m)/1e3,VZ=(p.z/m)/1e3;
        std::cout<<std::left<<std::setw(10)<<B.get_name()
                 <<SEP<<std::right<<std::setw(W)<<val(X)
                 <<SEP<<std::setw(W)<<val(Y)
                 <<SEP<<std::setw(W)<<val(Z)
                 <<SEP<<std::setw(W)<<VX
                 <<SEP<<std::setw(W)<<VY
                 <<SEP<<std::setw(W)<<VZ<<'\n';
    }
    return 0;
}
