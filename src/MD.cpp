#include "math.h"
#include "MD.hpp"
#include "config.hpp"
#include "molecule.hpp"
#include "randoms.hpp"
#include "utils.hpp"
#include "algebra.hpp"

Timer timer;

void Initialization(){
    timer.tik();
    stepCount = 0;
    ConstructMol();
    InitCoords();
    InitVels();
    InitAccels();
    InitAngCoords();
    InitAngVels();
    InitAngAccels();
    mkdir_fs(DataPath);
    EvalProps();
    save(false);
    printInfo();
}

void SingleStep(){
    stepCount++;
    PredictorStep();
    PredictorStepQ();
    GenSiteCoords();
    ComputeForces();
    ComputeTorqs();
    ComputeQa();
    ApplyThermostat();
    CorrectorStep();
    CorrectorStepQ();
    AdjustQuat();
    ApplyBoundaryCond();
    EvalProps();
    save(true);
    if(stepCount%stepAdjustTemp==0) AdjustTemp();
    if(stepCount%stepPrintInfo==0) printInfo();
}

bool NotFinished(){return stepCount<tstepNum;}

RMat ComputeInert(){
    double mass = 1.0/MassSites.size();
    RMat IMat; IMat.fill(0.);
    for(auto msite:MassSites){
        IMat[0] += msite[1]*msite[1] + msite[2]*msite[2];
        IMat[4] += msite[2]*msite[2] + msite[0]*msite[0];
        IMat[8] += msite[0]*msite[0] + msite[1]*msite[1];
        IMat[3] -= msite[0]*msite[1];
        IMat[6] -= msite[0]*msite[2];
        IMat[7] -= msite[1]*msite[2];
    }
    IMat[1] = IMat[3];
    IMat[2] = IMat[6];
    IMat[5] = IMat[7];
    IMat = IMat * mass;
    return IMat;
}
void ConstructMol(){
    Vec3d MassCenter = {0., 0., 0.};
    for(auto v:MassSites) MassCenter = MassCenter + v;
    MassCenter = MassCenter/double(MassSites.size());
    for(auto& v:LJSites) v = v - MassCenter;
    for(auto& v:MassSites) v = v - MassCenter;
    auto IMat = ComputeInert();
    if(!isDiag(IMat)){std::cout<<"IMat is not diagonal!\n"; exit(1);}
    mInert = getDiag(IMat);
}
void InitCoords(){
    std::array<int,3> ns;
    for(int i=0; i<3; i++)ns[i] = int(region[i]/ucell[i]);
    auto nmax = ns[0]*ns[1]*ns[2];
    if(nmax<MolNum){std::cout<<"Box region is too small!\n";exit(1);}
    int count = 0;
    for(int nz=0; nz<ns[2]; nz++){
        for(int ny=0; ny<ns[1]; ny++){
            for(int nx=0; nx<ns[0]; nx++){
                if(count<MolNum){
                    Vec3d c = {nx+0.5, ny+0.5, nz+0.5};
                    c = ElProd(c,ucell);
                    c = c - region*0.5;
                    mols[count].r = c;
                }
                count++;
            }
        }
    }
}

void InitVels(){
    Vec3d vSum = {0., 0., 0.};
    for(int n=0; n<MolNum; n++){
        mols[n].rv = VRand() * velMag;
        vSum = vSum + mols[n].rv;
    }
    vSum = vSum/double(MolNum);
    #pragma omp parallel for
    for(int n=0; n<MolNum; n++) mols[n].rv = mols[n].rv-vSum;
}

void InitAccels(){
    #pragma omp parallel for
    for(int n=0; n<MolNum; n++){
        mols[n].ra.fill(0.);
        mols[n].ra1.fill(0.);
        mols[n].ra2.fill(0.);
    }
}

void InitAngCoords(){
    for(int n=0; n<MolNum; n++){
        auto e = VRand();
        Vec3d eAng;
        eAng[0] = std::atan2(e[0],e[1]);
        eAng[1] = std::acos(e[2]);
        eAng[2] = 2. * M_PI * RandR();
        mols[n].q = EulerToQuat(eAng);
    }
}

Vec4d EulerToQuat(const Vec3d& ang){
    auto a1 = 0.5 * ang[1];
    auto a2 = 0.5 * (ang[0]-ang[2]);
    auto a3 = 0.5 * (ang[0]+ang[2]);
    Vec4d q = {std::sin(a1)*std::cos(a2), std::sin(a1)*std::sin(a2), std::cos(a1)*std::sin(a3), std::cos(a1)*std::cos(a3)};
    return q;
}

void InitAngVels(){
    for(int n=0; n<MolNum; n++){
        auto e = VRand();
        Vec4d qe = {e[0], e[1], e[2], 0.};
        mols[n].qv = QuatMul(mols[n].q, qe);
        auto f = velMag/std::sqrt(mInert*ElProd(e,e));
        mols[n].qv = mols[n].qv * f;
    }
}

void InitAngAccels(){
    #pragma omp parallel for
    for(int n=0; n<MolNum; n++){
        mols[n].qa.fill(0.);
        mols[n].qa1.fill(0.);
        mols[n].qa2.fill(0.);
    }
}

double PairPotential(Vec3d r1, Vec3d r2){
    auto r = norm(r1-r2);
    auto v = std::pow(1./r,6);
    return 4.*(v*v-v);
}

Vec3d PairForce(Vec3d r1, Vec3d r2){
    auto dr = r1 - r2;
    auto R2 = 1.0/(dr * dr);
    auto R6 = std::pow(R2,3);
    auto factor = 48.0 * R2 * R6 * (R6 - 0.5);
    return dr*factor;
}

Vec3d ComputeW(int n){
    Vec4d qt, qvt;
    qvt = mols.at(n).qv;
    qvt[3] *= -1.0;
    qt = QuatMul(qvt, mols[n].q) * 2.0;
    Vec3d w = {qt[0], qt[1], qt[2]};
    return w;
}

void ComputeQa(){
    #pragma omp parallel for
    for(int n = 0; n < MolNum; n++){
        Vec4d qs;
        Vec3d w;
        w = ComputeW(n);
        qs[0] = (mols[n].torq[0] + (mInert[1]-mInert[2])*w[1]*w[2])/mInert[0];
        qs[1] = (mols[n].torq[1] + (mInert[2]-mInert[0])*w[2]*w[0])/mInert[1];
        qs[2] = (mols[n].torq[2] + (mInert[0]-mInert[1])*w[0]*w[1])/mInert[2];
        qs[3] = mols[n].qv * mols[n].qv * (-2.0);
        mols[n].qa = QuatMul(mols[n].q, qs) * 0.5;
    }
}

RMat BuildRotMat(const Vec4d& q, bool transpose){
    double p[10], s;
    int k=0,k1,k2;
    for(k2 = 0; k2 < 4; k2++){
        for(k1 = k2; k1 < 4; k1++,k++)p[k] = 2.0 * q[k1]*q[k2];
    }
    RMat R;
    R[0] = p[0] + p[9] - 1.0;
    R[4] = p[4] + p[9] - 1.0;
    R[8] = p[7] + p[9] - 1.0;
    s = transpose?1.0:-1.0;
    R[1] = p[1] + s * p[8];
    R[3] = p[1] - s * p[8];
    R[2] = p[2] - s * p[6];
    R[6] = p[2] + s * p[6];
    R[5] = p[5] + s * p[3];
    R[7] = p[5] - s * p[3];
    return R;
}

void GenSiteCoords(){
    #pragma omp parallel for
    for(int n=0; n<MolNum; n++){
        RMat R;
        Vec3d t;
        R = BuildRotMat(mols[n].q, true);
        for(int j=0; j<siteNum; j++){
            t = MV(R, LJSites[j]);
            fsites[n*siteNum+j].r = mols[n].r + t;
        }
    }
}

void ComputeForces(){
    int totNum = siteNum * MolNum;
    #pragma omp parallel for
    for(int i=0; i<totNum; i++)fsites[i].f.fill(0.);
    for(int m1=0; m1<MolNum-1; m1++){
        for(int m2=m1+1; m2<MolNum; m2++){
            auto dr = mols[m1].r - mols[m2].r;
            if(dr*dr<rc2){
                auto ms1 = m1 * siteNum;
                auto ms2 = m2 * siteNum;
                for(int j1=0; j1<siteNum; j1++){
                    auto idx1 = ms1+j1;
                    for(int j2=0; j2<siteNum; j2++){
                        auto idx2 = ms2+j2;
                        auto f = PairForce(fsites[idx1].r,fsites[idx2].r);
                        fsites[idx1].f = fsites[idx1].f + f;
                        fsites[idx2].f = fsites[idx2].f - f;
                    }
                }
            }
        }
    }
}

void ComputeTorqs(){
    #pragma omp parallel for
    for(int n=0; n<MolNum; n++){
        RMat R;
        Vec3d dr, t, torqS;
        mols[n].ra.fill(0.0);
        torqS.fill(0.0);
        for(int j=0; j<siteNum; j++){
            auto sidx = n * siteNum + j;
            mols[n].ra = mols[n].ra + fsites[sidx].f;
            dr = fsites[sidx].r - mols[n].r;
            t = Cross(dr, fsites[sidx].f);
            torqS = torqS + t;
        }
        R = BuildRotMat(mols[n].q, false);
        mols[n].torq = MV(R,torqS);
    }
}
// predictor-corrector integration with k=4
template <size_t N>
void PCR4(arr<N>& r, arr<N>& r0, arr<N>& v, arr<N>& a, arr<N>& a1, arr<N>& a2, const Vec3d& c){
    r = r0 + v*dt + (a*c[0] + a1*c[1] + a2*c[2])*dt2;
}

template <size_t N>
void PCV4(arr<N>& r, arr<N>& r0, arr<N>& v, arr<N>& a, arr<N>& a1, arr<N>& a2, const Vec3d& c){
    v = (r-r0)/dt + (a*c[0] + a1*c[1] + a2*c[2])*dt;
}
void PR(int n){
    PCR4(mols[n].r, mols[n].r, mols[n].rv, mols[n].ra, mols[n].ra1, mols[n].ra2, pr);
}
void PRV(int n){
    PCV4(mols[n].r, mols[n].r0, mols[n].rv, mols[n].ra, mols[n].ra1, mols[n].ra2, pv);
}
void CR(int n){
    PCR4(mols[n].r, mols[n].r0, mols[n].rv0, mols[n].ra, mols[n].ra1, mols[n].ra2, cr);
}
void CRV(int n){
    PCV4(mols[n].r, mols[n].r0, mols[n].rv, mols[n].ra, mols[n].ra1, mols[n].ra2, cv);
}
void PQ(int n){
    PCR4(mols[n].q, mols[n].q, mols[n].qv, mols[n].qa, mols[n].qa1, mols[n].qa2, pr);
}
void PQV(int n){
    PCV4(mols[n].q, mols[n].q0, mols[n].qv, mols[n].qa, mols[n].qa1, mols[n].qa2, pv);
}
void CQ(int n){
    PCR4(mols[n].q, mols[n].q0, mols[n].qv0, mols[n].qa, mols[n].qa1, mols[n].qa2, cr);
}
void CQV(int n){
    PCV4(mols[n].q, mols[n].q0, mols[n].qv, mols[n].qa, mols[n].qa1, mols[n].qa2, cv);
}

void PredictorStep(){
    #pragma omp parallel for
    for(int n=0; n < MolNum; n++){
        mols[n].r0 = mols[n].r;
        mols[n].rv0 = mols[n].rv;
        PR(n); 
        PRV(n);
        mols[n].ra2 = mols[n].ra1;
        mols[n].ra1 = mols[n].ra;
    }
}

void CorrectorStep(){
    #pragma omp parallel for
    for(int n=0; n<MolNum; n++){
        CR(n);
        CRV(n);
    }
}

void PredictorStepQ(){
    #pragma omp parallel for
    for(int n=0; n<MolNum; n++){
        mols[n].q0 = mols[n].q;
        mols[n].qv0 = mols[n].qv;
        PQ(n);
        PQV(n);
        mols[n].qa2 = mols[n].qa1;
        mols[n].qa1 = mols[n].qa;
    }
}

void CorrectorStepQ(){
    #pragma omp parallel for
    for(int n=0; n<MolNum; n++){
        CQ(n);
        CQV(n);
    }
}

void ApplyBoundaryCond(){
    #pragma omp parallel for
    for(int n=0; n<MolNum; n++){
        for(int i=0; i<3; i++){
            if(mols[n].r[i]>=0.5*region[i])mols[n].r[i]-=region[i];
            else if(mols[n].r[i]<-0.5*region[i])mols[n].r[i]+=region[i];
        }
    }
}
void AdjustQuat(){
    #pragma omp parallel for
    for(int n = 0; n < MolNum; n++){
        auto qnorm = norm(mols[n].q);
        mols[n].q = mols[n].q/qnorm;
    }
}

double ComputeTK(){
    double Esum = 0.0;
    #pragma omp parallel for reduction(+:Esum)
    for(int n = 0; n < MolNum; n++) Esum += mols[n].rv*mols[n].rv;
    return 0.5*Esum;
}

double ComputeRK(){
    double Esum = 0.0;
    Vec3d w;
    #pragma omp parallel for reduction(+:Esum)
    for(int n = 0; n < MolNum; n++){
        w = ComputeW(n);
        Esum += mInert*ElProd(w,w);
    }
    return 0.5*Esum;
}

double ComputeV(){
    double vSum = 0.0;
    #pragma omp parallel for reduction(+:vSum)
    for(int m1=0; m1<MolNum-1; m1++){
        for(int m2=m1+1; m2<MolNum; m2++){
            auto dr = mols[m1].r - mols[m2].r;
            if(dr*dr<rc2){
                auto ms1 = m1 * siteNum;
                auto ms2 = m2 * siteNum;
                for(int j1=0; j1<siteNum; j1++){
                    auto idx1 = ms1+j1;
                    for(int j2=0; j2<siteNum; j2++){
                        auto idx2 = ms2+j2;
                        vSum += PairPotential(fsites[idx1].r, fsites[idx2].r);
                    }
                }
            }
        }
    }
    return vSum;
}

// Temperature Control
void ApplyThermostat(){
    double s1{0.0}, s2{0.0}, vFac;
    Vec3d w;
    #pragma omp parallel for reduction(+:s1,s2)
    for(int n=0; n<MolNum; n++){
        s1 += mols[n].rv * mols[n].ra;
        s2 += mols[n].rv * mols[n].rv;
        w = ComputeW(n);
        s1 += w * mols[n].torq;
        s2 += mInert * ElProd(w,w);
    }
    vFac = -s1/s2;
    #pragma omp parallel for
    for(int n=0; n<MolNum; n++){
        mols[n].ra = mols[n].ra + mols[n].rv * vFac;
        mols[n].qa = mols[n].qa + mols[n].qv * vFac;
    }
}

void AdjustTemp(){
    double ETK = ComputeTK();
    double vFac = velMag/std::sqrt(ETK/MolNum);
    double ERK = ComputeRK();
    double qvFac = velMag/std::sqrt(ERK/MolNum);
    #pragma omp parallel for
    for(int n=0; n<MolNum; n++){
        mols[n].rv = mols[n].rv * vFac;
        mols[n].qv = mols[n].qv * qvFac;
    }
}

void EvalProps(){
    ETransK = ComputeTK()/MolNum;
    ERotK = ComputeRK()/MolNum;
    EPot = ComputeV()/MolNum;
}

void save(bool isapp){
    std::ofstream outfile;
    save<double>(&ETransK,1,&outfile,ETKFile,isapp);
    save<double>(&ERotK,1,&outfile,ERKFile,isapp);
    save<double>(&EPot,1,&outfile,EVFile,isapp);
    for(auto site:fsites){
        save<double>(site.r.data(),3,&outfile,PosFile,isapp);
    }
}

void printInfo(){
    std::cout<<"\nStep:"<<stepCount<<", Time Passed:"<<timer.elapse()<<"ms"\
        <<"\nETransK:"<<ETransK<<", ERotK:"<<ERotK<<", EPot:"<<EPot<<"\n";
}