#include "micro_demag_pthread.hpp"
#include "../misc.hpp"


DemagFieldMultithread::DemagFieldMultithread (Mesh meshin, Material paramin, bool verbose, bool caching, int nthreads) : material(paramin),mesh(meshin), nthreads(nthreads){
    af::timer demagtimer = af::timer::start();
    if (caching == false){
        Nfft=N_cpp_alloc(mesh.n0_exp,mesh.n1_exp,mesh.n2_exp,mesh.dx,mesh.dy,mesh.dz);
        if (verbose) printf("time demag init [af-s]: %f\n", af::timer::stop(demagtimer));
    }
    else{
        std::string magafdir = setup_magafdir();
        unsigned long long maxsize = 2000000;
        unsigned long long reducedsize = 1000000;
        std::string nfft_id = "n0exp_"+std::to_string(mesh.n0_exp)+"_n1exp_"+std::to_string(mesh.n1_exp)+"_n2exp_"+std::to_string(mesh.n2_exp)+"_dx_"+std::to_string(1e9*mesh.dx)+"_dy_"+std::to_string(1e9*mesh.dy)+"_dz_"+std::to_string(1e9*mesh.dz);
        std::string path_to_nfft_cached = magafdir+nfft_id;
        int checkarray=-1;
        if (exists(path_to_nfft_cached)){
            try{
                checkarray = af::readArrayCheck(path_to_nfft_cached.c_str(), "");
            }
            catch (const af::exception& e){
                printf("Warning, af::readArrayCheck failed, omit reading array.\n%s\n", e.what());
                //std::cout << "Warning, af::readArrayCheck failed, omit reading array.\n"<< e.what() << std::endl;
            }
        }
        if(checkarray > -1){
            if (verbose) printf("Reading demag tensor from '%s'\n", path_to_nfft_cached.c_str());
            //if (verbose) std::cout << "Reading demag tensor from '"<< path_to_nfft_cached << "'." << std::endl;
            Nfft = af::readArray(path_to_nfft_cached.c_str(), "");
        }
        else{
            Nfft=N_cpp_alloc(mesh.n0_exp,mesh.n1_exp,mesh.n2_exp,mesh.dx,mesh.dy,mesh.dz);
            unsigned long long magafdir_size_in_bytes = GetDirSize(magafdir);
            //if (verbose) std::cout << "current size of ~/.magnum.af.cache = " << magafdir_size_in_bytes << " bytes." << std::endl;
            if (magafdir_size_in_bytes > maxsize){
                if (verbose) printf("Maintainance: size of '%s' is %f GB > %f GB, removing oldest files until size < %f GB\n", magafdir.c_str(), (double)magafdir_size_in_bytes/1e6, (double)maxsize/1e6, (double)reducedsize/1e6);
                //if (verbose) std::cout << "Warning: ~/.magnum.af.cache is larger than 1GB, omitting to save demag tensor" << std::endl;
                remove_oldest_files_until_size(magafdir.c_str(), reducedsize, verbose);
                if (verbose) printf("Maintainance finished: '%s' has now %f GB\n", magafdir.c_str(), (double)GetDirSize(magafdir)/1e6);
            }
            if (GetDirSize(magafdir) < maxsize){
                try{
                    if (verbose) printf("Saving demag tensor to'%s'\n", path_to_nfft_cached.c_str());
                    //if (verbose) std::cout << "Saving demag tensor to'"<< path_to_nfft_cached << "'." << std::endl;
                    af::saveArray("", Nfft, path_to_nfft_cached.c_str());
                    if (verbose) printf("Saved demag tensor to'%s'\n", path_to_nfft_cached.c_str());
                }
                catch (const af::exception& e){
                    printf("Warning, af::saveArray failed, omit saving demag tensor.\n%s\n", e.what());
                    //std::cout << "Warning, af::saveArray failed, omit saving demag tensor.\n"<< e.what() << std::endl;
                }
            }
        }
        if (verbose) printf("time demag init [af-s]: %f\n", af::timer::stop(demagtimer));
        //if (verbose) std::cout<<"time demag init [af-s]: "<< af::timer::stop(demagtimer) <<std::endl;
    }
}

//Energy calculation
//Edemag=-mu0/2 integral(M . Hdemag) dx
double DemagFieldMultithread::E(const State& state){
    return -constants::mu0/2. * material.ms * afvalue(sum(sum(sum(sum(h(state)*state.m,0),1),2),3)) * mesh.dx * mesh.dy * mesh.dz;
}

double DemagFieldMultithread::E(const State& state, const af::array& h){
    return -constants::mu0/2. * material.ms * afvalue(sum(sum(sum(sum(h * state.m,0),1),2),3)) * mesh.dx * mesh.dy * mesh.dz;
}

void DemagFieldMultithread::print_Nfft(){
    af::print("Nfft=", Nfft);
}


af::array DemagFieldMultithread::h(const State&  state){
    timer_demagsolve = af::timer::start();
    // FFT with zero-padding of the m field
    af::array mfft;
    if (mesh.n2_exp == 1){
        if (state.Ms.isempty()) mfft=af::fftR2C<2>(material.ms * state.m,af::dim4(mesh.n0_exp,mesh.n1_exp));
        else mfft=af::fftR2C<2>(state.Ms * state.m,af::dim4(mesh.n0_exp,mesh.n1_exp));
    }
    else {
        if (state.Ms.isempty()) mfft=af::fftR2C<3>(material.ms * state.m,af::dim4(mesh.n0_exp,mesh.n1_exp,mesh.n2_exp));
        else  mfft=af::fftR2C<3>(state.Ms * state.m,af::dim4(mesh.n0_exp,mesh.n1_exp,mesh.n2_exp));
    }
  
    // Pointwise product
    af::array hfft=af::array (mesh.n0_exp/2+1,mesh.n1_exp,mesh.n2_exp,3,c64);
    hfft(af::span,af::span,af::span,0)= Nfft(af::span,af::span,af::span,0) * mfft(af::span,af::span,af::span,0) 
                                   + Nfft(af::span,af::span,af::span,1) * mfft(af::span,af::span,af::span,1)
                                   + Nfft(af::span,af::span,af::span,2) * mfft(af::span,af::span,af::span,2);
    hfft(af::span,af::span,af::span,1)= Nfft(af::span,af::span,af::span,1) * mfft(af::span,af::span,af::span,0) 
                                   + Nfft(af::span,af::span,af::span,3) * mfft(af::span,af::span,af::span,1)
                                   + Nfft(af::span,af::span,af::span,4) * mfft(af::span,af::span,af::span,2);
    hfft(af::span,af::span,af::span,2)= Nfft(af::span,af::span,af::span,2) * mfft(af::span,af::span,af::span,0)
                                   + Nfft(af::span,af::span,af::span,4) * mfft(af::span,af::span,af::span,1)
                                   + Nfft(af::span,af::span,af::span,5) * mfft(af::span,af::span,af::span,2);
  
    // IFFT reversing padding
    af::array h_field;
    if (mesh.n2_exp == 1){
        h_field=af::fftC2R<2>(hfft);
        if(material.afsync) af::sync();
        cpu_time += af::timer::stop(timer_demagsolve);
        return h_field(af::seq(0,mesh.n0_exp/2-1),af::seq(0,mesh.n1_exp/2-1));
    }
    else {
        h_field=af::fftC2R<3>(hfft);
        if(material.afsync) af::sync();
        cpu_time += af::timer::stop(timer_demagsolve);
        return h_field(af::seq(0,mesh.n0_exp/2-1),af::seq(0,mesh.n1_exp/2-1),af::seq(0,mesh.n2_exp/2-1),af::span);
    }
}

struct LoopInfo {
    LoopInfo(){}
    LoopInfo(int i0_start, int i0_end, int n0_exp, int n1_exp, int n2_exp,  double dx,  double dy,  double dz):
        i0_start(i0_start), i0_end(i0_end), n0_exp(n0_exp), n1_exp(n1_exp), n2_exp(n2_exp), dx(dx), dy(dy), dz(dz){}
    int i0_start;
    int i0_end;
    int n0_exp;
    int n1_exp;
    int n2_exp;
    double dx;
    double dy;
    double dz;
};

//https://stackoverflow.com/questions/22657770/using-c-11-multithreading-on-non-static-member-function

double* N = NULL;
//int step_i = 0;

//void DemagFieldMultithread::setup_N(Mesh mesh) 
//void* DemagFieldMultithread::setup_N(void* arg) 

double thread_newellg(double x, double y, double z);
double thread_newellf(double x, double y, double z);
double thread_Nxxg(int ix, int iy, int iz, double dx, double dy, double dz);
double thread_Nxxf(int ix, int iy, int iz, double dx, double dy, double dz);

void* setup_N(void* arg) 
{ 
    //LoopInfo* lptr= arg;
    //LoopInfo loopinfo = *lptr;
    
    //LoopInfo* loopinfo = static_cast<LoopInfo>(arg);
    LoopInfo* loopinfo = static_cast<LoopInfo*>(arg);
    //std::cout << "called setup_N with loopinfo->>n0_exp " << loopinfo->n0_exp << std::endl;
    //std::cout << "called setup_N with loopinfo->>n0_exp " << loopinfo->n1_exp << std::endl;
    //std::cout << "called setup_N with loopinfo->>n0_exp " << loopinfo->n2_exp << std::endl;
    //std::cout << "called setup_N with loopinfo->>n0_exp " << loopinfo->dx << std::endl;
    //std::cout << "called setup_N with loopinfo->>n0_exp " << loopinfo->dy << std::endl;
    //std::cout << "called setup_N with loopinfo->>n0_exp " << loopinfo->dz << std::endl;
    //int core = step_i++; 

    for (int i0 = loopinfo->i0_start; i0 < loopinfo->i0_end; i0++){
        const int j0 = (i0 + loopinfo->n0_exp/2) % loopinfo->n0_exp - loopinfo->n0_exp/2;
        for (int i1 = 0; i1 < loopinfo->n1_exp; i1++){
            const int j1 = (i1 + loopinfo->n1_exp/2) % loopinfo->n1_exp - loopinfo->n1_exp/2;
            for (int i2 = 0; i2 < loopinfo->n2_exp; i2++ ){
                const int j2 = (i2 + loopinfo->n2_exp/2) % loopinfo->n2_exp - loopinfo->n2_exp/2;
                const int idx = 6*(i2+loopinfo->n2_exp*(i1+loopinfo->n1_exp*i0));
                N[idx+0] = thread_Nxxf(j0, j1, j2, loopinfo->dx, loopinfo->dy, loopinfo->dz);
                N[idx+1] = thread_Nxxg(j0, j1, j2, loopinfo->dx, loopinfo->dy, loopinfo->dz);
                N[idx+2] = thread_Nxxg(j0, j2, j1, loopinfo->dx, loopinfo->dz, loopinfo->dy);
                N[idx+3] = thread_Nxxf(j1, j2, j0, loopinfo->dy, loopinfo->dz, loopinfo->dx);
                N[idx+4] = thread_Nxxg(j1, j2, j0, loopinfo->dy, loopinfo->dz, loopinfo->dx);
                N[idx+5] = thread_Nxxf(j2, j0, j1, loopinfo->dz, loopinfo->dx, loopinfo->dy);
            }
        }
    }
    //std::cout << "finished all loops " << std::endl;
    //TODO//delete [] loopinfo;
    //loopinfo=NULL;
    
    //for (int i = core * MAX / 4; i < (core + 1) * MAX / 4; i++)  
    //    for (int j = 0; j < MAX; j++)  
    //        for (int k = 0; k < MAX; k++)  
    //            matC[i][j] += matA[i][k] * matB[k][j]; 
} 

af::array DemagFieldMultithread::N_cpp_alloc(int n0_exp, int n1_exp, int n2_exp, double dx, double dy, double dz){
    pthread_t threads[nthreads];
    //std::thread t[nthreads];
    //std::cout << "n0_exp= " <<  n0_exp << ", n_0_exp/nthreads= " << n0_exp/nthreads <<  std::endl;

    //LoopInfo loopinfo(i0_start, i0_end, n0_exp, n1_exp, n2_exp);
    //LoopInfo* loopinfo = malloc(sizeof(LoopInfo)); 
    //LoopInfo* loopinfo = (LoopInfo *)malloc(sizeof(LoopInfo *)); 
    
    //LoopInfo* loopinfo = NULL;

    //LoopInfo* loopinfo = new LoopInfo[1]; 
    //LoopInfo* loopinfo = new LoopInfo[nthreads]; 
    
    //LoopInfo* loopinfo = NULL;
    //loopinfo = new LoopInfo()[0]; 
    ////why not 1? //loopinfo = new LoopInfo[1]; 
    
    //LoopInfo* test_loopinfo = new LoopInfo(0, n0_exp, n0_exp, n1_exp, n2_exp, dx, dy, dz); 
    //LoopInfo* test_loopinfo = new LoopInfo()[10];
    struct LoopInfo loopinfo[nthreads];
    for (int i = 0; i < nthreads; i++){
        int start = i * n0_exp/nthreads;
        int end = (i +1) * n0_exp/nthreads;
        loopinfo[i]=LoopInfo(start, end, n0_exp, n1_exp, n2_exp, dx, dy, dz);
        //test_loopinfo[i].i0_start = 0;
        //test_loopinfo[i].i0_end = n0_exp;
        //test_loopinfo[i].n0_exp = n0_exp;
        //test_loopinfo[i].n1_exp = n1_exp;
        //test_loopinfo[i].n2_exp = n2_exp;
        //test_loopinfo[i].dx = dx;
        //test_loopinfo[i].dy = dz;
        //test_loopinfo[i].dz = dy;
    }
    //RUNNING//LoopInfo* loopinfo = new LoopInfo(0, n0_exp, n0_exp, n1_exp, n2_exp, dx, dy, dz); 

    
    //loopinfo->i0_start = 0;
    //loopinfo->i0_end   = n0_exp;
    //loopinfo->n0_exp   = n0_exp;
    //loopinfo->n1_exp   = n1_exp;
    //loopinfo->n2_exp   = n2_exp;
    //loopinfo->dx = dx;
    //loopinfo->dy = dz;
    //loopinfo->dz = dy;
    //std::cout << "setupt with loopinfo->n0_exp" << loopinfo->n0_exp << std::endl;

    N = new double[n0_exp*n1_exp*n2_exp*6];

    //printf("\nStarting treads:\n\n");
    for (int i = 0; i < nthreads; i++){
        //t[i] = std::thread(&this->setup_N, this, this->mesh);
        //t[i] = std::thread(this->setup_N,  this->mesh);
        //t[i] = std::thread(setup_N, i);
        int iret = pthread_create( &threads[i], NULL, setup_N, (void*) &loopinfo[i]);
        //int iret = pthread_create( &threads[i], NULL, setup_N, loopinfo);
        //TODO//int iret = pthread_create( &threads[i], NULL, setup_N, (void*) &loopinfo);
        if(iret)
        {
            fprintf(stderr,"Error - pthread_create() return code: %d\n",iret);
            exit(EXIT_FAILURE);
        }
     }

    //for (int i = 0; i < 1; i++){
    //printf("\n Waiting for threads:\n\n");
    for (int i = 0; i < nthreads; i++){
        //t[i].join();
    /* Waiting for all threads to be joined.*/
        int iret = pthread_join( threads[i], NULL);
        if(iret)
        {
            fprintf(stderr,"Error - pthread_join return code: %d\n",iret);
            exit(EXIT_FAILURE);
        }
     }
    //delete [] loopinfo;
    //loopinfo = NULL;

    //double N [n0_exp*n1_exp*n2_exp*6];
    af::array Naf(6,n2_exp,n1_exp,n0_exp,N);
    Naf=af::reorder(Naf,3,2,1,0);
    delete [] N;
    N = NULL;
    if (n2_exp == 1){
        Naf = af::fftR2C<2>(Naf);
    }
    else {
        Naf = af::fftR2C<3>(Naf);
    }
    return Naf;
}


double thread_newellf(double x, double y, double z){
    x=fabs(x);
    y=fabs(y);
    z=fabs(z);
    const double R = sqrt(pow(x,2) + pow(y,2) + pow(z,2));
    const double xx = pow(x,2);
    const double yy = pow(y,2);
    const double zz = pow(z,2);
  
    double result = 1.0 / 6.0 * (2.0*xx - yy - zz) * R;
    if(xx + zz > 0) result += y / 2.0 * (zz - xx) * asinh(y / (sqrt(xx + zz)));
    if(xx + yy > 0) result += z / 2.0 * (yy - xx) * asinh(z / (sqrt(xx + yy)));
    if(x  *  R > 0) result += - x*y*z * atan(y*z / (x * R));
    return result;
}

double thread_newellg(double x, double y, double z){
    z=fabs(z);
    const double R = sqrt(pow(x,2) + pow(y,2) + pow(z,2));
    const double xx = pow(x,2);
    const double yy = pow(y,2);
    const double zz = pow(z,2);
  
    double result = - x*y * R / 3.0;
    if(xx + yy > 0) result += x*y*z * asinh(z / (sqrt(xx + yy)));
    if(yy + zz > 0) result += y / 6.0 * (3.0 * zz - yy) * asinh(x / (sqrt(yy + zz)));
    if(xx + zz > 0) result += x / 6.0 * (3.0 * zz - xx) * asinh(y / (sqrt(xx + zz)));
    if(z  *  R > 0) result += - pow(z,3) / 6.0     * atan(x*y / (z * R));
    if(y  *  R!= 0) result += - z * yy / 2.0 * atan(x*z / (y * R));
    if(x  *  R!= 0) result += - z * xx / 2.0 * atan(y*z / (x * R));
    return result;
}

double thread_Nxxf(int ix, int iy, int iz, double dx, double dy, double dz){
    double x = dx*ix;
    double y = dy*iy;
    double z = dz*iz;
    double result = 8.0 * thread_newellf( x,    y,    z   ) \
           - 4.0 * thread_newellf( x+dx, y,    z   ) \
           - 4.0 * thread_newellf( x-dx, y,    z   ) \
           - 4.0 * thread_newellf( x,    y+dy, z   ) \
           - 4.0 * thread_newellf( x,    y-dy, z   ) \
           - 4.0 * thread_newellf( x,    y   , z+dz) \
           - 4.0 * thread_newellf( x,    y   , z-dz) \
           + 2.0 * thread_newellf( x+dx, y+dy, z   ) \
           + 2.0 * thread_newellf( x+dx, y-dy, z   ) \
           + 2.0 * thread_newellf( x-dx, y+dy, z   ) \
           + 2.0 * thread_newellf( x-dx, y-dy, z   ) \
           + 2.0 * thread_newellf( x+dx, y   , z+dz) \
           + 2.0 * thread_newellf( x+dx, y   , z-dz) \
           + 2.0 * thread_newellf( x-dx, y   , z+dz) \
           + 2.0 * thread_newellf( x-dx, y   , z-dz) \
           + 2.0 * thread_newellf( x   , y+dy, z+dz) \
           + 2.0 * thread_newellf( x   , y+dy, z-dz) \
           + 2.0 * thread_newellf( x   , y-dy, z+dz) \
           + 2.0 * thread_newellf( x   , y-dy, z-dz) \
           - 1.0 * thread_newellf( x+dx, y+dy, z+dz) \
           - 1.0 * thread_newellf( x+dx, y+dy, z-dz) \
           - 1.0 * thread_newellf( x+dx, y-dy, z+dz) \
           - 1.0 * thread_newellf( x+dx, y-dy, z-dz) \
           - 1.0 * thread_newellf( x-dx, y+dy, z+dz) \
           - 1.0 * thread_newellf( x-dx, y+dy, z-dz) \
           - 1.0 * thread_newellf( x-dx, y-dy, z+dz) \
           - 1.0 * thread_newellf( x-dx, y-dy, z-dz);
    return - result / (4.0 * M_PI * dx * dy * dz);
}

double thread_Nxxg(int ix, int iy, int iz, double dx, double dy, double dz){
    double x = dx*ix;
    double y = dy*iy;
    double z = dz*iz;
    double result = 8.0 * thread_newellg( x,    y,    z   ) \
                  - 4.0 * thread_newellg( x+dx, y,    z   ) \
                  - 4.0 * thread_newellg( x-dx, y,    z   ) \
                  - 4.0 * thread_newellg( x,    y+dy, z   ) \
                  - 4.0 * thread_newellg( x,    y-dy, z   ) \
                  - 4.0 * thread_newellg( x,    y   , z+dz) \
                  - 4.0 * thread_newellg( x,    y   , z-dz) \
                  + 2.0 * thread_newellg( x+dx, y+dy, z   ) \
                  + 2.0 * thread_newellg( x+dx, y-dy, z   ) \
                  + 2.0 * thread_newellg( x-dx, y+dy, z   ) \
                  + 2.0 * thread_newellg( x-dx, y-dy, z   ) \
                  + 2.0 * thread_newellg( x+dx, y   , z+dz) \
                  + 2.0 * thread_newellg( x+dx, y   , z-dz) \
                  + 2.0 * thread_newellg( x-dx, y   , z+dz) \
                  + 2.0 * thread_newellg( x-dx, y   , z-dz) \
                  + 2.0 * thread_newellg( x   , y+dy, z+dz) \
                  + 2.0 * thread_newellg( x   , y+dy, z-dz) \
                  + 2.0 * thread_newellg( x   , y-dy, z+dz) \
                  + 2.0 * thread_newellg( x   , y-dy, z-dz) \
                  - 1.0 * thread_newellg( x+dx, y+dy, z+dz) \
                  - 1.0 * thread_newellg( x+dx, y+dy, z-dz) \
                  - 1.0 * thread_newellg( x+dx, y-dy, z+dz) \
                  - 1.0 * thread_newellg( x+dx, y-dy, z-dz) \
                  - 1.0 * thread_newellg( x-dx, y+dy, z+dz) \
                  - 1.0 * thread_newellg( x-dx, y+dy, z-dz) \
                  - 1.0 * thread_newellg( x-dx, y-dy, z+dz) \
                  - 1.0 * thread_newellg( x-dx, y-dy, z-dz);
    result = - result / (4.0 * M_PI * dx * dy * dz);
    return result;
}
