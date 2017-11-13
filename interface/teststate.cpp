#include "teststate.hpp"
#include<cstdio>
testState::testState (Mesh mesh_in, Param param_in, long int aptr):
              mesh(mesh_in),param(param_in)
{
  printf("%p a CPP array address \n", (void*)aptr);
  printf("%p a CPP array address \n", (void**)aptr);
  //af::info();
  //af_array a = *(void **)aptr;
  void **a = (void **)aptr;
  af::array *A = new af::array(*a);
  printf("%p b.get() CPP array address \n", A->get());
  std::cout<<"segtest0"<<std::endl;
  af::print("ARRAY:", *A);
  //af_print(*A);
  std::cout<<"segtest1"<<std::endl;
  m=*A;//TODO check!!!
  af::print("ARRAY m:", m);
  //std::cout<<"segtest3"<<std::endl;
}

void testState::printn0(){
  std::cout<<"mesh.n0="<<mesh.n0  <<std::endl;
}
//#include "teststate.hpp"
//#include<cstdio>
//testState::testState (Mesh mesh_in, Param param_in, long int aptr):
//              mesh(mesh_in),param(param_in)
//{
//  printf("%p a CPP array address \n", (void*)aptr);
//  //af::info();
//  //af_array a = *(void **)aptr;
//  void **a = (void **)aptr;
//  af::array *A = new af::array(*a);
//  printf("%p b.get() CPP array address \n", A->get());
//  std::cout<<"segtest0"<<std::endl;
//  //af::print("ARRAY:", *A);
//  std::cout<<"segtest1"<<std::endl;
//  //m=*A;//TODO check!!!
//  //af::print("ARRAY m:", m);
//  std::cout<<"segtest3"<<std::endl;
//}
//
//void testState::printn0(){
//  std::cout<<"mesh.n0="<<mesh.n0  <<std::endl;
//}
