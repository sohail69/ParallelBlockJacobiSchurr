/*==============================================================*\
  This is the pure Solid Mechanics code solving the active strain
  static solid mechanics problem using mixed U-P taylor hood
  elements, using 27-node Hexahedra coupled to 8-node hexahedra
  using the simple pressure boundary condition
\*==============================================================*/

#include <iomanip>
#include <iostream>
#include <ctime>
#include "InterfaceC.h"
#include "Interface_Classes.hpp"

using namespace std;
//Update an Increment
void UpdateIncrement(double* u, double* u0, double* u1 
                   , int nSize, int Lstep, int nLsteps){ //Lsteps start with 0
    for(int I=0; I<nSize; I++) 
      u[I]=u0[I] + (u1[I] - u0[I])*( double(Lstep+1)/double(nLsteps) );
};

//Traction Stress update (Simple)
void UpdateTStressFunc1(double* Stress, double t, int *NTOTStress
                      , int *DIM, int *nloaded){
  double press = ((t < 0.2)&&(t> 0.055)) ?(945.38*t*t - 245.54*t):(5.0E+03);
  press = press - (5.0E+03);
  if(*nloaded>0) quick_pressure_(Stress,&press,NTOTStress,DIM,nloaded);
};
void UpdateTStressFunc(double* Stress, double *press, int *NTOTStress
                      , int *DIM, int *nloaded){
  if(*nloaded>0) quick_pressure_(Stress,press,NTOTStress,DIM,nloaded);
};

int main(){
  int DIM = 3; //Cube geometry or any 3D uncomment
  int nst, nloaded = 0, npri, npri2, nskip, NUMPE;
  int nels_pp, nn_pp, neq_pp, NEQHeat, NEQSolid, ntots, nod, nodFace;
  int NTOTStress, NTOTHeat, NTOTSolid, NODOFSolid;
  npri  = 25; //ParaView Output 
  npri2 = 20; //Solids Solver and Statistics

  //Tissue cuboid sample testCase
  char input1[50]  = "mesh/CuboidMesh ";
  char input2[50]  = "mesh/CuboidMeshSolids ";
  char input3[50]  = "mesh/CuboidMeshTract ";
  char input4[50]  = "mesh/CuboidMeshHeat ";
  char input5[50]  = "mesh/CuboidMeshCISD ";

  char  output1[50] = "results/TestCase  ";
  char output11[50] = "results/TestCaseCol  ";
  char  output2[50] = "results/TestCaseDisp  ";
  char  output3[50] = "results/TestCaseVolt  ";
  char  output4[50] = "results/TestCaseCstrain  ";
  char  output5[50] = "results/TestCasePress  ";

  double *fibre, *diff, *Stress, *Stress0, *Stress1;
  double *astrain, *astrain0, *astrain1, *astemp_pp;
  double *Km_pp, *Mmat_pp, *Rm_pp, *u_pp, *du_pp;
  double *disp_pp, *press_pp;
  double xcentre[DIM] = {0.0,0.0,0.0};

  //Linear and Non-linear solver parameters
  int nTstep = 10001;
  double dt = 0.1, t = 0.0;

  //Linear and Non-linear solver parameters
  int nlsteps = 30, nllimit = 8, llimit = 70, liters = 0;
  double lerr, nlerror, nltol = 3.0E-08, ltol1 = 1.0E-21, ltol2 = 1.0E-08;

  //=================
  //Mesh and problem definition
  //=================
  Mesh                   mesh(input1, 17, DIM, npri);
  Heat_problem<Mesh>     heat(&mesh);
  Solid_problemUP<Mesh>  solids(&mesh,3);

  //General problem classes related to FEA 
  NUMPE      = (&mesh)->return_NUMPE();
  nodFace    = (&mesh)->nodFace;
  nod        = (&mesh)->nod;
  nels_pp    = (&mesh)->return_NELS_PP();
  nn_pp      = (&mesh)->nn_pp;
  NEQHeat    = (&heat)->return_NEQ_PP();
  NEQSolid   = (&solids)->return_NEQ_PP();
  NTOTHeat   = (&heat)->return_NTOTS();
  NODOFSolid = (&solids)->return_NODOF();
  NTOTSolid  = (&solids)->return_NTOTS();
  nst        = DIM*(DIM+1)/2;
  NTOTStress = nst*nodFace;

  //Heat Array Allocations
  astrain1 = new double[NEQHeat];
  astrain0 = new double[NEQHeat];
  astrain  = new double[NEQHeat];
  diff     = new double[DIM*nels_pp];
  fibre    = new double[DIM*DIM*nels_pp];
  astemp_pp = new double[NTOTHeat*nels_pp];


  //Solid Array Allocations
  Km_pp    = new double[NTOTSolid*NTOTSolid*nels_pp];
  Mmat_pp  = new double[NTOTSolid*NTOTSolid*nels_pp];
  Rm_pp    = new double[NEQSolid];
  u_pp     = new double[NEQSolid];
  du_pp    = new double[NEQSolid];
  disp_pp  = new double[DIM*nn_pp];
  press_pp = new double[nn_pp];
  for(int i = 0; i < NEQSolid; i++){
    Rm_pp[i] = u_pp[i] = du_pp[i] = 0.0;
    if(i < (DIM*nn_pp) ) disp_pp[i] = 0.0;
    if(i < nn_pp) press_pp[i] = 0.0;
  };

  //Heat problem set-up BCs, materials etc..
  Read_Fibre_Diffusion<Mesh>(input4,fibre,diff,DIM,&mesh);

  //Solids problem set-up BCs, materials etc..
  solids.Read_Set_Matprops(input2, &mesh);
  solids.Read_Set_DirchletBCs(input2, (mesh.nbnd), &mesh);
  
  //Have to manually input number of Traction nodes
  solids.Read_Set_TractionBCs(input3, 121, &mesh); 

  nloaded=0; nloaded = solids.return_LoadedFaces();
  if(nloaded > 0) Stress  = new double[NTOTStress*nloaded];
  if(nloaded > 0) Stress0 = new double[NTOTStress*nloaded];
  if(nloaded > 0) Stress1 = new double[NTOTStress*nloaded];
  if(nloaded > 0) for(int I=0; I<(NTOTStress*nloaded); I++) Stress[I] = Stress0[I] = Stress1[I]=0.0;

  //=================
  //Output Ensi Case files
  //=================
  mesh.ENSI_GEO_output(output1,16);
  mesh.ENSI_Partition(output1,16);
  mesh.ENSI_Traction((&solids)->gg_Face);

  double pressures1, pressures0;
  int torp = 0;
  nlsteps = 3;//12;
  nllimit = 1;//12;
  double Pressures[5]   = {5.000E+03, 12.092E+03, 15.960E+03, 14.932E+03, 5.000E+03};
  double Strains[5]     = {1.236E-03, 1.218E-03, 8.075E-02, 8.916E-02, 8.399E-02};
  int LoadStepCounts[5] = {1, 1, 120, 15, 15};
  for(int Lcases=0; Lcases<nlsteps; Lcases++){
    for(int I=0; I<NEQHeat; I++) astrain0[I]= ((Lcases != 0) ? Strains[Lcases-1] : 0.00);
    for(int I=0; I<NEQHeat; I++) astrain1[I]= Strains[Lcases];
    if(NUMPE==1) cout << " Non-Linear Newton Iterations :" << endl;
    pressures0 =  ((Lcases != 0) ?  Pressures[Lcases-1]*(1.00E-04) : 0.00);
    pressures1 = Pressures[Lcases]*(1.00E-04);
    UpdateTStressFunc(Stress0,&pressures0,&NTOTStress,&DIM,&nloaded);  //Simple Update
    UpdateTStressFunc(Stress1,&pressures1,&NTOTStress,&DIM,&nloaded);
    nlsteps = LoadStepCounts[Lcases];

    for(int Lsteps = 0; Lsteps<nlsteps; Lsteps++){//LoadSteps
      if(NUMPE==1) cout<<" Load-Step :"<<(Lsteps+1)<< " of "<<nlsteps<< endl;
      UpdateIncrement(Stress,Stress0,Stress1,(NTOTStress*nloaded),Lsteps,nlsteps);
      UpdateIncrement(astrain,astrain0,astrain1,NEQHeat,Lsteps,nlsteps);
      mesh.GATHER(astrain, astemp_pp, &heat);

      for(int NIters = 0; NIters<nllimit; NIters++){//NewtonSteps
        double one=1.0, mone=-1.0;
        double du_norm = 0.0, R_norm = 0.0;
        solids.Update_LinearSystem(Km_pp,Rm_pp,u_pp,astemp_pp,fibre,Stress,&mesh);
        R_norm = mesh.norm_p(Rm_pp, NEQSolid);

        if(mesh.IsItConverged(&R_norm, &nltol)) break;

        //Solve Linear system
        // 1-SSOR, 2-Chol, 3-LDL^T, 4-LU, 5-LDU, 6-ILU preconditioner
/*
        liters=0; lerr=0.0;
        mesh.Set_LsolverPORDER(153);
        mesh.Set_LsolverLIMIT(4);
        mesh.Set_LsolverSOLVER(2);
        mesh.Set_LsolverTOL(ltol2);
        mesh.Solve_LinearSystem(Km_pp,Mmat_pp,du_pp,Rm_pp,&liters,&lerr,&solids,4);
        du_norm = mesh.norm_p(du_pp, NEQSolid);

        //Output some sim stats
        if(NUMPE==1) cout << setw(14) << "Newton Iteration :"
                          << setw(8)  << NIters  << setw(8)  << liters
                          << setw(14) << du_norm << setw(14) << R_norm
                          << setw(14) << lerr    << endl;
*/
        mesh.Increment(u_pp,u_pp,du_pp,&one,&mone,&NEQSolid);
      }//NewtonSteps

      torp = Lcases+1;
      split_u_p_(disp_pp,press_pp,u_pp,&DIM,&nn_pp,&NEQSolid);
      mesh.ENSI_Data_output(output4,23,astrain, &torp, &heat,1);
      mesh.ENSI_Data_output(output2,20,disp_pp,nn_pp,DIM,&torp,1); 
      mesh.ENSI_Data_output(output5,21,press_pp,nn_pp,1,&torp,1);
    }//LoadSteps
  }//Load cases	

  //Clean up and finalisation
  if(nloaded > 0) delete[] Stress, Stress0, Stress1;
  delete[] astrain1, astrain0, astrain;
  delete[] diff, fibre, astemp_pp, disp_pp, press_pp;
  delete[] Mmat_pp, Km_pp, Rm_pp, u_pp, du_pp;
  solids.Finalise();
  heat.Finalise();
  mesh.Finalise();
  return 0;
}
