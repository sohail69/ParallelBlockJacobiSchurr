#include <iostream>
#include <iomanip>
#include "InterfaceC.h"
/*-----------------------------------------------*\
This file contains the C/C++ classes needed for
the finite element integration and mesh routines
needed for the electromechanical problem described
it interfaces with the PARAFEM libraries, the 
extension libraries in src and the Ref_element 
libraries effectively making these classes, wrapper
classes if you want to see implementation of the
linear-solver, integration schemes and mesh 
manipulations look at the src files, PARAFEM 
libraries and the Ref_element files. To implement
a new element:
1. Implement the element in a similar
   manner seen in Ref_elements
2. Make element integration algorithm in
   InterfaceF.f90 (same as others)
3. Write an Interface function in the InterfaceC.h
4. Make a concrete class based on abstract class
   problem and template MESH
5. test it in C++ code
Examples of this can be seen bellow

Author: Sohail Rathore
\*-----------------------------------------------*/


using namespace std;
/*-----------------------------------------------*\
!                                                 !
!             Abstract problem class              !
!                                                 !
\*-----------------------------------------------*/
class Problem{
  public:
    int *gg_pp, *gg_Face, *MASK;
    int nst=0, nloadedFace=0;

    virtual int return_NODOF()  = 0;
    virtual int return_NTOTS()  = 0;
    virtual int return_NEQ_PP() = 0;
    virtual void set_MASK(int* mask) = 0;
    virtual int return_LoadedFaces() = 0;
};

/*-----------------------------------------------*\
!                                                 !
!         Finite element Mesh comm class          !
!                                                 !
\*-----------------------------------------------*/
class Mesh{
  public:
    int numpe, npes, npes_pp;                     //MPI interfacing data
    int partitioner, meshgen;                     //Mesh partitioner and generation
    int iel_start, ieq_start;

    char element[15];                             //Element type
    int ncolour;                                  //Element colouring
    int nbnd, nn, nels, nn_pp, nels_pp, np_types; //problem size
    int ndim, nod, nip, nodFace, nipFace, nFace;  //Isoparametric element
    int *gnum_pp, *etype_pp, *gg_colour;          //Reference steering arrays
    double *coord_pp;                             //Element nodal coordinates

  private:
    int solver = 1, llimit = 3, ell = 220;      //Linear solver parameters
    double ltol = 1.0E-012;                     //Linear solver error tolerance

    //Mesh reading routines
    void Read_Job(char FName[50], int nlen, int npri){
      read_jobinitialise_(FName, &nlen, element, &partitioner, &npri
                        , &numpe, &meshgen, &npes, &nod, &nip, &nbnd
                        , &nFace, &nodFace, &nipFace,  &iel_start
                        , &nels_pp, &nels, &nn, &np_types);
    };


    void Read_Mesh(char FName[50]){
      if(etype_pp == NULL) etype_pp = new int[nels_pp];
      if(gnum_pp  == NULL) gnum_pp  = new int[nod*nels_pp];
      if(coord_pp == NULL) coord_pp = new double[nod*ndim*nels_pp];
      read_mesh_(FName, element, coord_pp, gnum_pp, etype_pp
               , &nod, &ndim, &nn, &nels_pp, &nn_pp, &iel_start
               , &ieq_start, &meshgen, &numpe, &npes,  &npes_pp);
    };


  public:
    //Mesh constructor-reader
    Mesh(char argv[50], int nlen, int DIM, int npri){
      gnum_pp   = NULL;
      etype_pp  = NULL;
      coord_pp  = NULL;
      gg_colour = NULL;
      ndim = DIM;
      Read_Job(argv, nlen, npri);
      Read_Mesh(argv);
      /* Preconditioner element colouring */
/*
      if(gg_colour==NULL) gg_colour = new int[nels_pp];
      preconditioner_colouring_(gg_colour,&ncolour,gnum_pp,&nod,&nn_pp,&nels_pp
                               ,&nn,&nels,&npes,&numpe);
*/
    };

    //Sizing routines
    int return_DIM(){return ndim;};
    int return_NOD(){return nod;};
    int return_NELS_PP(){return nels_pp;};
    int return_NN_PP(){return nn_pp;};
    int return_NUMPE(){return numpe;};

    //Read in BCs
    void Read_BCs(char FName[50], int* boundary_N, int nbnd1, int nodof){
      read_boundarynodes_(FName, boundary_N, &nbnd1, &nodof,  &numpe);
    };


    void Find_BoundaryNodes(int* gg_pp, int* boundary_N, int nbnd1, Problem *prob){
      int nodof  = prob->return_NODOF();
      int ntots  = prob->return_NTOTS();
      boundarysetnf_(gg_pp, gnum_pp, boundary_N, &ndim, &nodof
                   , &ntots, &nod, &nbnd1, &nels_pp, &npes, &numpe);
    };

    void Find_BoundaryNodes(int* gg_pp,int* boundary_N,int *nbnd1,int *ndof,int *nodof)
    {
      boundarysetnf_(gg_pp, gnum_pp, boundary_N, &ndim, nodof
                   , ndof, &nod, nbnd1, &nels_pp, &npes, &numpe);
    };


    void Find_BoundaryFaces(int* gg_Face, int *nloadedFace, int* boundary_N
                          , int nbnd1, int nodof){
      int boundaryID = 3;
      boundaryfacedetection_(gg_Face, nloadedFace, gnum_pp, boundary_N, &boundaryID
                           , element, &ndim, &nodof, &nod, &nFace, &nodFace, &nbnd1
                           , &nels_pp, &npes, &numpe);
    };

    //Read in Material Properties
    void Read_Materials(char FName[50], double* Matprops, Problem *prob){
      int nmat = 1;
      read_materialprops_(FName, &numpe, Matprops, &nmat, &np_types);
    };

    //Read in partitioned real data
    void Read_DataR(char FName[50], double* Data_pp, int dof, int n_pp){
      read_datar_(FName, &numpe, &npes, Data_pp, &n_pp, &dof);
    };


    //Read in partitioned integer data
    void Read_DataI(char FName[50], int* Data_pp, int dof, int n_pp){
      read_datai_(FName, &numpe, &npes, Data_pp, &n_pp, &dof);
    };


    //Linear Setting Solver routines
    void Set_LsolverTOL(double tol){ltol = tol;};
    void Set_LsolverLIMIT(int limit){llimit = limit;};
    void Set_LsolverPORDER(int PORDER){ell = PORDER;};
    void Set_LsolverSOLVER(int solvers){solver = solvers;};


    //Linear Solver (CG, BICGSTABL and GMRES) preconditioned
    void Solve_LinearSystem(double* Amat, double* Mmat, double* xvec, double* bvec
                           ,int *iters,double *error,Problem *prob, int precon){
      int ntots  = prob->return_NTOTS();
      int neq_pp = prob->return_NEQ_PP();
      int nodof  = prob->return_NODOF();
      int *MASK; MASK = new int[nodof*nod];
      prob->set_MASK(MASK);
      linearsolve_(Amat, Mmat, xvec, bvec, MASK, gg_colour, &ncolour
                  , &ntots, &nod, &nodof, &nels_pp, &neq_pp, &nn_pp, &ltol
                  , &llimit, iters, &ell, error, &solver, &precon);
      delete[] MASK;
    };



    //Real data gather operation
    void GATHER(double* x, double* pmul, Problem *prob){
      int ntots  = prob->return_NTOTS();
      int neq_pp = prob->return_NEQ_PP();
      int nodof  = prob->return_NODOF();
      int *MASK; MASK = new int[nodof*nod];
      prob->set_MASK(MASK);
      gather_m_(x, pmul, MASK, &ntots, &nodof
               , &nod, &nels_pp, &neq_pp, &nn_pp);
      delete[] MASK;
    };


    //Real data scatter operation
    void SCATTER(double* x, double* pmul, Problem *prob){
      int ntots  = prob->return_NTOTS();
      int neq_pp = prob->return_NEQ_PP();
      int nodof  = prob->return_NODOF();
      int *MASK; MASK = new int[nodof*nod];
      prob->set_MASK(MASK);
      scatter_m_(x, pmul, MASK, &ntots, &nodof
                , &nod, &nels_pp, &neq_pp, &nn_pp);
      delete[] MASK;
    };


    //Real EBE matrix vector multiplication
    void PARAMATVEC(double* storA, double* x, double* b, Problem *prob){
      int ntots  = prob->return_NTOTS();
      int neq_pp = prob->return_NEQ_PP();
      int nodof  = prob->return_NODOF();
      int *MASK; MASK = new int[nodof*nod];
      prob->set_MASK(MASK);
      para_matvec_(storA,x,b,MASK,&nels_pp,&nn_pp,&neq_pp,&ntots,&nod,&nodof);
      delete[] MASK;
    }


    //Increment a vector
    void Increment(double* unew, double* uold, double* du, double *a, double *b, int *n_pp){
      increment_(unew, uold, du, a, b, n_pp);
    };


    //Calculates a vector increment
    void Calc_increment(double* du, double a, double*  unew, double*  uold, int n_pp){
      double one = 1.00, Mone = -1.00;
      increment_(unew, uold, du, &one, &Mone, &n_pp);
      double aInv=1.0/a;
      scalar_vector_product_(du,  &aInv, du, &n_pp);
    };

    //Vector Norm
    double norm_p(double* r_pp, int nlength){
      double norm = 0.0;
      norm_pp_(&norm, r_pp, &nlength);
      return norm;
    };

    void SUM_PP(double *total, double *partial){
      sum_pp_(total,partial);
    };

    //Error stop condition
    bool IsItConverged(double *error, double *rtol){
      bool ConvergenceTest;
      int IsConverged=0;
      stop_cond_pp_(&IsConverged,error,rtol);
      ConvergenceTest = (IsConverged != 0);
      return ConvergenceTest;
    };


    //Output data routines
    void ENSI_GEO_output(char argv[50], int nlen){
      ensi_geo_output_(argv, element, &nlen, coord_pp, gnum_pp, &numpe
                      , &npes, &ndim, &nod, &nn, &nels_pp, &nn_pp);
    };


    void ENSI_Data_output(char argv[50], int nlen, double* xnew_pp
                        , int *j, Problem *prob, int var){
      int ntots  = prob->return_NTOTS();
      int neq_pp = prob->return_NEQ_PP();
      int nodof  = prob->return_NODOF();
      ensi_data_output_(argv, xnew_pp, &numpe, &npes, &nlen, &neq_pp
                      , &nn_pp, &nodof, &ndim, j, &var, &nod, element);
    };


    void ENSI_Data_output(char argv[50], int nlen, double* xnew_pp
                        , int nn_pp, int nodof, int *j, int var){
      int neq_pp = nodof*nn_pp;
      ensi_data_output_(argv, xnew_pp, &numpe, &npes, &nlen, &neq_pp
                      , &nn_pp, &nodof, &ndim, j, &var, &nod, element);
    };


    void ENSI_Partition(char argv[50], int nlen){
      ensi_partition_data_(argv,&numpe,&npes,&nlen,&nn_pp,&ndim);
    }

    void ENSI_Traction(int* gg_Face){
      ensi_traction_data_(element,gg_Face,&numpe,&npes,&ndim,&nFace,&nodFace
                         ,&nod,&nels_pp,&nn_pp);
    }

    //Finalise MPI process
    void Finalise(){
      if(etype_pp  != NULL) delete[] etype_pp;
      if(gnum_pp   != NULL) delete[] gnum_pp;
      if(gg_colour != NULL) delete[] gg_colour;
      if(coord_pp  != NULL) delete[] coord_pp;
      finalise_();
    };
};

/*-----------------------------------------------*\
!                                                 !
!          Broad Data reader functions            !
!                                                 !
\*-----------------------------------------------*/
template <typename MESH>
void Read_Fibre_Diffusion(char FName[50], double* Fibre
                        , double* Diff, int ndim, MESH *mesh){
  int dof = (ndim+1)*ndim;
  double *Data_pp; Data_pp = new double[dof*(mesh->nels_pp)];
  mesh->Read_DataR(FName, Data_pp, dof, mesh->nels_pp);
  fibre_diffusion_(Fibre, Diff, Data_pp, &dof, &ndim, &(mesh->nels_pp));
  delete[] Data_pp;
};

template <typename MESH>
void Read_Cell_Data(char FName[50], int* CellIDS, int* StimDefs, MESH *mesh){
  int dof = 2;
  int *Data_pp; Data_pp = new int[dof*(mesh->nn_pp)];
  int ndim = mesh->ndim;
  mesh->Read_DataI(FName, Data_pp, dof, (mesh->nn_pp));
  cisd_(CellIDS, StimDefs, Data_pp, &dof, &ndim, &(mesh->nn_pp));
  delete[] Data_pp;
};


/*-----------------------------------------------*\
!                                                 !
!          Heat equation problem class            !
!                                                 !
\*-----------------------------------------------*/
template <typename MESH>
class Heat_problem: public Problem{
  private:
    int nFace, nodof = 1, nels_pp;
    int nmask = 0, ntots = 0, neq_pp = 0, nr = 0;
    int *gg_pp, *gg_Face, *MASK;
    double *val_pp;
    double theta = 0.5;
  public:
    Heat_problem(MESH *mesh){
      int nod   = mesh->return_NOD();
      int nn_pp = mesh->return_NN_PP();
      ntots  = nodof*nod;
      neq_pp = nodof*nn_pp;
      nmask  = nod*nodof;
      nels_pp = mesh->return_NELS_PP();
      nFace   = mesh->nFace;
      gg_pp   = NULL;
      gg_Face = NULL;
      MASK    = NULL;
      val_pp  = NULL;
      MASK    = new int[nmask];
      calcheatmasks_(MASK, &nod, &nodof);
    };


    //Sizing routines
    virtual int return_NODOF(){return nodof;};
    virtual int return_NTOTS(){return ntots;};
    virtual int return_NEQ_PP(){return neq_pp;};
    virtual void set_MASK(int* mask){
      for(int i = 0; i < nmask; i++) mask[i] = MASK[i];
    };
    virtual int return_LoadedFaces(){return 0;};

    //Linear System routines;
    void Update_LinearSystem(double* KA_pp, double* KB_pp, double* utemp
                           , double* v_pp, double* b_pp, double* fibre
                           , double* diff, double *dtim, MESH *mesh){
      Update_Stiffness(KA_pp, KB_pp, utemp, fibre, diff, dtim, mesh);
      Update_RHS(KB_pp, v_pp, b_pp, mesh);
    };


    void Update_Stiffness(double* KA_pp, double* KB_pp, double* utemp
                        , double* fibre, double* diff, double *dtim
                        , MESH *mesh){
      //Heat-element integration
      heat_integration_(KA_pp, KB_pp, (mesh->element), (mesh->coord_pp)
                     , utemp, diff, fibre, &theta, dtim, &nels_pp, &(mesh->ndim)
                     , &(mesh->nod), &ntots, &(mesh->nip));
    };


    void Update_RHS(double* KB_pp, double* v_pp, double* b_pp, MESH *mesh){
      mesh->PARAMATVEC(KB_pp, v_pp, b_pp, this);
    };


    void Read_Set_BCS(char FName, int* boundary_N, Mesh *mesh){
      if(gg_pp==NULL) gg_pp = new int[ntots*nels_pp];
      if(gg_Face==NULL) gg_Face = new int[(nFace+1)*nels_pp];
    };


    //Finalise solver
    void Finalise(){
      if(gg_pp    != NULL) delete[] gg_pp;
      if(gg_Face  != NULL) delete[] gg_Face;
      if(MASK     != NULL) delete[] MASK;
      if(val_pp   != NULL) delete[] val_pp;
    };
};

/*-----------------------------------------------*\
!                                                 !
!      Solid mechanics Mixed-UP problem class     !
!    for mixed Hexahedra 27-Node-U and 8-Node-P   !
!                                                 !
\*-----------------------------------------------*/
template <typename MESH>
class Solid_problemUP: public Problem{
  private:
    int nTract, nDirch;
    int nFace, nodof, nodofU, nels_pp, nprop=2;
    int nmask=0, ntots=0, ndofU=0, neq_pp=0, nr=0;

    double *matprops, *val_pp;
  public:
    int *gg_pp, *gg_Face, *MASK;
    int nst=0, nloadedFace=0;

    Solid_problemUP(MESH *mesh, int mat){
      int ndim  = mesh->return_DIM();
      int nod   = mesh->return_NOD();
      int nn_pp = mesh->return_NN_PP();
      nels_pp   = mesh->return_NELS_PP();
      nodofU    = ndim;
      nodof     = ndim+1;
      nst       = ((ndim+1)*ndim)/2;
      ntots     = nodof*nod;
      ndofU     = nodofU*nod;
      neq_pp    = nodof*nn_pp;
      nmask     = nodof*nod;
      nFace     = mesh->nFace;
      gg_pp     = NULL;
      gg_Face   = NULL;
      MASK      = NULL;
      matprops  = NULL;
      val_pp    = NULL;
      MASK      = new int[nmask];
      calcsolidsmasks_(MASK, &nod, &nodof);
    };


    //Sizing routines
    virtual int return_NODOF(){return nodof;};
    virtual int return_NTOTS(){return ntots;};
    virtual int return_NEQ_PP(){return neq_pp;};
    virtual void set_MASK(int* mask){
      for(int i = 0; i < nmask; i++) mask[i] = MASK[i];
    };
    virtual int return_LoadedFaces(){return nloadedFace;};

    //Linear System routines;
    void Update_LinearSystem(double* Km_pp, double* r_pp, double* u_pp, double* astrain
                           , double* fibre, double* Stress, MESH *mesh){
      //Solid-mechanics element integration
      double *rtemp, *utemp;
      rtemp  = new double[ntots*nels_pp];
      utemp  = new double[ntots*nels_pp];
      mesh->GATHER(u_pp, utemp, this);
      solid_integration_up_(rtemp, Km_pp, utemp, astrain, fibre
                       , (mesh->coord_pp), gg_pp, gg_Face, val_pp, Stress
                       , matprops, &nels_pp, &ntots, &(mesh->ndim), &nst
                       , &(mesh->nip), &(mesh->nod), &nodof, &(mesh->nFace)
                       , &(mesh->nodFace), &(mesh->nipFace), &nloadedFace
                       , &nr, &nprop, (mesh->element));
      mesh->SCATTER(r_pp, rtemp, this);
      delete[] rtemp, utemp;
    };

    void Update_active_strain(double* F0Inv, double* astrain, double* fibre
                            , MESH *mesh){
      integration_active_strain_(F0Inv, astrain, fibre, (mesh->coord_pp)
                               , (mesh->element), &(mesh->ndim), &(mesh->nip)
                               , &(mesh->nod), &nels_pp);
    };

    void Read_Set_Matprops(char argv[50], Mesh *mesh){
      double *matprops1 = new double[nprop*(mesh->np_types)];
      mesh->Read_Materials(argv, matprops1, this);
      if(matprops != NULL)   delete[] matprops;
      if(matprops == NULL)   matprops = new double[nprop];
      for(int I=0; I<nprop; I++) matprops[I] = matprops1[I];
      delete[] matprops1;
    };

    void Read_Set_DirchletBCs(char argv[50], int nbnd, Mesh *mesh){
      nDirch = nbnd;
      if(gg_pp ==NULL) gg_pp  = new    int[ndofU*nels_pp];
      if(val_pp==NULL) val_pp = new double[ndofU*nels_pp];
      int *boundary_N; boundary_N = new int[nDirch*(nodofU+1)];
      mesh->Read_BCs(argv, boundary_N, nDirch, nodofU);
      mesh->Find_BoundaryNodes(gg_pp, boundary_N, &nbnd, &ndofU, &nodofU);
      for(int i = 0; i < ndofU*nels_pp; i++) val_pp[i] = 0.0;
      delete[] boundary_N;
    };

    void Read_Set_TractionBCs(char argv[50], int nbnd, Mesh *mesh){
      nTract = nbnd;
      if(gg_Face==NULL) gg_Face = new int[(nFace+1)*nels_pp];
      int *boundary_N; boundary_N = new int[nTract*(nodofU+1)];
      mesh->Read_BCs(argv, boundary_N, nTract, nodofU);
      mesh->Find_BoundaryFaces(gg_Face, &nloadedFace, boundary_N, nTract, nodofU);
      delete[] boundary_N;
    };

    //Finalise solver
    void Finalise(){
      if(gg_pp    != NULL) delete[] gg_pp;
      if(gg_Face  != NULL) delete[] gg_Face;
      if(matprops != NULL) delete[] matprops;
      if(MASK     != NULL) delete[] MASK;
      if(val_pp   != NULL) delete[] val_pp;
    };
};
