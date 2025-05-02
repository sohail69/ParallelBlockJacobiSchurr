#ifndef INTERFACEC_H
#define INTERFACEC_H 
//--------------------------------------------------------------------
// This is a file that interfaces the Fortran90/95 functions
// in the InterfaceF.f90 file to make them available as C functions
// effectively making this a function wrapper file
// when converting between the languages there are some rules:
//
// 1. A subroutine in FORTRAN90 is a void function in C/C++
// 2. FORTRAN90 variables are references/pointers in C/C++
// 3. All FORTRAN90 functions/subroutines are in lower case with 
//    an additional "_" at the end of the function name.
// 4. FORTRAN90 does array operations better so don't bother doing
//    them in C/C++
// 5. If FORTRAN90 is used in this way all global Allocations MUST be done
//    in C/C++
// 6. Local Allocations of arrays in FORTRAN90 must be deallocated properly
//    and not be globally accessible outside of a subroutine 
//    (GATHER-SCATTER based module exception (Just leave that one alone))
//
// Author: Sohail Rathore
//--------------------------------------------------------------------



extern "C"{
  //--------------------------------------------------------------------
  // Mesh and distributed mathematics related functions
  //--------------------------------------------------------------------
  //Job reader
  void read_jobinitialise_(char job_name[50], int *nlen, char element[15]
                         , int *partitioner, int *npri, int *numpes, int *mesh
                         , int *npess, int *nod, int *nip, int *nbnd, int *nFace
                         , int *nodFace, int *nipFace, int *iels_start, int *nel_pp
                         , int *nels, int *nn, int *np_types);


  //Mesh reader
  void read_mesh_(char job_name[50], char element[15], double* coord_pp
                , int* gnum_pp, int* etype_pp, int *nod, int *ndim, int *nn
                , int *nel_pp, int *neqs_pp, int *iels_start, int *ieqs_start
                , int *meshgen, int *numpes, int *npess, int *npes_pp);


  //Boundary node reader
  void read_boundarynodes_(char job_name[50], int* boundary_N, int *nbnd
                         , int *nodof, int *numpes);


  //Material property reader
  void read_materialprops_(char job_name[50], int *numpe, double* matprops
                         , int *nmat, int *np_types);


  //Real data reader
  void read_datar_(char JOB_NAME[50], int *numpes, int *npess, double* data_pp
                ,  int *n_pp, int *dof);


  //Integer data reader
  void read_datai_(char JOB_NAME[50], int *numpes, int *npess, int* data_pp
                ,  int *n_pp, int *dof);


  void boundarysetnf_(int* gg_pp, int* gnum_pp, int* boundary_N, int *ndim
                   , int *nodof, int *ntots, int *nod, int *nbnd, int *nel_pp
                   , int *npess, int *numpes);


  //detects boundary faces for flux boundary conditions
  void boundaryfacedetection_(int* gg_Face, int *nloadedFace, int* gnum_pp
                            , int* boundary_N, int *boundaryID, char element[15]
                            , int *ndim, int *nodof, int *nod, int *nFace
                            , int *nodFace, int *nbnd, int *nel_pp, int *npess
                            , int *numpes);


  //Multi-field Gather operation
  void gather_m_(double* x, double* pmul, int* MASK, int *ntots, int *nodof
               , int *nod, int *nel, int *neqs_pp, int *nn_pp);


  //Multi-field Scatter operation
  void scatter_m_(double* x, double* pmul, int* MASK, int *ntots, int *nodof
                , int *nod, int *nel, int *neqs_pp, int *nn_pp);


  //Multi-field Matrix vector product
  void para_matvec_(double* storA, double* x, double* b
                  , int* MASK, int *nel, int *nn_pp, int *neqs_pp
                  , int *ntots, int *nod, int *nodof);


  //Increment vector
  void increment_(double* unew, double* uold, double* du, double *a, double *b, int *n_pp);
  void dot_product_pp_(double *dproduct, double* u_pp, double* v_pp, int *n_pp);

  //Scalar divides a Vector
  void scalar_vector_product_(double* unew,double *a,double* uold,int *n_pp);

  //Calculate Vector 2-Norm
  void norm_pp_(double *twonorm, double* u_pp, int *neqs_pp);

  //Parallel scalar global sum
  void sum_pp_(double *Total, double *partial);
  
  //Parallel integer maxval
  void integer_maxval_pp_(int *globalMax, int *localMax);

  //Parallel real8/double maxval
  void real8_maxval_pp_(double *globalMax, double *localMax);

  //Parallel real8/double vector-maxval
  void real8_maxvalvec_pp_(double *globalMax, double* Vector, int *NEQs_pp);

  //Parallel real8/double vector-ABS-maxval
  void real8_absmaxvalvec_pp_(double *globalMax, double* Vector, int *NEQs_pp);

  //Output ENSI geometry file
  void ensi_geo_output_(char argv[50], char element[15], int *nlen
                      , double* gcoord_pp, int* gnum_pp, int *numpes
                      , int *npess, int *ndim, int *nod, int *nn
                      , int *nel_pp, int *nn_pp);

  //Traction load surface ENSI Output
  void ensi_traction_data_(char element[15], int* gg_Face, int *numpes
                         , int *npess, int *ndim, int *nFace, int *nodFace
                          , int *nod, int *nel_pp, int *nn_pp);

  //Output ENSI data file
  void ensi_data_output_(char argv[50], double* xnew_pp, int *numpes
                       , int *npess, int *nlen, int *neq_p, int *nn_pp
                       , int *nodof, int *ndim, int *j, int *var
                       , int *nod, char element[15]);


  void ensi_partition_data_(char argv[50], int *numpes, int *npess, int *nlen
                          , int *nn_pp, int *ndim);

  void ensi_elmcolour_output_(char argv[50], int* gg_colour, int *numpes
                             ,int *npess, int *nlen, int *nel_pp, int *ndim
                             ,int *nod, char element[15]);

  //Finalise parafem and MPI processes
  void finalise_();
  

  //--------------------------------------------------------------------
  // Linear solver related functions related to Mesh
  //--------------------------------------------------------------------
  // ELement colouring for preconditioners
  void preconditioner_colouring_(int* gg_colour, int *ncolour, int* gg_pp
                               , int *nod, int *nn_pp, int *nel_pp
                               , int *nn, int *nels, int *npes, int *numpes);


  //Linear solver (CG, BICGSTAB(L), GMRESR(L) )
  void linearsolve_(double* A_mat, double* M_mat, double* x_vec, double* b_vec
                  , int* NodalMask, int* gg_colour, int *ncolours
                  , int *ntots, int *nod, int *nodof, int *nel_pp, int *neqs_pp
                  , int *nn_pp, double *ltol, int *limit, int *iters, int *ell
                  , double *error, int *solver, int *precon);

  //Error stop condition
  void stop_cond_pp_(int *IsConverged, double *error, double *rtol);

  //--------------------------------------------------------------------
  // Static displacement-element Solid mechanics 
  //--------------------------------------------------------------------
  //Element nodal mask for solid element
  void calcsolidsmasks_(int* SolidsMask, int *nod, int *nodof);

  //Integrates the active strain components
  void integration_active_strain_(double* F0Inv, double* astrain, double* fibre
                                , double* coord, char element[15], int *ndim
                                , int *nip, int *nod, int *nel_pp);

  //Quick and dirty pressure
  void quick_pressure_(double* Stress, double *pressure, int *ntotsStress, int *ndim
                     , int *nloadedFace);

  //--------------------------------------------------------------------
  // Static displacement-pressure mixed-element Solid mechanics 
  //--------------------------------------------------------------------
  //Solid-mechanics element integration Jacobian-residual
  void solid_integration_up_(double* Residual, double* StoreKE, double* utemp
                        , double* astrain, double* fibre, double* coord
                        , int* gg_pp, int* gg_Face, double* val_pp, double* Stress
                        , double* MATPROP, int *nel_pp, int *ntots, int *ndim
                        , int *nst, int *nip, int *nod, int *nodof, int *nFace
                        , int *nodFace, int *nipFace, int *nloadedFace, int *nr
                        , int *np_types, int *nprop, int *material
                        , char element[15]);


  //--------------------------------------------------------------------
  // ALE Orthotropic Heat equation related functions
  //--------------------------------------------------------------------
  //Element nodal mask for heat element
  void calcheatmasks_(int* HeatMask, int *nod, int *nodof);


  //Heat equation element integration
  void heat_integration_(double* StorKA, double* StorKB, char element[15]
                       , double* coord, double* utemp, double* diff
                       , double* fibre, double *theta, double *dtim
                       , int *nel_pp, int *ndim, int *nod, int *ntots
                       , int *nip);

  //--------------------------------------------------------------------
  // Data separators
  //--------------------------------------------------------------------
  //Cell ID and stim defs
  void cisd_(int* CellID, int* StimDef, int* data_pp, int *dof
           , int *ndim, int *nn_pp);

  //Diffusion and Fibre orientation
  void fibre_diffusion_(double* Fibre, double* Diff, double* data_pp
                      , int *dof, int *ndim, int *nels_pp);

  void split_u_p_(double* u_pp, double* pressure, double* data_pp
                , int *nodof, int *nn_pp, int *neq_pp);


  //--------------------------------------------------------------------
  // PreCICE Adapter routines
  //--------------------------------------------------------------------
  void initialise_adapter_(const char* dataName, char readItCheckp[50]
                         , char writeInitialData[50], char writeItCheckp[50]
                         , const char* participantName, const char* config
                         , int *rank, int *commsize, int *dimensions);


  void set_datapoints_(double* gcoord_pp, int* gg_Face, int* vertexIDs
                     , char element[15], int *ndim, int *nFace, int *nodFace
                     , int *nel_pp, int *meshID, int *VertexSize, int *nloaded);


  void write_data_(double* utemp, int* gg_Face, int* vertexIDs
                 , char element[15], int *nodof, int *ndofU, int *nFace
                 , int *nodFace, int *nel_pp, int *displID, int *VertexSize
                 , int *nloaded, int *ndim);


  void read_data_(double* stress_pp, int* gg_Face, int* vertexIDs
                , char element[15], int *nodof, int *ndofS, int *nFace
                , int *nodFace, int *nel_pp, int *stressID, int *VertexSize
                , int *nloaded, int *ndim);

  //--------------------------------------------------------------------
  // Simulation Timing
  //--------------------------------------------------------------------
  void wtime_(double *time);
};
#endif
