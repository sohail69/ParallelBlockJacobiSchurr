
struct UP_MixedField{
  //Element variables
  u_Esol(ndim*nodU,nel_pp)
  p_Esol(nodP,nel_pp)

  //Global entity variables
  u_Gsol(ndim*nn_ppU)
  p_Gsol(nodP,nel_pp)
};

struct ActiveStrainNeoHookean{
//Material Inputs
g_val(nodU,nel_pp)
fibre(ndim,ndim,nodU,nel_pp)



};

struct TaylorHoodHexBasis{
N_im(nodU,nip)
dN_ikm(ndim,nodU,nip)
H_n(nodP,nip)
dH_kn(ndim,nodP,nip)
Wi(nip)
};


struct HierachicalHexMesh{
//Vertex coords
xCoord(ndim,nVerts,nel_pp)

//Forward Entity maps
ggEVerts(nEVerts, nel_pp);
ggEEdges(nEEdges, nel_pp);
ggEFaces(nEFaces, nel_pp);
ggEVolms(nEVolms, nel_pp);

//Inverse Entity Maps
ggGVerts(nelMax,nGVerts);
ggGEdges(nelMax,nGEdges);
ggGFaces(nelMax,nGFaces);
ggGVolms(nelMax,nGVolms);

};


template<typename Vars, typename MatProps, typename FEBasis, typename Mesh>
void Integrate_solid_elms(Vars Var, Matprops Mat, FEBasis FE, Mesh mesh, ){


}






#pragma omp target private(Iel_start, Iel_end, I)
template<unsigned IELstart, unsigned IELend>
void Integrate_solid_elms(){

  for(int Igauss=0; Igauss<nip; Igauss++){//Integrate over Gauss-Points
    



  }

}


//Outputs
elm_res(ntot,nel_pp)
elm_jac(ntot,ntot,nel_pp)
