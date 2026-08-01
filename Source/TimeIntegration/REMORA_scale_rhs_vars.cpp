#include <REMORA.H>

using namespace amrex;

void
REMORA::scale_rhs_vars ()
{
    // These two routines are exact mathematical inverses that exist only so
    // AMReX can interpolate ru/rv/ru2d/rv2d conservatively between levels. In
    // floating point they are NOT inverses: cff is not a power of two, so a
    // value round-tripped through them comes back as fl(fl(x*cff)/cff), which
    // differs from x by an ulp in a large fraction of cells.
    //
    // ru/rv are recomputed from scratch every step so the round-trip is
    // invisible there, but ru2d/rv2d are the AB3 memory of the 3-D forcing:
    // what step n+1 reads is what step n wrote, perturbed. ROMS has no analogue,
    // and the perturbation lands on nothing but the barotropic solution, once
    // per baroclinic step, 26784 times a month. With a single level there is
    // nothing to interpolate, so skip it.
    if (finest_level == 0) return;
    for (int lev=0; lev<=finest_level;lev++) {
        MultiFab& mf_cons = *cons_new[lev];
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(mf_cons, TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            Array4<Real const> const& pm   = vec_pm[lev]->array(mfi);
            Array4<Real const> const& pn   = vec_pn[lev]->array(mfi);
            Array4<Real      > const& ru   = vec_ru[lev]->array(mfi);
            Array4<Real      > const& rv   = vec_rv[lev]->array(mfi);
            Array4<Real      > const& ru2d = vec_ru2d[lev]->array(mfi);
            Array4<Real      > const& rv2d = vec_rv2d[lev]->array(mfi);

            Box ubx = mfi.grownnodaltilebox(0,IntVect(NGROW,NGROW,0));
            Box vbx = mfi.grownnodaltilebox(1,IntVect(NGROW,NGROW,0));
            Box ubx2d = ubx; ubx2d.makeSlab(2,0);
            Box vbx2d = vbx; vbx2d.makeSlab(2,0);

            ParallelFor(ubx, 2, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                Real cff = (pm(i,j,0)+pm(i-1,j,0)) * (pn(i,j,0)+pn(i-1,j,0));
                ru(i,j,k,n) = ru(i,j,k,n) / cff;
            });

            ParallelFor(vbx, 2, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                Real cff = (pm(i,j,0)+pm(i,j-1,0)) * (pn(i,j,0)+pn(i,j-1,0));
                rv(i,j,k,n) = rv(i,j,k,n) / cff;
            });

            ParallelFor(ubx2d, 2, [=] AMREX_GPU_DEVICE (int i, int j, int , int n)
            {
                Real cff = (pm(i,j,0)+pm(i-1,j,0)) * (pn(i,j,0)+pn(i-1,j,0));
                ru2d(i,j,0,n) = ru2d(i,j,0,n) / cff;
            });

            ParallelFor(vbx2d, 2, [=] AMREX_GPU_DEVICE (int i, int j, int , int n)
            {
                Real cff = (pm(i,j,0)+pm(i,j-1,0)) * (pn(i,j,0)+pn(i,j-1,0));
                rv2d(i,j,0,n) = rv2d(i,j,0,n) / cff;
            });
        }
    }
}

void
REMORA::scale_rhs_vars_inv ()
{
    // These two routines are exact mathematical inverses that exist only so
    // AMReX can interpolate ru/rv/ru2d/rv2d conservatively between levels. In
    // floating point they are NOT inverses: cff is not a power of two, so a
    // value round-tripped through them comes back as fl(fl(x*cff)/cff), which
    // differs from x by an ulp in a large fraction of cells.
    //
    // ru/rv are recomputed from scratch every step so the round-trip is
    // invisible there, but ru2d/rv2d are the AB3 memory of the 3-D forcing:
    // what step n+1 reads is what step n wrote, perturbed. ROMS has no analogue,
    // and the perturbation lands on nothing but the barotropic solution, once
    // per baroclinic step, 26784 times a month. With a single level there is
    // nothing to interpolate, so skip it.
    if (finest_level == 0) return;
    for (int lev=0; lev<=finest_level;lev++) {
        MultiFab& mf_cons = *cons_new[lev];
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(mf_cons, TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            Array4<Real const> const& pm   = vec_pm[lev]->array(mfi);
            Array4<Real const> const& pn   = vec_pn[lev]->array(mfi);
            Array4<Real      > const& ru   = vec_ru[lev]->array(mfi);
            Array4<Real      > const& rv   = vec_rv[lev]->array(mfi);
            Array4<Real      > const& ru2d = vec_ru2d[lev]->array(mfi);
            Array4<Real      > const& rv2d = vec_rv2d[lev]->array(mfi);

            Box ubx = mfi.grownnodaltilebox(0,IntVect(NGROW,NGROW,0));
            Box vbx = mfi.grownnodaltilebox(1,IntVect(NGROW,NGROW,0));
            Box ubx2d = ubx; ubx2d.makeSlab(2,0);
            Box vbx2d = vbx; vbx2d.makeSlab(2,0);

            ParallelFor(ubx, 2, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                Real cff = (pm(i,j,0)+pm(i-1,j,0)) * (pn(i,j,0)+pn(i-1,j,0));
                ru(i,j,k,n) = ru(i,j,k,n) * cff;
            });

            ParallelFor(vbx, 2, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                Real cff = (pm(i,j,0)+pm(i,j-1,0)) * (pn(i,j,0)+pn(i,j-1,0));
                rv(i,j,k,n) = rv(i,j,k,n) * cff;
            });

            ParallelFor(ubx2d, 2, [=] AMREX_GPU_DEVICE (int i, int j, int , int n)
            {
                Real cff = (pm(i,j,0)+pm(i-1,j,0)) * (pn(i,j,0)+pn(i-1,j,0));
                ru2d(i,j,0,n) = ru2d(i,j,0,n) * cff;
            });

            ParallelFor(vbx2d, 2, [=] AMREX_GPU_DEVICE (int i, int j, int , int n)
            {
                Real cff = (pm(i,j,0)+pm(i,j-1,0)) * (pn(i,j,0)+pn(i,j-1,0));
                rv2d(i,j,0,n) = rv2d(i,j,0,n) * cff;
            });
        }
    }

}
