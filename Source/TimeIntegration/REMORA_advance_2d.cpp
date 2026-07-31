#include <REMORA.H>
#include <REMORA_LoopBounds.H>

using namespace amrex;
/** Nonlinear shallow-water rpimitive equations predictor (Leap-frog) and
 * corrector (Adams-Moulton) time-stepping engine. Corresponds to Nonlinear/step2d_LF_AM3.h
 * in ROMS.
 *
 * @param[in   ] lev        level of refinement (coarsest level is 0)
 * @param[in   ] mf_rhoS    density perturbation
 * @param[in   ] mf_rhoA    vertically-averaged density
 * @param[inout] mf_ru2d    RHS contributions to 2D u-momentum
 * @param[inout] mf_rv2d    RHS contribtuions to 2D v-momentum
 * @param[inout] mf_rufrc   before first predictor, vertical integral of 3D RHS for uvel, converted to forcing terms
 * @param[inout] mf_rvfrc   before first predictor, vertical integral of 3D RHS for vvel, converted to forcing term
 * @param[inout] mf_Zt_avg1 average of sea surface height over all fast steps
 * @param[inout] mf_DU_avg1 time-averaged u-flux for 2D equations
 * @param[inout] mf_DU_avg2 time-averaged u-flux for 3D equation coupling
 * @param[inout] mf_DV_avg1 time-averaged v-flux for 2D equations
 * @param[inout] mf_DV_avg2 time-averaged v-flux for 3D equation coupling
 * @param[inout] mf_rubar   RHS of vertically integrated u-momentum
 * @param[inout] mf_rvbar   RHS of vertically integrated v-momentum
 * @param[inout] mf_rzeta   RHS of sea surface height
 * @param[inout] mf_ubar    vertically integrated u-momentum
 * @param[inout] mf_vbar    vertically integrated v-momentum
 * @param[inout] mf_zeta    Sea-surface height
 * @param[in   ] mf_h       Bathymetry
 * @param[in   ] mf_pm      1 / dx
 * @param[in   ] mf_pn      1 / dy
 * @param[in   ] mf_fcor    Coriolis factor
 * @param[inout] mf_visc2_p Harmonic viscosity at psi points
 * @param[inout] mf_visc2_r Harmoic viscosity at rho points
 * @param[in   ] mf_mskr    Land-sea mask at rho-points
 * @param[in   ] mf_msku    Land-sea mask at u-points
 * @param[in   ] mf_mskv    Land-sea mask at v-points
 * @param[in   ] mf_mskp    Land-sea mask at psi-points
 * @param[in   ] dtfast_lev Length of current barotropic step
 * @param[in   ] predictor_2d_step Is this a predictor step?
 * @param[in   ] first_2d_step Is this the first barotropic step?
 * @param[in   ] my_iif     Which barotropic predictor-corrector pair?
 * @param[inout] next_indx1 Cached index for
 */

void
REMORA::advance_2d (int lev,
                   MultiFab const* mf_rhoS,
                   MultiFab const* mf_rhoA,
                   MultiFab      * mf_ru2d,
                   MultiFab      * mf_rv2d,
                   MultiFab      * mf_rufrc,
                   MultiFab      * mf_rvfrc,
                   MultiFab      * mf_Zt_avg1,
                   std::unique_ptr<MultiFab>& mf_DU_avg1,
                   std::unique_ptr<MultiFab>& mf_DU_avg2,
                   std::unique_ptr<MultiFab>& mf_DV_avg1,
                   std::unique_ptr<MultiFab>& mf_DV_avg2,
                   std::unique_ptr<MultiFab>& mf_rubar,
                   std::unique_ptr<MultiFab>& mf_rvbar,
                   std::unique_ptr<MultiFab>& mf_rzeta,
                   std::unique_ptr<MultiFab>& mf_ubar,
                   std::unique_ptr<MultiFab>& mf_vbar,
                   MultiFab      * mf_zeta,
                   MultiFab const* mf_h,
                   MultiFab const* mf_pm,
                   MultiFab const* mf_pn,
                   MultiFab const* mf_fcor,
                   MultiFab const* mf_visc2_p,
                   MultiFab const* mf_visc2_r,
                   MultiFab const* mf_mskr,
                   MultiFab const* mf_msku,
                   MultiFab const* mf_mskv,
                   MultiFab const* mf_mskp,
                   Real dtfast_lev,
                   bool predictor_2d_step,
                   bool first_2d_step, int my_iif,
                   int & next_indx1)
{
    BL_PROFILE("REMORA::advance2d()");
    int iic = istep[lev];
    const int nnew  = 0;
    const int nstp  = 0;
    int ntfirst = 0;

    int knew = 3;
    // These closed forms replace ROMS's stateful indx1/next_indx1 toggling
    // (main3d.F). They agree with ROMS only when nfast is ODD, because that is
    // what makes indx1 flip parity exactly once per baroclinic step. With an
    // even nfast the two desynchronise after the first step and every
    // subsequent barotropic substep reads the wrong time level. NDTFAST=44
    // gives nfast=61, so Moana is safe. The parity is asserted once in
    // set_weights rather than here, which runs every barotropic substep.
    int krhs = (my_iif + iic) % 2 + 1;
    int kstp = my_iif <=1 ? iic % 2 + 1 : (iic % 2 + my_iif % 2 + 1) % 2 + 1;
    int indx1 = krhs;
    if (predictor_2d_step) {
        next_indx1 = 3 - indx1;
    } else {
        knew = next_indx1;
        kstp = 3 - knew;
        krhs = 3;
        //If it's not the auxiliary time step, set indx1 to next_indx1
        // NOTE: should this ever not execute?
        // Include indx1 updates for diagnostic purposes?
        //        if (my_iif<nfast+1)
        //            indx1=next_indx1;
    }
    int ptsk = 3-kstp;
    knew-=1;
    krhs-=1;
    kstp-=1;
    // Include indx1 updates for diagnostic purposes?
    //indx1-=1;
    ptsk-=1;
    auto ba = mf_h->boxArray();
    auto dm = mf_h->DistributionMap();

    MultiFab mf_DUon(convert(ba,IntVect(1,0,0)),dm,1,IntVect(NGROW,NGROW,0));
    MultiFab mf_DVom(convert(ba,IntVect(0,1,0)),dm,1,IntVect(NGROW,NGROW,0));

    int ncomp = 0;
    int fomn_comp = ncomp++;
    int Drhs_comp = ncomp++;
    int Dnew_comp = ncomp++;
    int zwrk_comp = ncomp++;
    int gzeta_comp = ncomp++;
    int gzeta2_comp = ncomp++;
    int gzetaSA_comp = ncomp++;
    int Dstp_comp = ncomp++;
    int rhs_ubar_comp = ncomp++;
    int rhs_vbar_comp = ncomp++;
    int rhs_zeta_comp = ncomp++;
    int zeta_new_comp = ncomp++;

    MultiFab mf(ba,dm,ncomp,IntVect(NGROW+1,NGROW+1,0));

    for ( MFIter mfi(*mf_rhoS, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real      > const& ubar = mf_ubar->array(mfi);
        Array4<Real      > const& vbar = mf_vbar->array(mfi);
        Array4<Real      > const& zeta = mf_zeta->array(mfi);
        Array4<Real const> const& h    = mf_h->const_array(mfi);

        Array4<Real const> const& pm   = mf_pm->const_array(mfi);
        Array4<Real const> const& pn   = mf_pn->const_array(mfi);

        Box bx = mfi.tilebox();
        Box gbx = mfi.growntilebox();
        Box gbx1 = mfi.growntilebox(IntVect(NGROW-1,NGROW-1,0));
        Box gbx2 = mfi.growntilebox(IntVect(NGROW,NGROW,0));
        Box xgbx2 = mfi.grownnodaltilebox(0, IntVect(NGROW,NGROW,0));
        Box ygbx2 = mfi.grownnodaltilebox(1, IntVect(NGROW,NGROW,0));

        Box tbxp1 = bx;
        Box tbxp2 = bx;
        Box tbxp3 = bx;
        tbxp1.grow(IntVect(NGROW-1,NGROW-1,0));
        tbxp2.grow(IntVect(NGROW,NGROW,0));
        tbxp3.grow(IntVect(NGROW+1,NGROW+1,0));

        Box bxD   = bx  ;   bxD.makeSlab(2,0);
        Box gbxD  = gbx ;  gbxD.makeSlab(2,0);
        Box gbx1D = gbx1; gbx1D.makeSlab(2,0);
        Box gbx2D = gbx2; gbx2D.makeSlab(2,0);

        Box tbxp2D = tbxp2;
        tbxp2D.makeSlab(2,0);

        // step2d work arrays
        FArrayBox fab_Drhs(makeSlab(tbxp3,2,0),1,The_Async_Arena());
        auto Drhs=fab_Drhs.array();

        auto DUon = mf_DUon.array(mfi);
        auto DVom = mf_DVom.array(mfi);

        ParallelFor(makeSlab(tbxp3,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            Drhs(i,j,0)=zeta(i,j,0,krhs)+h(i,j,0);
        });

        ParallelFor(makeSlab(xgbx2,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            Real on_u = two / (pn(i,j,0)+pn(i-1,j,0));
            Real cff1= Real(0.5) * on_u *(Drhs(i,j,0)+Drhs(i-1,j,0));
            DUon(i,j,0)=ubar(i,j,0,krhs)*cff1;
        });

        ParallelFor(makeSlab(ygbx2,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            Real om_v = two / (pm(i,j,0)+pm(i,j-1,0));
            Real cff1= Real(0.5) * om_v * (Drhs(i,j,0)+Drhs(i,j-1,0));
            DVom(i,j,0)=vbar(i,j,0,krhs)*cff1;
         });
    }

    // These are needed to pass the tests with bathymetry but I don't quite see why
    mf_DUon.FillBoundary(geom[lev].periodicity());
    mf_DVom.FillBoundary(geom[lev].periodicity());

#ifdef REMORA_USE_NETCDF
    if (solverChoice.do_m2_clim_nudg) {
        ubar_clim_data_from_file->update_interpolated_to_time(t_new[lev], lev, vec_ubar[lev].get(), geom, ref_ratio);
        vbar_clim_data_from_file->update_interpolated_to_time(t_new[lev], lev, vec_vbar[lev].get(), geom, ref_ratio);
    }
#endif

    for ( MFIter mfi(*mf_rhoS, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real const> const& rhoS = mf_rhoS->const_array(mfi);
        Array4<Real const> const& rhoA = mf_rhoA->const_array(mfi);
        Array4<Real const> const& h    = mf_h->const_array(mfi);

        Array4<Real      > const& rufrc   = mf_rufrc->array(mfi);
        Array4<Real      > const& rvfrc   = mf_rvfrc->array(mfi);
        Array4<Real      > const& Zt_avg1 = mf_Zt_avg1->array(mfi);
        Array4<Real      > const& ubar    = mf_ubar->array(mfi);
        Array4<Real      > const& vbar    = mf_vbar->array(mfi);
        Array4<Real      > const& zeta = mf_zeta->array(mfi);
        Array4<Real      > const& DU_avg1 = (mf_DU_avg1)->array(mfi);
        Array4<Real      > const& DU_avg2 = (mf_DU_avg2)->array(mfi);
        Array4<Real      > const& DV_avg1 = (mf_DV_avg1)->array(mfi);
        Array4<Real      > const& DV_avg2 = (mf_DV_avg2)->array(mfi);
        Array4<Real      > const& ru2d = (mf_ru2d)->array(mfi);
        Array4<Real      > const& rv2d = (mf_rv2d)->array(mfi);
        Array4<Real      > const& rubar = (mf_rubar)->array(mfi);
        Array4<Real      > const& rvbar = (mf_rvbar)->array(mfi);
        Array4<Real      > const& rzeta = (mf_rzeta)->array(mfi);
        Array4<Real const> const& visc2_p = mf_visc2_p->const_array(mfi);
        Array4<Real const> const& visc2_r = mf_visc2_r->const_array(mfi);

        Array4<Real const> const& pm   = mf_pm->const_array(mfi);
        Array4<Real const> const& pn   = mf_pn->const_array(mfi);
        Array4<Real const> const& fcor = mf_fcor->const_array(mfi);

        Array4<Real const> const& mskr = mf_mskr->const_array(mfi);
        Array4<Real const> const& msku = mf_msku->const_array(mfi);
        Array4<Real const> const& mskv = mf_mskv->const_array(mfi);
        Array4<Real const> const& mskp = mf_mskp->const_array(mfi);

        Box bx = mfi.tilebox();
        Box gbx = mfi.growntilebox();
        Box gbx1 = mfi.growntilebox(IntVect(NGROW-1,NGROW-1,0));
        Box gbx2 = mfi.growntilebox(IntVect(NGROW,NGROW,0));
        Box gbx3 = mfi.growntilebox(IntVect(NGROW+1,NGROW+1,0));
        Box xgbx2 = mfi.grownnodaltilebox(0, IntVect(NGROW,NGROW,0));
        Box ygbx2 = mfi.grownnodaltilebox(1, IntVect(NGROW,NGROW,0));

        Box xbxD = mfi.nodaltilebox(0);
        xbxD.makeSlab(2,0);

        Box ybxD = mfi.nodaltilebox(1);
        ybxD.makeSlab(2,0);

        // GPU-PARITY: ROMS's 2D momentum loops are `DO j=Jstr,Jend ;
        // DO i=IstrU,Iend` (u) and `DO j=JstrV,Jend ; DO i=Istr,Iend` (v), so
        // the domain-edge normal face is NOT stepped -- u2dbc/v2dbc set it. The
        // untrimmed nodaltilebox includes it. Every ParallelFor below that writes
        // rhs_ubar/rufrc/ru2d/ubar/rubar therefore takes the _adj box; see
        // Source/Utils/REMORA_LoopBounds.H for the term-by-term ROMS citations.
        Box xbxD_adj = roms_mom_box(xbxD, Geom(lev).Domain(), 0);
        Box ybxD_adj = roms_mom_box(ybxD, Geom(lev).Domain(), 1);

        Box tbxp1  = bx;  tbxp1.grow(IntVect(NGROW-1,NGROW-1,0));
        Box tbxp2  = bx;  tbxp2.grow(IntVect(NGROW,NGROW,0));
        Box tbxp3  = bx;  tbxp3.grow(IntVect(NGROW+1,NGROW+1,0));

        Box bxD   = bx;   bxD.makeSlab(2,0);
        Box gbxD  = gbx;  gbxD.makeSlab(2,0);
        Box gbx1D = gbx1; gbx1D.makeSlab(2,0);
        Box gbx2D = gbx2; gbx2D.makeSlab(2,0);

        Box tbxp2D = tbxp2;
        tbxp2D.makeSlab(2,0);

        auto fomn = mf.array(mfi,fomn_comp);
        auto Drhs = mf.array(mfi,Drhs_comp);
        auto Drhs_const = mf.const_array(mfi,Drhs_comp);
        auto Dnew = mf.array(mfi,Dnew_comp);
        auto zwrk = mf.array(mfi,zwrk_comp);
        auto gzeta = mf.array(mfi,gzeta_comp);
        auto gzeta2 = mf.array(mfi,gzeta2_comp);
        auto gzetaSA = mf.array(mfi,gzetaSA_comp);
        auto Dstp = mf.array(mfi,Dstp_comp);
        auto rhs_ubar = mf.array(mfi,rhs_ubar_comp);
        auto rhs_vbar = mf.array(mfi,rhs_vbar_comp);
        auto rhs_zeta = mf.array(mfi,rhs_zeta_comp);
        auto zeta_new = mf.array(mfi,zeta_new_comp);

        FArrayBox & fab_DUon=mf_DUon[mfi];
        FArrayBox & fab_DVom=mf_DVom[mfi];
        auto DUon=fab_DUon.array();
        auto DVom=fab_DVom.array();

        auto weight1 = vec_weight1.dataPtr();
        auto weight2 = vec_weight2.dataPtr();

        // Deliberately the UNTRIMMED box: this zeroing has no ROMS counterpart.
        // ROMS's rhs_ubar/rhs_vbar are stack scratch arrays that the
        // pressure-gradient loop ASSIGNS over IstrU:Iend; REMORA's live in `mf`,
        // a MultiFab freshly allocated (and therefore uninitialised) on every
        // call. Zeroing the two domain-edge faces as well keeps them
        // deterministic rather than garbage. Nothing reads them -- every
        // consumer below is on the _adj box -- so this cannot change the answer,
        // and it stops an uninitialised value reaching a plotfile.
        ParallelFor(xbxD, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            rhs_ubar(i,j,0)=zero;
        });

        ParallelFor(ybxD, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            rhs_vbar(i,j,0)=zero;
        });

        if (solverChoice.use_coriolis) {
            ParallelFor(tbxp2D, [=] AMREX_GPU_DEVICE (int i, int j, int  )
            {
                fomn(i,j,0) = fcor(i,j,0)*(one/(pm(i,j,0)*pn(i,j,0)));
            });
        }

        ParallelFor(makeSlab(tbxp3,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            Drhs(i,j,0)=zeta(i,j,0,krhs)+h(i,j,0);
        });

        if(predictor_2d_step)
        {
            if(first_2d_step) {
                Real cff2=(Real(-1.0)/Real(12.0))*weight2[my_iif+1];
                ParallelFor(makeSlab(gbx3,2,0),
                [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    Zt_avg1(i,j,0)=zero;
                });
                ParallelFor(makeSlab(xgbx2,2,0),
                [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    DU_avg1(i,j,0)=zero;
                    DU_avg2(i,j,0)=cff2*DUon(i,j,0);
                });
                ParallelFor(makeSlab(ygbx2,2,0),
                [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    DV_avg1(i,j,0)=zero;
                    DV_avg2(i,j,0)=cff2*DVom(i,j,0);
                });
            }
            else {
                Real cff1_wt1 = weight1[my_iif-1];
                Real cff2_wt1 = (Real(8.0)/Real(12.0))*weight2[my_iif]-
                                (one/Real(12.0))*weight2[my_iif+1];

                ParallelFor(makeSlab(gbx3,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    Zt_avg1(i,j,0) += cff1_wt1*zeta(i,j,0,krhs);
                });

                ParallelFor(makeSlab(xgbx2,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    DU_avg1(i,j,0) += cff1_wt1*DUon(i,j,0);
                    DU_avg2(i,j,0) += cff2_wt1*DUon(i,j,0);
                });

                ParallelFor(makeSlab(ygbx2,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    DV_avg1(i,j,0) += cff1_wt1*DVom(i,j,0);
                    DV_avg2(i,j,0) += cff2_wt1*DVom(i,j,0);
                });
            }
        }
        else {
            Real cff2_wt2;

            if (first_2d_step) {
                cff2_wt2=weight2[my_iif];
            } else {
                cff2_wt2=Real(5.0)/Real(12.0)*weight2[my_iif];
            }

            ParallelFor(makeSlab(xgbx2,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                DU_avg2(i,j,0)=DU_avg2(i,j,0)+cff2_wt2*DUon(i,j,0);
            });

            ParallelFor(makeSlab(ygbx2,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                DV_avg2(i,j,0)=DV_avg2(i,j,0)+cff2_wt2*DVom(i,j,0);
            });
        }
        //
        //  Do not perform the actual time stepping during the auxiliary
        //  (nfast(ng)+1) time step. Jump to next box
        //

        if (my_iif>=nfast) {
            continue; }
        //Load new free-surface values into shared array at both predictor
        //and corrector steps
        //
        //=======================================================================
        //  Time step free-surface equation.
        //=======================================================================
        //
        //  During the first time-step, the predictor step is Forward-Euler
        //  and the corrector step is Backward-Euler. Otherwise, the predictor
        //  step is Leap-frog and the corrector step is Adams-Moulton.
        //

        // todo: gzeta

        // ROMS step2d.F: fac = 1000/rho0. Hard-coding 1025 in the denominator
        // happens to be right for Moana and silently wrong for every other
        // rho0, with no error and no warning.
        Real fac=Real(1000.0)/solverChoice.rho0;

        if (my_iif==0) {
            Real cff1=dtfast_lev;

            ParallelFor(makeSlab(tbxp1,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                rhs_zeta(i,j,0) = (DUon(i,j,0)-DUon(i+1,j,0))+
                                  (DVom(i,j,0)-DVom(i,j+1,0));
                zeta_new(i,j,0) = (zeta(i,j,0,kstp)+ pm(i,j,0)*pn(i,j,0)*cff1*rhs_zeta(i,j,0)) * mskr(i,j,0);
                Dnew(i,j,0) = zeta_new(i,j,0)+h(i,j,0);

                //Pressure gradient terms:
                zwrk(i,j,0)=Real(0.5)*(zeta(i,j,0,kstp)+zeta_new(i,j,0));
                gzeta(i,j,0)=(fac+rhoS(i,j,0))*zwrk(i,j,0);
                gzeta2(i,j,0)=gzeta(i,j,0)*zwrk(i,j,0);
                gzetaSA(i,j,0)=zwrk(i,j,0)*(rhoS(i,j,0)-rhoA(i,j,0));
            });

        } else if (predictor_2d_step) {

            Real cff1=two * dtfast_lev;
            Real cff4=Real(4.0) / Real(25.0);
            Real cff5=one - two*cff4;

            ParallelFor(makeSlab(tbxp1,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                rhs_zeta(i,j,0)=(DUon(i,j,0)-DUon(i+1,j,0))+
                                (DVom(i,j,0)-DVom(i,j+1,0));
                zeta_new(i,j,0)=(zeta(i,j,0,kstp)+
                                pm(i,j,0)*pn(i,j,0)*cff1*rhs_zeta(i,j,0)) * mskr(i,j,0);
                Dnew(i,j,0)=zeta_new(i,j,0)+h(i,j,0);
                //Pressure gradient terms
                zwrk(i,j,0)=cff5*zeta(i,j,0,krhs)+
                    cff4*(zeta(i,j,0,kstp)+zeta_new(i,j,0));
                gzeta(i,j,0)=(fac+rhoS(i,j,0))*zwrk(i,j,0);
                gzeta2(i,j,0)=gzeta(i,j,0)*zwrk(i,j,0);
                gzetaSA(i,j,0)=zwrk(i,j,0)*(rhoS(i,j,0)-rhoA(i,j,0));
            });

        } else if (!predictor_2d_step) { //AKA if(corrector_2d_step)

            Real cff1=dtfast_lev * Real(5.0)/Real(12.0);
            Real cff2=dtfast_lev * Real(8.0)/Real(12.0);
            Real cff3=dtfast_lev * one/Real(12.0);
            Real cff4=two/Real(5.0);
            Real cff5=one-cff4;

            ParallelFor(makeSlab(tbxp1,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                Real cff=cff1*((DUon(i,j,0)-DUon(i+1,j,0))+
                               (DVom(i,j,0)-DVom(i,j+1,0)));
                zeta_new(i,j,0)=zeta(i,j,0,kstp)+
                    pm(i,j,0)*pn(i,j,0)*(cff+
                                         cff2*rzeta(i,j,0,kstp)-
                                         cff3*rzeta(i,j,0,ptsk));
                zeta_new(i,j,0) *= mskr(i,j,0);
                Dnew(i,j,0)=zeta_new(i,j,0)+h(i,j,0);
                //Pressure gradient terms
                zwrk(i,j,0)=cff5*zeta_new(i,j,0)+cff4*zeta(i,j,0,krhs);
                gzeta(i,j,0)=(fac+rhoS(i,j,0))*zwrk(i,j,0);
                gzeta2(i,j,0)=gzeta(i,j,0)*zwrk(i,j,0);
                gzetaSA(i,j,0)=zwrk(i,j,0)*(rhoS(i,j,0)-rhoA(i,j,0));
            });
        }

        //
        //  Load new free-surface values into shared array at both predictor
        //  and corrector steps.
        //
        //// zeta(knew) only valid at zeta_new, i.e. tbxp1
        ParallelFor(makeSlab(gbx1,2,0),
        [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            zeta(i,j,0,knew) = zeta_new(i,j,0);
        });

        //
        //  If predictor step, load right-side-term into shared array.
        //
        if (predictor_2d_step) {
            ParallelFor(makeSlab(gbx1,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                rzeta(i,j,0,krhs)=rhs_zeta(i,j,0);
            });
        }

        //
        //=======================================================================
        //  Compute right-hand-side for the 2D momentum equations.
        //=======================================================================
        //
/*
!
!-----------------------------------------------------------------------
!  Compute pressure gradient terms.
!-----------------------------------------------------------------------
!
*/
        Real cff1 = Real(0.5) * g;
        Real cff2 = one / Real(3.0);
        ParallelFor(xbxD_adj,
        [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            Real on_u = two / (pn(i,j,0)+pn(i-1,j,0));
            rhs_ubar(i,j,0)=cff1 * on_u *
                          ((    h(i-1,j,0) +     h(i,j,0))*
                           (gzeta(i-1,j,0) - gzeta(i,j,0))+
                           (    h(i-1,j,0) -     h(i,j,0))*
                           (      gzetaSA(i-1,j,0) + gzetaSA(i,j,0)+
                            cff2*(   rhoA(i-1,j,0) -    rhoA(i,j,0))*
                                 (   zwrk(i-1,j,0) -    zwrk(i,j,0)))+
                           (gzeta2(i-1,j,0)- gzeta2(i  ,j,0)));
        });

        ParallelFor(ybxD_adj,
        [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            Real om_v = two / (pm(i,j,0)+pm(i,j-1,0));
            rhs_vbar(i,j,0) = cff1*om_v *
                          ((    h(i,j-1,0) +     h(i,j,0))*
                           (gzeta(i,j-1,0) - gzeta(i,j,0))+
                           (    h(i,j-1,0) -     h(i,j,0))*
                           (gzetaSA(i,j-1,0)+ gzetaSA(i,j  ,0)+
                            cff2*(rhoA(i,j-1,0)- rhoA(i,j  ,0))*
                                 (zwrk(i,j-1,0)- zwrk(i,j  ,0)))+
                           (gzeta2(i,j-1,0)- gzeta2(i,j  ,0)));
        });

        // Advection terms for 2d ubar, vbar added to rhs_ubar and rhs_vbar
        //
        //-----------------------------------------------------------------------
        // rhs_uv_2d
        //-----------------------------------------------------------------------
        //
        Array4<Real const> const& ubar_const = mf_ubar->const_array(mfi);
        Array4<Real const> const& vbar_const = mf_vbar->const_array(mfi);

        rhs_uv_2d(lev,xbxD_adj, ybxD_adj, ubar_const, vbar_const, rhs_ubar, rhs_vbar, DUon, DVom, krhs);

        //-----------------------------------------------------------------------
        // Add Coriolis forcing
        //-----------------------------------------------------------------------
        if (solverChoice.use_coriolis) {
            // Coriolis terms for 2d ubar, vbar added to rhs_ubar and rhs_vbar
            //
            //-----------------------------------------------------------------------
            // coriolis
            //-----------------------------------------------------------------------
            //
            coriolis(xbxD_adj, ybxD_adj, ubar_const, vbar_const, rhs_ubar, rhs_vbar, Drhs, fomn, krhs, 0);
        }

        if (solverChoice.use_curvilinear_grid) {
            Array4<Real const> const& dndx = vec_dndx[lev]->const_array(mfi);
            Array4<Real const> const& dmde = vec_dmde[lev]->const_array(mfi);
            curvilinear(bxD, xbxD_adj, ybxD_adj, ubar_const, vbar_const, rhs_ubar, rhs_vbar, Drhs, dndx, dmde, krhs, 0);
        }

        //-----------------------------------------------------------------------
        //Add in horizontal harmonic viscosity.
        // Consider generalizing or copying uv3dmix, where Drhs is used instead of Hz and u=>ubar v=>vbar, drop dt terms
        //-----------------------------------------------------------------------
        uv3dmix(xbxD_adj, ybxD_adj, ubar, vbar, ubar, vbar, rhs_ubar, rhs_vbar,
                visc2_p, visc2_r, Drhs_const,
                pm, pn, mskp, krhs, nnew, zero);

#ifdef REMORA_USE_NETCDF
        if (solverChoice.do_m2_clim_nudg) {
            Array4<Real      > const& ubar_krhs    = mf_ubar->array(mfi, krhs);
            Array4<Real      > const& vbar_krhs    = mf_vbar->array(mfi, krhs);
            Array4<const Real> const& ubar_clim = ubar_clim_data_from_file->get_interpolated_mf(lev)->const_array(mfi);
            Array4<const Real> const& vbar_clim = vbar_clim_data_from_file->get_interpolated_mf(lev)->const_array(mfi);
            Array4<const Real> const& ubar_nudg_coeff = vec_nudg_coeff[BdyVars::ubar][lev]->const_array(mfi);
            Array4<const Real> const& vbar_nudg_coeff = vec_nudg_coeff[BdyVars::vbar][lev]->const_array(mfi);
            // Boxes are like this to match ROMS
            apply_clim_nudg(xbxD_adj, 1, 0, rhs_ubar, ubar_krhs, ubar_clim, ubar_nudg_coeff, Drhs_const, pm, pn);
            apply_clim_nudg(ybxD_adj, 0, 1, rhs_vbar, vbar_krhs, vbar_clim, vbar_nudg_coeff, Drhs_const, pm, pn);
        }
#endif

        //-----------------------------------------------------------------------
        // Coupling from 3d to 2d
        //-----------------------------------------------------------------------
        if (first_2d_step&&predictor_2d_step)
        {
            if (iic==ntfirst) {
                ParallelFor(xbxD_adj, [=] AMREX_GPU_DEVICE (int i, int j, int )
                {
                    rufrc(i,j,0)    -= rhs_ubar(i,j,0);
                    rhs_ubar(i,j,0) += rufrc(i,j,0);
                    ru2d(i,j,0,nstp)  = rufrc(i,j,0);
                });

                ParallelFor(ybxD_adj, [=] AMREX_GPU_DEVICE (int i, int j, int )
                {
                    rvfrc(i,j,0)    -= rhs_vbar(i,j,0);
                    rhs_vbar(i,j,0) += rvfrc(i,j,0);
                    rv2d(i,j,0,nstp)  = rvfrc(i,j,0);
                });

            } else if (iic==(ntfirst+1)) {

                ParallelFor(xbxD_adj, [=] AMREX_GPU_DEVICE (int i, int j, int )
                {
                    rufrc(i,j,0)=rufrc(i,j,0)-rhs_ubar(i,j,0);
                    rhs_ubar(i,j,0)=rhs_ubar(i,j,0)+Real(1.5)*rufrc(i,j,0)-Real(0.5)*ru2d(i,j,0,0);
                    ru2d(i,j,0,1)=rufrc(i,j,0);
                    Real r_swap= ru2d(i,j,0,1);
                    ru2d(i,j,0,1) = ru2d(i,j,0,0);
                    ru2d(i,j,0,0) = r_swap;
                });

                ParallelFor(ybxD_adj, [=] AMREX_GPU_DEVICE (int i, int j, int )
                {
                    rvfrc(i,j,0)=rvfrc(i,j,0)-rhs_vbar(i,j,0);
                    rhs_vbar(i,j,0)=rhs_vbar(i,j,0)+Real(1.5)*rvfrc(i,j,0)-Real(0.5)*rv2d(i,j,0,0);
                    rv2d(i,j,0,1)=rvfrc(i,j,0);
                    Real r_swap= rv2d(i,j,0,1);
                    rv2d(i,j,0,1) = rv2d(i,j,0,0);
                    rv2d(i,j,0,0) = r_swap;
                });

            } else {
                cff1=Real(23.0)/Real(12.0);
                cff2=Real(16.0)/Real(12.0);
                Real cff3= Real(5.0)/Real(12.0);

                ParallelFor(xbxD_adj, [=] AMREX_GPU_DEVICE (int i, int j, int )
                {
                    rufrc(i,j,0)=rufrc(i,j,0)-rhs_ubar(i,j,0);
                    rhs_ubar(i,j,0)=rhs_ubar(i,j,0)+
                        cff1*rufrc(i,j,0)-
                        cff2*ru2d(i,j,0,0)+
                        cff3*ru2d(i,j,0,1);
                    ru2d(i,j,0,1)=rufrc(i,j,0);
                    Real r_swap= ru2d(i,j,0,1);
                    ru2d(i,j,0,1) = ru2d(i,j,0,0);
                    ru2d(i,j,0,0) = r_swap;
                });

                ParallelFor(ybxD_adj, [=] AMREX_GPU_DEVICE (int i, int j, int )
                {
                    rvfrc(i,j,0)=rvfrc(i,j,0)-rhs_vbar(i,j,0);
                    rhs_vbar(i,j,0)=rhs_vbar(i,j,0)+
                          cff1*rvfrc(i,j,0)-
                          cff2*rv2d(i,j,0,0)+
                          cff3*rv2d(i,j,0,1);
                    rv2d(i,j,0,1)=rvfrc(i,j,0);

                    Real r_swap= rv2d(i,j,0,1);
                    rv2d(i,j,0,1) = rv2d(i,j,0,0);
                    rv2d(i,j,0,0) = r_swap;
                });
            }
        } else {
            ParallelFor(xbxD_adj, [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                rhs_ubar(i,j,0) += rufrc(i,j,0);
            });

            ParallelFor(ybxD_adj, [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                rhs_vbar(i,j,0) += rvfrc(i,j,0);
            });
        }

        //
        //=======================================================================
        //  Time step 2D momentum equations.
        //=======================================================================
        //
        //  Compute total water column depth.
        //
        ParallelFor(makeSlab(tbxp3,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
              Dstp(i,j,0)=zeta(i,j,0,kstp)+h(i,j,0);
        });

        //
        //  During the first time-step, the predictor step is Forward-Euler
        //  and the corrector step is Backward-Euler. Otherwise, the predictor
        //  step is Leap-frog and the corrector step is Adams-Moulton.
        //
        if (my_iif==0) {
            cff1=Real(0.5)*dtfast_lev;
            ParallelFor(xbxD_adj,
            [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                Real cff=(pm(i,j,0)+pm(i-1,j,0))*(pn(i,j,0)+pn(i-1,j,0));
                Real Dnew_avg =one/(Dnew(i,j,0)+Dnew(i-1,j,0));
                ubar(i,j,0,knew)=(ubar(i,j,0,kstp)*
                                 (Dstp(i,j,0)+Dstp(i-1,j,0))+
                                  cff*cff1*rhs_ubar(i,j,0))*Dnew_avg * msku(i,j,0);
            });
            ParallelFor(ybxD_adj,
            [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                Real cff=(pm(i,j,0)+pm(i,j-1,0))*(pn(i,j,0)+pn(i,j-1,0));
                Real Dnew_avg=one/(Dnew(i,j,0)+Dnew(i,j-1,0));
                vbar(i,j,0,knew)=(vbar(i,j,0,kstp)*
                                 (Dstp(i,j,0)+Dstp(i,j-1,0))+
                                  cff*cff1*rhs_vbar(i,j,0))*Dnew_avg * mskv(i,j,0);
            });

        } else if (predictor_2d_step) {

            cff1=dtfast_lev;
            ParallelFor(xbxD_adj,
            [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                Real cff=(pm(i,j,0)+pm(i-1,j,0))*(pn(i,j,0)+pn(i-1,j,0));
                Real Dnew_avg=one/(Dnew(i,j,0)+Dnew(i-1,j,0));
                ubar(i,j,0,knew)=(ubar(i,j,0,kstp)*
                                 (Dstp(i,j,0)+Dstp(i-1,j,0))+
                                  cff*cff1*rhs_ubar(i,j,0))*Dnew_avg * msku(i,j,0);
            });
            ParallelFor(ybxD_adj,
            [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                Real cff=(pm(i,j,0)+pm(i,j-1,0))*(pn(i,j,0)+pn(i,j-1,0));
                Real Dnew_avg=one/(Dnew(i,j,0)+Dnew(i,j-1,0));
                vbar(i,j,0,knew)=(vbar(i,j,0,kstp)*
                                 (Dstp(i,j,0)+Dstp(i,j-1,0))+
                                  cff*cff1*rhs_vbar(i,j,0))*Dnew_avg * mskv(i,j,0);
            });

        } else if ((!predictor_2d_step)) {

            cff1=Real(0.5)*dtfast_lev*Real(5.0)/Real(12.0);
            cff2=Real(0.5)*dtfast_lev*Real(8.0)/Real(12.0);
            Real cff3=Real(0.5)*dtfast_lev*one/Real(12.0);
            ParallelFor(xbxD_adj,
            [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                Real cff=(pm(i,j,0)+pm(i-1,j,0))*(pn(i,j,0)+pn(i-1,j,0));
                Real Dnew_avg=one/(Dnew(i,j,0)+Dnew(i-1,j,0));
                ubar(i,j,0,knew)=(ubar(i,j,0,kstp)*
                                 (Dstp(i,j,0)+Dstp(i-1,j,0))+
                                 cff*(cff1*rhs_ubar(i,j,0)+
                                      cff2*rubar(i,j,0,kstp)-
                                      cff3*rubar(i,j,0,ptsk)))*Dnew_avg * msku(i,j,0);
            });
            ParallelFor(ybxD_adj,
            [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                Real cff=(pm(i,j,0)+pm(i,j-1,0))*(pn(i,j,0)+pn(i,j-1,0));
                Real Dnew_avg=one/(Dnew(i,j,0)+Dnew(i,j-1,0));
                vbar(i,j,0,knew)=(vbar(i,j,0,kstp)*
                                 (Dstp(i,j,0)+Dstp(i,j-1,0))+
                                 cff*(cff1*rhs_vbar(i,j,0)+
                                      cff2*rvbar(i,j,0,kstp)-
                                      cff3*rvbar(i,j,0,ptsk)))*Dnew_avg * mskv(i,j,0);
            });
        }

        //store rhs_ubar and rhs_vbar to save later
        //
        //  If predictor step, load right-side-term into shared arrays for
        //  future use during the subsequent corrector step.
        //

        if (predictor_2d_step) {
            ParallelFor(xbxD_adj, [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                rubar(i,j,0,krhs)=rhs_ubar(i,j,0);
            });
            ParallelFor(ybxD_adj, [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                rvbar(i,j,0,krhs)=rhs_vbar(i,j,0);
            });
        }
    }

    // GPU-PARITY: ROMS applies NO 2D boundary condition and NO river momentum
    // source on the auxiliary (nfast+1) barotropic step. An earlier comment here
    // asserted the opposite -- "step2d calls zetabc/u2dbc/v2dbc every time" --
    // and the guard was widened to `my_iif<=nfast` on the strength of it. The
    // preprocessed ROMS source that was actually compiled says otherwise
    // (atlas external/moana_spec/, from step2d.F -> step2d_LF_AM3.h):
    //
    //   step2d.f90:887-889   "Do not perform the actual time stepping during
    //                         the auxiliary (nfast(ng)+1) time step."
    //                         IF (iif(ng).gt.nfast(ng)) RETURN
    //
    // and every routine this block is the analogue of sits AFTER that return:
    //
    //   step2d.f90:1005  LwSrc  mass point source
    //   step2d.f90:1020  CALL zetabc_tile
    //   step2d.f90:1665  CALL u2dbc_tile
    //   step2d.f90:1670  CALL v2dbc_tile
    //   step2d.f90:1695  LuvSrc river momentum point source
    //
    // It is deliberate, not incidental. main3d.f90:736-738 states what the
    // auxiliary step is for: "No actual time-stepping is performed during the
    // auxiliary (nfast+1) time-step. It is needed to finalize the fast-time
    // averaging of 2D fields". Everything before the return -- Drhs/DUon/DVom
    // and the Zt_avg1/DU_avg1/DV_avg1/DU_avg2/DV_avg2 accumulation
    // (step2d.f90:750-863) -- does run; nothing is time-stepped, so there is no
    // new state to apply a boundary condition to. The only boundary work ROMS
    // does on that step is the periodic exchange of the *averaged* fields
    // (step2d.f90:869-885), which is a no-op in this non-periodic domain.
    //
    // REMORA's loop is 0-based, `my_iif = 0..nfast` (REMORA_Advance.cpp:61,71
    // and REMORA_TimeStepML.cpp:107,109), i.e. `my_iif == iif-1`, so
    // `my_iif == nfast` IS ROMS's auxiliary step. The port already mirrors the
    // early return at line ~395 above (`if (my_iif>=nfast) continue;`) -- but
    // that `continue` only escapes the MFIter loop, and this block sits outside
    // it, so it was not covered. `my_iif<nfast` is the guard that matches.
    //
    // The corrector is skipped on the same step (advance_2d_onestep,
    // `my_iif < nfast_counter-1`), which matches main3d.f90:766 guarding the
    // corrector call with `IF (iif(ng).lt.(nfast(ng)+1))`.
    if (my_iif<nfast) {
        int know;
        Real dt2d;
        if (my_iif==0) {
            know = krhs;
            dt2d = dtfast_lev;
        } else if (predictor_2d_step) {
            know = krhs;
            dt2d = two * dtfast_lev;
        } else {
            know = kstp;
            dt2d = dtfast_lev;
        }

        MultiFab ubar_know(*vec_ubar[lev], make_alias, know, 1);
        MultiFab vbar_know(*vec_vbar[lev], make_alias, know, 1);
        MultiFab zeta_know(*vec_zeta[lev], make_alias, know, 1);
        // t_old, not t_new: ROMS evaluates open-boundary data at the time the
        // step STARTS. main3d.F:522 increments time(ng), then set_data
        // interpolates, then the history record is written, then the state is
        // stepped -- and since history record 0 carries ocean_time = T0, the
        // loop must be entered with time(ng) = T0 - dt. So the step T0 -> T0+dt
        // consumes boundary data at T0. Measured: a scan of
        // remora.bdy_time_shift has a sharp minimum at exactly -dt (u max
        // 1.91e-06 at -100 s, against 2.89e-05 at -90, 1.70e-05 at -110,
        // 2.89e-04 unshifted and 5.78e-04 at +100).
        FillPatch(lev, t_old[lev], *vec_ubar[lev], GetVecOfPtrs(vec_ubar), ubar_bc(), BdyVars::ubar,
                  knew, false,true, 0,know, dt2d, ubar_know);
        FillPatch(lev, t_old[lev], *vec_vbar[lev], GetVecOfPtrs(vec_vbar), vbar_bc(), BdyVars::vbar,
                  knew, false,true, 0,know, dt2d, vbar_know);
        FillPatch(lev, t_old[lev], *vec_zeta[lev], GetVecOfPtrs(vec_zeta), zeta_bc(), BdyVars::zeta,
                  knew, false,false, 0,know, dt2d, zeta_know);

#ifdef REMORA_USE_NETCDF
        if (solverChoice.do_rivers) {
            // t_old here but t_new + river_time_shift for river_source_transport
            // (advance_3d.cpp) and river_source_cons (prestep.cpp). That
            // inconsistency is pre-existing and is left visible rather than
            // silently harmonised: a scan of remora.river_time_shift decides which
            // is right, the same way it decided set_tides.
            river_source_transportbar->update_interpolated_to_time(t_old[lev]);
            int* river_direction_d = river_direction.data();
            for ( MFIter mfi(*mf_rhoS, TilingIfNotGPU()); mfi.isValid(); ++mfi )
            {
                Array4<const int > const& river_pos = vec_river_position[lev]->const_array(mfi);
                Array4<const Real> const& river_transportbar = river_source_transportbar->fab_interp->array();
                Array4<Real      > const& ubar    = mf_ubar->array(mfi);
                Array4<Real      > const& vbar    = mf_vbar->array(mfi);
                Array4<Real const> const& zeta    = mf_zeta->const_array(mfi);
                Array4<Real const> const& h       = mf_h->const_array(mfi);
                Array4<Real const> const& pm      = mf_pm->const_array(mfi);
                Array4<Real const> const& pn      = mf_pn->const_array(mfi);

                Box gbx1D = mfi.growntilebox(IntVect(NGROW-1,NGROW-1,0));
                gbx1D.makeSlab(2,0);

                ParallelFor(gbx1D, [=] AMREX_GPU_DEVICE (int i, int j, int )
                {
                    int iriver = river_pos(i,j,0);
                    if (iriver >= 0) {
                        if (river_direction_d[iriver] == 0) {
                            Real on_u = two / (pn(i,j,0)+pn(i-1,j,0));
                            Real cff = one / (on_u * Real(0.5) * (zeta(i-1,j,0,knew) + h(i-1,j,0) +
                                        zeta(i,j,0,knew) + h(i,j,0)));
                            ubar(i,j,0,knew) = river_transportbar(iriver,0,0) * cff;
                        } else {
                            Real om_v = two / (pm(i,j,0)+pm(i,j-1,0));
                            Real cff = one / (om_v * Real(0.5) * (zeta(i,j-1,0,knew) + h(i,j-1,0) +
                                        zeta(i,j,0,knew) + h(i,j,0)));
                            vbar(i,j,0,knew) = river_transportbar(iriver,0,0) * cff;
                        }
                    }
                });
            }
        }
        FillPatchNoBC(lev, t_old[lev], *vec_ubar[lev], GetVecOfPtrs(vec_ubar), BdyVars::ubar,
                  knew, false,false);
        FillPatchNoBC(lev, t_old[lev], *vec_vbar[lev], GetVecOfPtrs(vec_vbar), BdyVars::vbar,
                  knew, false,false);
#endif
    }
}
