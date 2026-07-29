#include "REMORA.H"

using namespace amrex;

#ifdef REMORA_USE_NETCDF
/*
 * @param[in   ] lev           level to operate on
 * @param[inout] mf_to_fill    data on which to apply BCs
 * @param[in   ] mf_mask       land-sea mask
 * @param[in   ] time          current time
 * @param[in   ] bccomp        index into both domain_bcs_type_bcr and bc_extdir_vals for icomp=0
 * @param[in   ] bdy_var_type  which netcdf boundary data to fill from
 * @param[in   ] icomp_to_fill component to update
 * @param[in   ] icomp_calc    component to reference from on RHS
 * @param[in   ] mf_calc       data for RHS of calculation
 * @param[in   ] dt_calc       time step for the calculation
 */

void
REMORA::fill_from_bdyfiles (int lev, MultiFab& mf_to_fill, const MultiFab& mf_mask, const Real time, const int bccomp,
                            const int bdy_var_type, const int icomp_to_fill, const int icomp_calc, const MultiFab& mf_calc, const Real dt_calc)
{
    // Which variable are we filling
    int ivar = bdy_var_type;

    //
    // Note that "domain" is mapped onto the type of box the data is in
    //
    Box domain = geom[lev].Domain();

    const auto& mf_index_type = mf_to_fill.boxArray().ixType();
    domain.convert(mf_index_type);

    const auto& dom_lo = amrex::lbound(domain);
    const auto& dom_hi = amrex::ubound(domain);

    int ncomp;

    // If we are doing the scalars then do salt as well as temp
    if (ivar == BdyVars::t) {
        ncomp = 2;
    } else {
        ncomp = 1;
    }

    // This must be true for the logic below to work
    AMREX_ALWAYS_ASSERT(Temp_comp == 0);
    AMREX_ALWAYS_ASSERT(Salt_comp == 1);

    const Real eps= Real(1.0e-20);
    const bool null_mf_calc = (!mf_calc.ok());

    // ------------------------------------------------------------------------
    // Tidal forcing (ROMS ADD_FSOBC / ADD_M2OBC, set_tides.F lines 704-770 and
    // 849-931): the tidal elevation and tidal currents are ADDED to the
    // free-surface and 2D momentum open-boundary values that were just
    // interpolated from the boundary file. Only zeta, ubar and vbar are affected.
    //
    // The addition is done where the boundary value is consumed rather than by
    // mutating the stored interpolated data, so that repeated calls within a step
    // cannot accumulate tide on tide.
    //
    // Index map (ROMS -> REMORA, both rho and u/v indices shift by -1):
    //   zeta_west (j)  += 0.5*(Etide(Istr-1,j) + Etide(Istr  ,j)) -> 0.5*(E(dlo.x-1,j)+E(dlo.x  ,j))
    //   zeta_east (j)  += 0.5*(Etide(Iend  ,j) + Etide(Iend+1,j)) -> 0.5*(E(dhi.x  ,j)+E(dhi.x+1,j))
    //   zeta_south(i)  += 0.5*(Etide(i,Jstr-1) + Etide(i,Jstr  )) -> 0.5*(E(i,dlo.y-1)+E(i,dlo.y  ))
    //   zeta_north(i)  += 0.5*(Etide(i,Jend  ) + Etide(i,Jend+1)) -> 0.5*(E(i,dhi.y  )+E(i,dhi.y+1))
    //   ubar_west (j)  += Utide(Istr  ,j)      -> Utide(dlo.x  ,j)   [western u face]
    //   ubar_east (j)  += Utide(Iend+1,j)      -> Utide(dhi.x+1,j)   [eastern u face]
    //   ubar_south(i)  += Utide(i,Jstr-1)      -> Utide(i,dlo.y-1)   [u row below domain]
    //   ubar_north(i)  += Utide(i,Jend+1)      -> Utide(i,dhi.y+1)   [u row above domain]
    //   vbar_west (j)  += Vtide(Istr-1,j)      -> Vtide(dlo.x-1,j)   [v column west of domain]
    //   vbar_east (j)  += Vtide(Iend+1,j)      -> Vtide(dhi.x+1,j)   [v column east of domain]
    //   vbar_south(i)  += Vtide(i,Jstr  )      -> Vtide(i,dlo.y  )   [southern v face]
    //   vbar_north(i)  += Vtide(i,Jend+1)      -> Vtide(i,dhi.y+1)   [northern v face]
    // where dlo/dhi are the bounds of the CELL-CENTERED domain.
    //
    // The Flather condition for ubar/vbar also consumes the boundary free surface,
    // which in ROMS is the same (tide-carrying) BOUNDARY%zeta_* array, so the same
    // zeta tide is added to bry_val_zeta.
    // ------------------------------------------------------------------------
    const bool do_tides   = solverChoice.use_tides;
    const bool tide_zeta  = do_tides && (bccomp == zeta_bc());
    const bool tide_ubar  = do_tides && (bccomp == ubar_bc());
    const bool tide_vbar  = do_tides && (bccomp == vbar_bc());
    const bool add_tides  = tide_zeta || tide_ubar || tide_vbar;

    // Bounds of the cell-centered domain (note: `domain` above has been converted
    // to the nodality of the variable being filled)
    const auto& dlo = amrex::lbound(geom[lev].Domain());
    const auto& dhi = amrex::ubound(geom[lev].Domain());

    for (int icomp = 0; icomp < ncomp; icomp++) // This is to do both temp and salt if doing scalars
    {
        // If we're doing zeta, ubar, or vbar, then calc_arr only has a single component
        // corresponding to the component to be used in calculating the boundary
        // value. Since we access icomp + icomp_to_fill_calc, we need icomp_to_fill_calc to be zero.
        // If it's another variable, either we aren't using calc_arr
        // or the components correspond to salt, temp, etc so we leave it as is.
        int icomp_to_fill_calc = (bccomp == zeta_bc() || bccomp == ubar_bc() ||
                              bccomp == vbar_bc()) ? 0 : icomp_to_fill;

        boundary_series[lev][ivar+icomp]->update_interpolated_to_time(time);

        const auto& bdatxlo = boundary_series[lev][ivar+icomp]->xlo_dat_interp.const_array();
        const auto& bdatxhi = boundary_series[lev][ivar+icomp]->xhi_dat_interp.const_array();
        const auto& bdatylo = boundary_series[lev][ivar+icomp]->ylo_dat_interp.const_array();
        const auto& bdatyhi = boundary_series[lev][ivar+icomp]->yhi_dat_interp.const_array();

        const auto& bx_bdatxlo = boundary_series[lev][ivar+icomp]->xlo_dat_interp.box();
        const auto& bx_bdatxhi = boundary_series[lev][ivar+icomp]->xhi_dat_interp.box();
        const auto& bx_bdatylo = boundary_series[lev][ivar+icomp]->ylo_dat_interp.box();
        const auto& bx_bdatyhi = boundary_series[lev][ivar+icomp]->yhi_dat_interp.box();

        if (domain_bcs_type[bccomp+icomp].lo(0) == REMORABCType::flather ||
            domain_bcs_type[bccomp+icomp].hi(0) == REMORABCType::flather ||
            domain_bcs_type[bccomp+icomp].lo(1) == REMORABCType::flather ||
            domain_bcs_type[bccomp+icomp].hi(1) == REMORABCType::flather) {
            boundary_series[lev][BdyVars::zeta]->update_interpolated_to_time(time);
        }
        const auto& bdatxlo_zeta = domain_bcs_type[bccomp+icomp].lo(0) == REMORABCType::flather ?
                                   boundary_series[lev][BdyVars::zeta]->xlo_dat_interp.const_array() : Array4<Real>();
        const auto& bdatxhi_zeta = domain_bcs_type[bccomp+icomp].hi(0) == REMORABCType::flather ?
                                   boundary_series[lev][BdyVars::zeta]->xhi_dat_interp.const_array() : Array4<Real>();
        const auto& bdatylo_zeta = domain_bcs_type[bccomp+icomp].lo(1) == REMORABCType::flather ?
                                   boundary_series[lev][BdyVars::zeta]->ylo_dat_interp.const_array() : Array4<Real>();
        const auto& bdatyhi_zeta = domain_bcs_type[bccomp+icomp].hi(1) == REMORABCType::flather ?
                                   boundary_series[lev][BdyVars::zeta]->yhi_dat_interp.const_array() : Array4<Real>();

        const bool apply_west  = (domain_bcs_type[bccomp+icomp].lo(0) == REMORABCType::clamped) ||
                                 (domain_bcs_type[bccomp+icomp].lo(0) == REMORABCType::flather) ||
                                 (domain_bcs_type[bccomp+icomp].lo(0) == REMORABCType::chapman) ||
                                 (domain_bcs_type[bccomp+icomp].lo(0) == REMORABCType::orlanski_rad_nudge);
        const bool apply_east  = (domain_bcs_type[bccomp+icomp].hi(0) == REMORABCType::clamped) ||
                                 (domain_bcs_type[bccomp+icomp].hi(0) == REMORABCType::flather) ||
                                 (domain_bcs_type[bccomp+icomp].hi(0) == REMORABCType::chapman) ||
                                 (domain_bcs_type[bccomp+icomp].hi(0) == REMORABCType::orlanski_rad_nudge);
        const bool apply_south = (domain_bcs_type[bccomp+icomp].lo(1) == REMORABCType::clamped) ||
                                 (domain_bcs_type[bccomp+icomp].lo(1) == REMORABCType::flather) ||
                                 (domain_bcs_type[bccomp+icomp].lo(1) == REMORABCType::chapman) ||
                                 (domain_bcs_type[bccomp+icomp].lo(1) == REMORABCType::orlanski_rad_nudge);
        const bool apply_north = (domain_bcs_type[bccomp+icomp].hi(1) == REMORABCType::clamped) ||
                                 (domain_bcs_type[bccomp+icomp].hi(1) == REMORABCType::flather) ||
                                 (domain_bcs_type[bccomp+icomp].hi(1) == REMORABCType::chapman) ||
                                 (domain_bcs_type[bccomp+icomp].hi(1) == REMORABCType::orlanski_rad_nudge);

        const bool cell_centered = (mf_index_type[0] == 0 and mf_index_type[1] == 0);

        const Real obcfac = solverChoice.obcfac;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        // Currently no tiling in order to get the logic right
        for (MFIter mfi(mf_to_fill,false); mfi.isValid(); ++mfi)
        {
            Box mf_box(mf_to_fill[mfi.index()].box());

            // Compute intersections of the FAB to be filled and the bdry data boxes
            Box xlo = bx_bdatxlo & mf_box;
            Box xhi = bx_bdatxhi & mf_box;
            Box ylo = bx_bdatylo & mf_box;
            Box yhi = bx_bdatyhi & mf_box;

            xlo.setSmall(0,lbound(mf_box).x);
            xhi.setBig  (0,ubound(mf_box).x);
            ylo.setSmall(1,lbound(mf_box).y);
            yhi.setBig  (1,ubound(mf_box).y);

            Box xlo_ylo = xlo & ylo;
            Box xlo_yhi = xlo & yhi;
            Box xhi_ylo = xhi & ylo;
            Box xhi_yhi = xhi & yhi;

            Box xlo_edge = xlo; xlo_edge.setSmall(0,ubound(xlo).x); xlo_edge.setBig(0,ubound(xlo).x);
            Box xhi_edge = xhi; xhi_edge.setSmall(0,lbound(xhi).x); xhi_edge.setBig(0,lbound(xhi).x);
            Box ylo_edge = ylo; ylo_edge.setSmall(1,ubound(ylo).y); ylo_edge.setBig(1,ubound(ylo).y);
            Box yhi_edge = yhi; yhi_edge.setSmall(1,lbound(yhi).y); yhi_edge.setBig(1,lbound(yhi).y);

            Box xlo_ghost = xlo; xlo_ghost.setBig(0,ubound(xlo).x-1);
            Box xhi_ghost = xhi; xhi_ghost.setSmall(0,lbound(xhi).x+1);
            Box ylo_ghost = ylo; ylo_ghost.setBig(1,ubound(ylo).y-1);
            Box yhi_ghost = yhi; yhi_ghost.setSmall(1,lbound(yhi).y+1);

            // The box arithmetic here is settled, not guessed:
            // NCTimeSeriesBoundary.cpp:129 builds yhi_bx as the single row
            // y = hi[1]+1, so yhi_edge is exactly dom_hi.y+1 -- ROMS's
            // t(i,Jend+1) -- and yhi_ghost is everything beyond it.
            const Array4<Real>& dest_arr = mf_to_fill.array(mfi);
            const Array4<const Real>& mask_arr = mf_mask.array(mfi);
            const Array4<const Real>& calc_arr = (!null_mf_calc) ? mf_calc.array(mfi) : Array4<amrex::Real>();
            const Array4<const Real>& h_arr = vec_h[lev]->const_array(mfi);
            const Array4<const Real>& zeta_arr = vec_zeta[lev]->const_array(mfi);
            const Array4<const Real>& pm = vec_pm[lev]->const_array(mfi);
            const Array4<const Real>& pn = vec_pn[lev]->const_array(mfi);

            const Array4<const Real>& msku = vec_msku[lev]->const_array(mfi);
            const Array4<const Real>& mskv = vec_mskv[lev]->const_array(mfi);

            // Tidal elevation / currents for this box (zero unless remora.tides)
            const Array4<const Real>& Etide = vec_Etide[lev]->const_array(mfi);
            const Array4<const Real>& Utide = vec_Utide[lev]->const_array(mfi);
            const Array4<const Real>& Vtide = vec_Vtide[lev]->const_array(mfi);

            // ROMS PRESS_COMPENSATE (u2dbc_im.F/v2dbc_im.F): the Flather
            // condition compares interior zeta corrected by the inverse
            // barometer, fac*(Pair_a+Pair_b - 2*OneAtm), against the
            // (IB-free) boundary zeta. Pair is in millibar.
            const bool press_comp = solverChoice.press_compensate;
            const Real fac_pc    = Real(100.0)/(g*solverChoice.rho0);
            const Real OneAtm_pc = Real(1013.25);
            const Array4<const Real> Pair_bc = press_comp ?
                vec_Pair[lev]->const_array(mfi) : Array4<const Real>{};

            const Array4<const Real> nudg_coeff_out = vec_nudg_coeff[bdy_var_type][lev]->const_array(mfi);

            //
            // We are inside a loop over components so we do one at a time here
            //
            Vector<BCRec> bcrs(1);
            amrex::setBC(mf_box, domain, bccomp+icomp, 0, 1, domain_bcs_type, bcrs);

            // xlo: ori = 0
            // ylo: ori = 1
            // zlo: ori = 2
            // xhi: ori = 3
            // yhi: ori = 4
            // zhi: ori = 5

            auto bcr = bcrs[0];

            // Even though we don't loop over xlo itself, this is the right condition to check, since xlo_edge will always be the same for each grid,
            // but if the grid doesn't include the low x-boundary, the xlo box will be invalid and the execution will be skipped.
            if (!xlo.isEmpty() && apply_west) {
                ParallelFor(grow(xlo_edge,IntVect(0,-1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    // ROMS set_tides.f90:704-719 (zeta_west), 849-868 (ubar/vbar_west)
                    Real tide_zeta_val = add_tides ?
                        Real(0.5) * (Etide(dlo.x-1,j,0) + Etide(dlo.x,j,0)) : zero;
                    Real tide_val = tide_zeta ? tide_zeta_val :
                                   (tide_ubar ? Utide(dlo.x  ,j,0) :
                                   (tide_vbar ? Vtide(dlo.x-1,j,0) : zero));

                    Real bry_val = bdatxlo(ubound(xlo).x,j,k,0) + tide_val;
                    if (bcr.lo(0) == REMORABCType::clamped) {
                        dest_arr(i,j,k,icomp+icomp_to_fill) = bry_val * mask_arr(i,j,0);
                    } else if (bcr.lo(0) == REMORABCType::flather) {
                        Real bry_val_zeta = bdatxlo_zeta(ubound(xlo).x-1,j,k,0) + tide_zeta_val;
                        Real cff = one / (Real(0.5) * (h_arr(dom_lo.x-1,j,0) + zeta_arr(dom_lo.x-1,j,0,icomp_calc)
                                                     + h_arr(dom_lo.x,j,0) + zeta_arr(dom_lo.x,j,0,icomp_calc)));
                        Real Cx = std::sqrt(g * cff);
                        Real zsum = zeta_arr(dom_lo.x-1,j,0,icomp_calc) + zeta_arr(dom_lo.x,j,0,icomp_calc);
                        if (press_comp) {
                            zsum += fac_pc * (Pair_bc(dom_lo.x-1,j,0) + Pair_bc(dom_lo.x,j,0) - two*OneAtm_pc);
                        }
                        dest_arr(i,j,k,icomp+icomp_to_fill) = (bry_val
                                - Cx * (Real(0.5) * zsum
                                    - bry_val_zeta)) * mask_arr(i,j,0);
                    } else if (bcr.lo(0) == REMORABCType::chapman) {
                        Real cff = dt_calc * Real(0.5) * (pm(dom_lo.x,j-mf_index_type[1],0) + pm(dom_lo.x,j,0));
                        Real cff1 = std::sqrt(g * Real(0.5) * (h_arr(dom_lo.x,j-mf_index_type[1],0)
                                    + zeta_arr(dom_lo.x,j-mf_index_type[1],0,icomp_calc) + h_arr(dom_lo.x,j,0)
                                    + zeta_arr(dom_lo.x,j,0,icomp_calc)));
                        Real Cx = cff * cff1;
                        Real cff2 = one / (one + Cx);
                        dest_arr(i,j,k,icomp+icomp_to_fill) = cff2 * (dest_arr(dom_lo.x-1,j,k,icomp_calc)
                                + Cx * dest_arr(dom_lo.x,j,k,icomp+icomp_to_fill)) * mask_arr(i,j,0);
                    } else if (bcr.lo(0) == REMORABCType::orlanski_rad_nudge) {
                        Real grad_lo_im1   = (calc_arr(dom_lo.x+mf_index_type[0]-1,j  ,k,icomp+icomp_to_fill_calc) - calc_arr(dom_lo.x-1+mf_index_type[0],j-1,k,icomp+icomp_to_fill_calc));
                        Real grad_lo       = (calc_arr(dom_lo.x+mf_index_type[0]  ,j  ,k,icomp+icomp_to_fill_calc) - calc_arr(dom_lo.x  +mf_index_type[0],j-1,k,icomp+icomp_to_fill_calc));
                        Real grad_lo_imjp1 = (calc_arr(dom_lo.x+mf_index_type[0]-1,j+1,k,icomp+icomp_to_fill_calc) - calc_arr(dom_lo.x-1+mf_index_type[0],j  ,k,icomp+icomp_to_fill_calc));
                        Real grad_lo_jp1   = (calc_arr(dom_lo.x+mf_index_type[0]  ,j+1,k,icomp+icomp_to_fill_calc) - calc_arr(dom_lo.x  +mf_index_type[0],j  ,k,icomp+icomp_to_fill_calc));
                        if (cell_centered) {
                                // ROMS masks each gradient with the vmask at
                                // ITS OWN (i,j), not at the boundary point:
                                // t3dbc_im.f90 west uses vmask(Istr-1,j) for
                                // grad(Istr-1,j) and vmask(Istr,j) for
                                // grad(Istr,j), and the j+1 pair likewise.
                                // A single mskv(i,j) is right only for the
                                // first of the four.
                                const int ib = dom_lo.x-1+mf_index_type[0];
                                const int ii = dom_lo.x  +mf_index_type[0];
                                grad_lo_im1   *= mskv(ib,j  ,0);
                                grad_lo       *= mskv(ii,j  ,0);
                                grad_lo_imjp1 *= mskv(ib,j+1,0);
                                grad_lo_jp1   *= mskv(ii,j+1,0);
                        }
                        Real dTdt = calc_arr(dom_lo.x+mf_index_type[0],j,k,icomp+icomp_to_fill_calc) - dest_arr(dom_lo.x+mf_index_type[0]  ,j,k,icomp+icomp_to_fill);
                        Real dTdx = dest_arr(dom_lo.x+mf_index_type[0],j,k,icomp+icomp_to_fill) - dest_arr(dom_lo.x+mf_index_type[0]+1,j,k,icomp+icomp_to_fill);
                        Real tau;
                        Real nudg_coeff_out_local = (nudg_coeff_out(i-mf_index_type[0],j-mf_index_type[1],k) +
                                                     nudg_coeff_out(i,j,k)) * Real(0.5);
                        if (dTdt*dTdx < zero) {
                            tau = nudg_coeff_out_local * obcfac * dt_calc;
                            dTdt = zero;
                        } else {
                            tau = nudg_coeff_out_local * dt_calc;
                        }
                        Real dTde = (dTdt * (grad_lo+grad_lo_jp1) > zero) ? grad_lo : grad_lo_jp1;
                        Real cff = std::max(dTdx*dTdx+dTde*dTde,eps);
                        Real Cx = dTdt * dTdx;
                        // Tangential radiation term. ROMS t3dbc_im.f90:629-639
                        // (and u3dbc/v3dbc identically) carries
                        //   -MAX(Ce,0)*grad(ghost,j) - MIN(Ce,0)*grad(ghost,j+1)
                        // in the numerator, with Ce clamped to +/-cff. Omitting
                        // it leaves the boundary purely normally-radiating: it
                        // agrees while the field is along-boundary uniform
                        // (dTdt ~ 0 => Ce ~ 0) and diverges as soon as a
                        // gradient develops along the edge.
                        Real Ce = std::min(cff, std::max(dTdt*dTde, -cff));
                        dest_arr(i,j,k,icomp+icomp_to_fill) = (cff * calc_arr(dom_lo.x-1+mf_index_type[0],j,k,icomp+icomp_to_fill_calc) + Cx * dest_arr(dom_lo.x+mf_index_type[0],j,k,icomp+icomp_to_fill)
                                - std::max(Ce,Real(0.0)) * grad_lo_im1
                                - std::min(Ce,Real(0.0)) * grad_lo_imjp1) / (cff+Cx);
                        dest_arr(i,j,k,icomp+icomp_to_fill) = mask_arr(i,j,0) * (dest_arr(dom_lo.x-1+mf_index_type[0],j,k,icomp+icomp_to_fill) + tau * (bry_val - calc_arr(dom_lo.x-1+mf_index_type[0],j,k,icomp+icomp_to_fill_calc)));
                    }
                });
                ParallelFor(grow(xlo_ghost,IntVect(0,-1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    dest_arr(i,j,k,icomp+icomp_to_fill) = dest_arr(ubound(xlo).x,j,k,icomp+icomp_to_fill);
                });
            }

            // See comment on xlo
            if (!xhi.isEmpty() && apply_east) {
                ParallelFor(grow(xhi_edge,IntVect(0,-1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    // ROMS set_tides.f90:721-736 (zeta_east), 870-889 (ubar/vbar_east)
                    Real tide_zeta_val = add_tides ?
                        Real(0.5) * (Etide(dhi.x,j,0) + Etide(dhi.x+1,j,0)) : zero;
                    Real tide_val = tide_zeta ? tide_zeta_val :
                                   (tide_ubar ? Utide(dhi.x+1,j,0) :
                                   (tide_vbar ? Vtide(dhi.x+1,j,0) : zero));

                    Real bry_val = bdatxhi(lbound(xhi).x,j,k,0) + tide_val;
                    if (bcr.hi(0) == REMORABCType::clamped) {
                        dest_arr(i,j,k,icomp+icomp_to_fill) = bry_val * mask_arr(i,j,0);
                    } else if (bcr.hi(0) == REMORABCType::flather) {
                        Real bry_val_zeta = bdatxhi_zeta(lbound(xhi).x,j,k,0) + tide_zeta_val;
                        Real cff = one / (Real(0.5) * (h_arr(dom_hi.x-1,j,0) + zeta_arr(dom_hi.x-1,j,0,icomp_calc)
                                                     + h_arr(dom_hi.x,j,0) + zeta_arr(dom_hi.x,j,0,icomp_calc)));
                        Real Cx = std::sqrt(g * cff);
                        Real zsum = zeta_arr(dom_hi.x-1,j,0,icomp_calc) + zeta_arr(dom_hi.x,j,0,icomp_calc);
                        if (press_comp) {
                            zsum += fac_pc * (Pair_bc(dom_hi.x-1,j,0) + Pair_bc(dom_hi.x,j,0) - two*OneAtm_pc);
                        }
                        dest_arr(i,j,k,icomp+icomp_to_fill) = (bry_val
                                + Cx * (Real(0.5) * zsum
                                    - bry_val_zeta)) * mask_arr(i,j,0);
                    } else if (bcr.hi(0) == REMORABCType::chapman) {
                        Real cff = dt_calc * Real(0.5) * (pm(dom_hi.x,j-mf_index_type[1],0) + pm(dom_hi.x,j,0));
                        Real cff1 = std::sqrt(g * Real(0.5) * (h_arr(dom_hi.x,j-mf_index_type[1],0)
                                    + zeta_arr(dom_hi.x,j-mf_index_type[1],0,icomp_calc) + h_arr(dom_hi.x,j,0)
                                    + zeta_arr(dom_hi.x,j,0,icomp_calc)));
                        Real Cx = cff * cff1;
                        Real cff2 = one / (one + Cx);
                        dest_arr(i,j,k,icomp+icomp_to_fill) = cff2 * (dest_arr(dom_hi.x+1,j,k,icomp_calc)
                                + Cx * dest_arr(dom_hi.x,j,k,icomp+icomp_to_fill)) * mask_arr(i,j,0);
                    } else if (bcr.hi(0) == REMORABCType::orlanski_rad_nudge) {
                        Real grad_hi      = (calc_arr(dom_hi.x-mf_index_type[0]  ,j  ,k,icomp+icomp_to_fill_calc) - calc_arr(dom_hi.x-mf_index_type[0]  ,j-1,k,icomp+icomp_to_fill_calc));
                        Real grad_hi_ip1  = (calc_arr(dom_hi.x-mf_index_type[0]+1,j  ,k,icomp+icomp_to_fill_calc) - calc_arr(dom_hi.x-mf_index_type[0]+1,j-1,k,icomp+icomp_to_fill_calc));
                        Real grad_hi_jp1  = (calc_arr(dom_hi.x-mf_index_type[0]  ,j+1,k,icomp+icomp_to_fill_calc) - calc_arr(dom_hi.x-mf_index_type[0]  ,j  ,k,icomp+icomp_to_fill_calc));
                        Real grad_hi_ijp1 = (calc_arr(dom_hi.x-mf_index_type[0]+1,j+1,k,icomp+icomp_to_fill_calc) - calc_arr(dom_hi.x-mf_index_type[0]+1,j  ,k,icomp+icomp_to_fill_calc));
                        if (cell_centered) {
                            // See comment on xlo: each gradient takes the
                            // vmask at its own (i,j).
                            const int ii = dom_hi.x-mf_index_type[0];
                            const int ib = dom_hi.x-mf_index_type[0]+1;
                            grad_hi      *= mskv(ii,j  ,0);
                            grad_hi_ip1  *= mskv(ib,j  ,0);
                            grad_hi_jp1  *= mskv(ii,j+1,0);
                            grad_hi_ijp1 *= mskv(ib,j+1,0);
                        }
                        Real dTdt = calc_arr(dom_hi.x-mf_index_type[0],j,k,icomp+icomp_to_fill_calc) - dest_arr(dom_hi.x-mf_index_type[0]  ,j,k,icomp+icomp_to_fill);
                        Real dTdx = dest_arr(dom_hi.x-mf_index_type[0],j,k,icomp+icomp_to_fill) - dest_arr(dom_hi.x-mf_index_type[0]-1,j,k,icomp+icomp_to_fill);
                        Real tau;
                        Real nudg_coeff_out_local = (nudg_coeff_out(i-mf_index_type[0],j-mf_index_type[1],k) +
                                                     nudg_coeff_out(i,j,k)) * Real(0.5);
                        if (dTdt*dTdx < zero) {
                            tau = nudg_coeff_out_local * obcfac * dt_calc;
                            dTdt = zero;
                        } else {
                            tau = nudg_coeff_out_local * dt_calc;
                        }
                        if (dTdt * dTdx < zero) dTdt = zero;
                        Real dTde = (dTdt * (grad_hi + grad_hi_jp1) > zero) ? grad_hi : grad_hi_jp1;
                        Real cff = std::max(dTdx*dTdx + dTde*dTde,eps);
                        Real Cx = dTdt * dTdx;
                        Real Ce = std::min(cff, std::max(dTdt*dTde, -cff));
                        dest_arr(i,j,k,icomp+icomp_to_fill) = (cff * calc_arr(dom_hi.x+1-mf_index_type[0],j,k,icomp+icomp_to_fill_calc) + Cx * dest_arr(dom_hi.x-mf_index_type[0],j,k,icomp+icomp_to_fill)
                                - std::max(Ce,Real(0.0)) * grad_hi_ip1
                                - std::min(Ce,Real(0.0)) * grad_hi_ijp1) * mask_arr(i,j,0) / (cff+Cx);
                        dest_arr(i,j,k,icomp+icomp_to_fill) = mask_arr(i,j,0) * (dest_arr(dom_hi.x+1-mf_index_type[0],j,k,icomp+icomp_to_fill) + tau * (bry_val - calc_arr(dom_hi.x+1-mf_index_type[0],j,k,icomp+icomp_to_fill_calc)));
                    }
                });
                ParallelFor(grow(xhi_ghost,IntVect(0,-1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    dest_arr(i,j,k,icomp+icomp_to_fill) = dest_arr(lbound(xhi).x,j,k,icomp+icomp_to_fill);
                });
            }

            // See comment on xlo
            if (!ylo.isEmpty() && apply_south) {
                ParallelFor(grow(ylo_edge,IntVect(-1,0,0)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    // ROMS set_tides.f90:738-753 (zeta_south), 891-910 (ubar/vbar_south)
                    Real tide_zeta_val = add_tides ?
                        Real(0.5) * (Etide(i,dlo.y-1,0) + Etide(i,dlo.y,0)) : zero;
                    Real tide_val = tide_zeta ? tide_zeta_val :
                                   (tide_ubar ? Utide(i,dlo.y-1,0) :
                                   (tide_vbar ? Vtide(i,dlo.y  ,0) : zero));

                    Real bry_val = bdatylo(i,ubound(ylo).y,k,0) + tide_val;
                    if (bcr.lo(1) == REMORABCType::clamped) {
                        dest_arr(i,j,k,icomp+icomp_to_fill) = bry_val * mask_arr(i,j,0);
                    } else if (bcr.lo(1) == REMORABCType::flather) {
                        Real bry_val_zeta = bdatylo_zeta(i,ubound(ylo).y-1,k,0) + tide_zeta_val;
                        Real cff = one / (Real(0.5) * (h_arr(i,dom_lo.y-1,0) + zeta_arr(i,dom_lo.y-1,0,icomp_calc)
                                                     + h_arr(i,dom_lo.y,0) + zeta_arr(i,dom_lo.y,0,icomp_calc)));
                        Real Ce = std::sqrt(g * cff);
                        Real zsum = zeta_arr(i,dom_lo.y-1,0,icomp_calc) + zeta_arr(i,dom_lo.y,0,icomp_calc);
                        if (press_comp) {
                            zsum += fac_pc * (Pair_bc(i,dom_lo.y-1,0) + Pair_bc(i,dom_lo.y,0) - two*OneAtm_pc);
                        }
                        dest_arr(i,j,k,icomp+icomp_to_fill) = (bry_val
                                - Ce * (Real(0.5) * zsum
                                    - bry_val_zeta)) * mask_arr(i,j,0);
                    } else if (bcr.lo(1) == REMORABCType::chapman) {
                        Real cff = dt_calc * Real(0.5) * (pn(i-mf_index_type[0],dom_lo.y,0) + pn(i,dom_lo.y,0));
                        Real cff1 = std::sqrt(g * Real(0.5) * (h_arr(i-mf_index_type[0],dom_lo.y,0) +
                                    zeta_arr(i-mf_index_type[0],dom_lo.y,0,icomp_calc) + h_arr(i,dom_lo.y,0)
                                    + zeta_arr(i,dom_lo.y,0,icomp_calc)));
                        Real Ce = cff * cff1;
                        Real cff2 = one / (one + Ce);
                        dest_arr(i,j,k,icomp+icomp_to_fill) = cff2 * (dest_arr(i,dom_lo.y-1,k,icomp_calc)
                                + Ce * dest_arr(i,dom_lo.y,k,icomp+icomp_to_fill)) * mask_arr(i,j,0);
                    } else if (bcr.lo(1) == REMORABCType::orlanski_rad_nudge) {
                        Real grad_lo       = (calc_arr(i  ,dom_lo.y+mf_index_type[1],  k,icomp+icomp_to_fill_calc) - calc_arr(i-1,dom_lo.y+mf_index_type[1]  ,k,icomp+icomp_to_fill_calc));
                        Real grad_lo_jm1   = (calc_arr(i  ,dom_lo.y+mf_index_type[1]-1,k,icomp+icomp_to_fill_calc) - calc_arr(i-1,dom_lo.y+mf_index_type[1]-1,k,icomp+icomp_to_fill_calc));
                        Real grad_lo_ip1   = (calc_arr(i+1,dom_lo.y+mf_index_type[1]  ,k,icomp+icomp_to_fill_calc) - calc_arr(i  ,dom_lo.y+mf_index_type[1]  ,k,icomp+icomp_to_fill_calc));
                        Real grad_lo_ipjm1 = (calc_arr(i+1,dom_lo.y+mf_index_type[1]-1,k,icomp+icomp_to_fill_calc) - calc_arr(i  ,dom_lo.y+mf_index_type[1]-1,k,icomp+icomp_to_fill_calc));
                        if (cell_centered) {
                            // See comment on xlo: each gradient takes the
                            // umask at its own (i,j).
                            const int jj = dom_lo.y+mf_index_type[1];
                            const int jb = dom_lo.y+mf_index_type[1]-1;
                            grad_lo       *= msku(i  ,jj,0);
                            grad_lo_jm1   *= msku(i  ,jb,0);
                            grad_lo_ip1   *= msku(i+1,jj,0);
                            grad_lo_ipjm1 *= msku(i+1,jb,0);
                        }
                        Real dTdt = calc_arr(i,dom_lo.y+mf_index_type[1],k,icomp+icomp_to_fill_calc) - dest_arr(i,dom_lo.y  +mf_index_type[1],k,icomp+icomp_to_fill);
                        Real dTde = dest_arr(i,dom_lo.y+mf_index_type[1],k,icomp+icomp_to_fill) - dest_arr(i,dom_lo.y+1+mf_index_type[1],k,icomp+icomp_to_fill);
                        Real tau;
                        Real nudg_coeff_out_local = (nudg_coeff_out(i-mf_index_type[0],j-mf_index_type[1],k) +
                                                     nudg_coeff_out(i,j,k)) * Real(0.5);
                        if (dTdt*dTde < zero) {
                            tau = nudg_coeff_out_local * obcfac * dt_calc;
                            dTdt = zero;
                        } else {
                            tau = nudg_coeff_out_local * dt_calc;
                        }
                        if (dTdt * dTde < zero) dTdt = zero;
                        Real dTdx = (dTdt * (grad_lo + grad_lo_ip1) > zero) ? grad_lo : grad_lo_ip1;
                        Real cff = std::max(dTdx*dTdx + dTde*dTde, eps);
                        Real Ce = dTdt*dTde;
                        // Tangential term, ROMS t3dbc_im.f90:884-895. On the
                        // y-normal edges the roles swap: Ce is the unclamped
                        // normal coefficient and Cx is the clamped tangential
                        // one, applied to the ghost-ROW gradients.
                        Real Cx = std::min(cff, std::max(dTdt*dTdx, -cff));
                        dest_arr(i,j,k,icomp+icomp_to_fill) = (cff * calc_arr(i,dom_lo.y-1+mf_index_type[1],k,icomp+icomp_to_fill_calc) + Ce * dest_arr(i,dom_lo.y+mf_index_type[1],k,icomp+icomp_to_fill)
                                - std::max(Cx,Real(0.0)) * grad_lo_jm1
                                - std::min(Cx,Real(0.0)) * grad_lo_ipjm1) / (cff+Ce);
                        dest_arr(i,j,k,icomp+icomp_to_fill) = mask_arr(i,j,0) * (dest_arr(i,dom_lo.y-1+mf_index_type[1],k,icomp+icomp_to_fill) + tau * (bry_val - calc_arr(i,dom_lo.y-1+mf_index_type[1],k,icomp+icomp_to_fill_calc)));
                    }
                });
                ParallelFor(grow(ylo_ghost,IntVect(-1,0,0)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    dest_arr(i,j,k,icomp+icomp_to_fill) = dest_arr(i,ubound(ylo).y,k,icomp+icomp_to_fill);
                });
            }

            // See comment on xlo
            if (!yhi.isEmpty() && apply_north) {
                ParallelFor(grow(yhi_edge,IntVect(-1,0,0)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    // ROMS set_tides.f90:755-770 (zeta_north), 912-931 (ubar/vbar_north)
                    Real tide_zeta_val = add_tides ?
                        Real(0.5) * (Etide(i,dhi.y,0) + Etide(i,dhi.y+1,0)) : zero;
                    Real tide_val = tide_zeta ? tide_zeta_val :
                                   (tide_ubar ? Utide(i,dhi.y+1,0) :
                                   (tide_vbar ? Vtide(i,dhi.y+1,0) : zero));

                    Real bry_val = bdatyhi(i,lbound(yhi).y,k,0) + tide_val;
                    if (bcr.hi(1) == REMORABCType::clamped) {
                        dest_arr(i,j,k,icomp+icomp_to_fill) = bry_val * mask_arr(i,j,0);
                    } else if (bcr.hi(1) == REMORABCType::flather) {
                        Real bry_val_zeta = bdatyhi_zeta(i,lbound(yhi).y,k,0) + tide_zeta_val;
                        Real cff = one / (Real(0.5) * (h_arr(i,dom_hi.y-1,0) + zeta_arr(i,dom_hi.y-1,0,icomp_calc)
                                                     + h_arr(i,dom_hi.y,0) + zeta_arr(i,dom_hi.y,0,icomp_calc)));
                        Real Ce = std::sqrt(g * cff);
                        Real zsum = zeta_arr(i,dom_hi.y-1,0,icomp_calc) + zeta_arr(i,dom_hi.y,0,icomp_calc);
                        if (press_comp) {
                            zsum += fac_pc * (Pair_bc(i,dom_hi.y-1,0) + Pair_bc(i,dom_hi.y,0) - two*OneAtm_pc);
                        }
                        dest_arr(i,j,k,icomp+icomp_to_fill) = (bry_val
                                + Ce * (Real(0.5) * zsum
                                    - bry_val_zeta)) * mask_arr(i,j,0);
                    } else if (bcr.hi(1) == REMORABCType::chapman) {
                        Real cff = dt_calc * Real(0.5) * (pn(i-mf_index_type[0],dom_hi.y,0) + pn(i,dom_hi.y,0));
                        Real cff1 = std::sqrt(g * Real(0.5) * (h_arr(i-mf_index_type[0],dom_hi.y,0)
                                                          + zeta_arr(i-mf_index_type[0],dom_hi.y,0,icomp_calc) +
                                                            h_arr(i,dom_hi.y,0) + zeta_arr(i,dom_hi.y,0,icomp_calc)));
                        Real Ce = cff * cff1;
                        Real cff2 = one / (one + Ce);
                        dest_arr(i,j,k,icomp+icomp_to_fill) = cff2 * (dest_arr(i,dom_hi.y+1,k,icomp_calc)
                                + Ce * dest_arr(i,dom_hi.y,k,icomp+icomp_to_fill)) * mask_arr(i,j,0);
                    } else if (bcr.hi(1) == REMORABCType::orlanski_rad_nudge) {
                        Real grad_hi      = calc_arr(i  ,dom_hi.y-mf_index_type[1]  ,k,icomp+icomp_to_fill_calc) - calc_arr(i-1,dom_hi.y-mf_index_type[1]  ,k,icomp+icomp_to_fill_calc);
                        Real grad_hi_jp1  = calc_arr(i  ,dom_hi.y-mf_index_type[1]+1,k,icomp+icomp_to_fill_calc) - calc_arr(i-1,dom_hi.y-mf_index_type[1]+1,k,icomp+icomp_to_fill_calc);
                        Real grad_hi_ip1  = calc_arr(i+1,dom_hi.y-mf_index_type[1]  ,k,icomp+icomp_to_fill_calc) - calc_arr(i  ,dom_hi.y-mf_index_type[1]  ,k,icomp+icomp_to_fill_calc);
                        Real grad_hi_ijp1 = calc_arr(i+1,dom_hi.y-mf_index_type[1]+1,k,icomp+icomp_to_fill_calc) - calc_arr(i  ,dom_hi.y-mf_index_type[1]+1,k,icomp+icomp_to_fill_calc);
                        if (cell_centered) {
                            // See comment on xlo: each gradient takes the
                            // umask at its own (i,j).
                            const int jj = dom_hi.y-mf_index_type[1];
                            const int jb = dom_hi.y-mf_index_type[1]+1;
                            grad_hi      *= msku(i  ,jj,0);
                            grad_hi_jp1  *= msku(i  ,jb,0);
                            grad_hi_ip1  *= msku(i+1,jj,0);
                            grad_hi_ijp1 *= msku(i+1,jb,0);
                        }
                        Real dTdt = calc_arr(i,dom_hi.y-mf_index_type[1],k,icomp+icomp_to_fill_calc) - dest_arr(i,dom_hi.y  -mf_index_type[1],k,icomp+icomp_to_fill);
                        Real dTde = dest_arr(i,dom_hi.y-mf_index_type[1],k,icomp+icomp_to_fill) - dest_arr(i,dom_hi.y-1-mf_index_type[1],k,icomp+icomp_to_fill);
                        Real tau;
                        Real nudg_coeff_out_local = (nudg_coeff_out(i-mf_index_type[0],j-mf_index_type[1],k) +
                                                     nudg_coeff_out(i,j,k)) * Real(0.5);
                        if (dTdt*dTde < zero) {
                            tau = nudg_coeff_out_local * obcfac * dt_calc;
                            dTdt = zero;
                        } else {
                            tau = nudg_coeff_out_local * dt_calc;
                        }
                        if (dTdt * dTde < zero) dTdt = zero;
                        Real dTdx = (dTdt * (grad_hi + grad_hi_ip1) > zero) ? grad_hi : grad_hi_ip1;
                        Real cff = std::max(dTdx*dTdx + dTde*dTde, eps);
                        Real Ce = dTdt*dTde;
                        Real Cx = std::min(cff, std::max(dTdt*dTdx, -cff));
                        dest_arr(i,j,k,icomp+icomp_to_fill) = (cff*calc_arr(i,dom_hi.y+1-mf_index_type[1],k,icomp+icomp_to_fill_calc) + Ce*dest_arr(i,dom_hi.y-mf_index_type[1],k,icomp+icomp_to_fill)
                                - std::max(Cx,Real(0.0)) * grad_hi_jp1
                                - std::min(Cx,Real(0.0)) * grad_hi_ijp1) * mask_arr(i,j,0) / (cff+Ce);
                        dest_arr(i,j,k,icomp+icomp_to_fill) = mask_arr(i,j,0) * (dest_arr(i,dom_hi.y+1-mf_index_type[1],k,icomp+icomp_to_fill) + tau * (bry_val - calc_arr(i,dom_hi.y+1-mf_index_type[1],k,icomp+icomp_to_fill_calc)));
                    }
                });
                ParallelFor(grow(yhi_ghost,IntVect(-1,0,0)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    dest_arr(i,j,k,icomp+icomp_to_fill) = dest_arr(i,lbound(yhi).y,k,icomp+icomp_to_fill);
                });

            }
            // If we've applied boundary conditions to either side, update the corner
            if (!xlo_ylo.isEmpty() && (apply_west || apply_south)) {
                ParallelFor(xlo_ylo, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    dest_arr(i,j,k,icomp+icomp_to_fill) = Real(0.5) * (dest_arr(i,dom_lo.y+mf_index_type[1],k,icomp+icomp_to_fill)
                                                               + dest_arr(dom_lo.x+mf_index_type[0],j,k,icomp+icomp_to_fill));
                });
            }
            if (!xlo_yhi.isEmpty() && (apply_west || apply_north)) {
                ParallelFor(xlo_yhi, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    dest_arr(i,j,k,icomp+icomp_to_fill) = Real(0.5) * (dest_arr(i,dom_hi.y-mf_index_type[1],k,icomp+icomp_to_fill)
                                                               + dest_arr(dom_lo.x+mf_index_type[0],j,k,icomp+icomp_to_fill));
                });
            }
            if (!xhi_ylo.isEmpty() && (apply_east || apply_south)) {
                ParallelFor(xhi_ylo, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    dest_arr(i,j,k,icomp+icomp_to_fill) = Real(0.5) * (dest_arr(i,dom_lo.y+mf_index_type[1],k,icomp+icomp_to_fill)
                                                               + dest_arr(dom_hi.x-mf_index_type[0],j,k,icomp+icomp_to_fill));
                });
            }
            if (!xhi_yhi.isEmpty() && (apply_east || apply_north)) {
                ParallelFor(xhi_yhi, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    dest_arr(i,j,k,icomp+icomp_to_fill) = Real(0.5) * (dest_arr(i,dom_hi.y-mf_index_type[1],k,icomp+icomp_to_fill)
                                                               + dest_arr(dom_hi.x-mf_index_type[0],j,k,icomp+icomp_to_fill));
                });
            }
        } // mfi
    } // icomp
}
#endif
