#include <cmath>
#include <REMORA_DataStruct.H>
#include <REMORA.H>
#include <REMORA_prob_common.H>

using namespace amrex;

/**
 * @param[in   ] lev     level to operate on
 */
void REMORA::set_weights (int /*lev*/) {

    Real gamma, scale;
    Real wsum, shift, cff;

    //HACK should possibly store fixed_ndtfast elsewhere
    int ndtfast=fixed_ndtfast_ratio>0 ? fixed_ndtfast_ratio : static_cast<int>(fixed_fast_dt / fixed_dt);

    //From mod_scalars
    Real Falpha = two;
    Real Fbeta = Real(4.0);
    Real Fgamma = Real(0.284);

    vec_weight1.resize(2*ndtfast+1);
    vec_weight2.resize(2*ndtfast+1);

    auto weight1 = vec_weight1.dataPtr();
    auto weight2 = vec_weight2.dataPtr();

//
//=======================================================================
//  Compute time-averaging filter for barotropic fields.
//=======================================================================
//
//  Initialize both sets of weights to zero.
//
    nfast=0;
    for(int i=1;i<=2*ndtfast;i++) {
        weight1[i-1]=zero;
        weight2[i-1]=zero;
    }
//
//-----------------------------------------------------------------------
//  Power-law shape filters.
//-----------------------------------------------------------------------
//
//  The power-law shape filters are given by:
//
//     F(xi)=xi^Falpha*(1-xi^Fbeta)-Fgamma*xi
//
//  where xi=scale*i/ndtfast; and scale, Falpha, Fbeta, Fgamma, and
//  normalization are chosen to yield the correct zeroth-order
//  (normalization), first-order (consistency), and second-order moments,
//  resulting in overall second-order temporal accuracy for time-averaged
//  barotropic motions resolved by baroclinic time step.
//
    scale=(Falpha+one)*(Falpha+Fbeta+one) /
        ((Falpha+two)*(Falpha+Fbeta+two)*Real(ndtfast));
    //
    //  Find center of gravity of the primary weighting shape function and
    //  iteratively adjust "scale" to place the  centroid exactly at
    //  "ndtfast".
    //
    gamma = Fgamma*max(zero, one-Real(10.0)/Real(ndtfast));

    for (int iter=1;iter<=16;iter++) {
        nfast=0;
        for(int i=1;i<=2*ndtfast;i++) {
            cff=scale*Real(i);

            weight1[i-1]=Real(pow(cff,Falpha)-pow(cff,(Falpha+Fbeta)))-gamma*cff;

            if (weight1[i-1] > zero) {
                nfast=i;
            }

            if ( (nfast>0) && (weight1[i-1] < zero) ) {
                weight1[i-1] = zero;
            }
        }
        wsum  = zero;
        shift = zero;
        for(int i=1;i<=nfast;i++) {
            wsum=wsum+weight1[i-1];
            shift=shift+weight1[i-1]*Real(i);
        }
        // GPU-PARITY: this is `scale = (scale*shift) / (wsum*ndtfast)`, NOT
        // `scale *= shift/(wsum*ndtfast)`. set_weights.F:95 writes
        //     scale=scale*shift/(wsum*REAL(ndtfast(ng),dp))
        // which Fortran evaluates left to right, so the numerator is rounded
        // once as a product before the division. The compound-assignment form
        // rounds the quotient first and then multiplies -- a different result,
        // and it is fed back into `scale` through sixteen iterations. Measured:
        // with the two forms compiled under matched IEEE flags, the compound
        // form left 116 of the 122 final weight fields differing from ROMS in
        // the last bits; this form leaves 0. See the note at the end of this
        // file.
        scale = scale*shift/(wsum*Real(ndtfast));
    }
//
//-----------------------------------------------------------------------
//  Post-processing of primary weights.
//-----------------------------------------------------------------------
//
//  Although it is assumed that the initial settings of the primary
//  weights has its center of gravity "reasonably close" to NDTFAST,
//  it may be not so according to the discrete rules of integration.
//  The following procedure is designed to put the center of gravity
//  exactly to NDTFAST by computing mismatch (NDTFAST-shift) and
//  applying basically an upstream advection of weights to eliminate
//  the mismatch iteratively. Once this procedure is complete primary
//  weights are normalized.
//
//  Find center of gravity of the primary weights and subsequently
//  calculate the mismatch to be compensated.
//
    for (int iter=1;iter<=ndtfast;iter++) {
        wsum  = zero;
        shift = zero;
        for(int i=1;i<=nfast;i++) {
            wsum=wsum+weight1[i-1];
            shift=shift+Real(i)*weight1[i-1];
        }
        shift=shift/wsum;
        cff=Real(ndtfast)-shift;
        //
        //  Apply advection step using either whole, or fractional shifts.
        //  Notice that none of the four loops here is reversible.
        //
        if (cff > one) {
            nfast=nfast+1;
            for (int i=nfast;i>=2;i--) {
                weight1[i-1]=weight1[i-1-1];
            }
            weight1[1-1] = zero;
        } else if (cff> zero) {
            wsum=one-cff;
            for (int i=nfast;i>=2;i--) {
                weight1[i-1]=wsum*weight1[i-1]+cff*weight1[i-1-1];
            }
            weight1[1-1]=wsum*weight1[1-1];
        } else if (cff < Real(-1.0)) {
            nfast=nfast-1;
            for (int i=1;i<=nfast;i++) {
                weight1[i-1]=weight1[i+1-1];
            }
            weight1[nfast+1-1] = zero;
        } else if (cff < zero) {
            wsum=one+cff;
            for (int i=1;i<=nfast-1;i++) {
                weight1[i-1]=wsum*weight1[i-1]-cff*weight1[i+1-1];
            }
            weight1[nfast-1]=wsum*weight1[nfast-1];
        }
    }

    //  Set SECONDARY weights assuming that backward Euler time step is used
    //  for free surface.  Notice that array weight2[i] is assumed to
    //  have all-zero status at entry in this segment of code.
    for(int j=1;j<=nfast;j++) {
        cff=weight1[j-1];
        for(int i=1;i<=j;i++) {
            weight2[i-1]=weight2[i-1]+cff;
        }
    }

    //
    //  Normalize both set of weights.
    //
    wsum = zero;
    cff  = zero;
    for(int i=1;i<=nfast;i++) {
        wsum=wsum+weight1[i-1];
        cff=cff+weight2[i-1];
    }

    wsum = one / wsum;
    cff  = one / cff;

    for(int i=1;i<=nfast;i++) {
        weight1[i-1]=wsum*weight1[i-1];
        weight2[i-1]=cff*weight2[i-1];
    }

    // advance_2d derives the barotropic time indices krhs/kstp from closed
    // forms in my_iif and iic, replacing ROMS's stateful indx1/next_indx1
    // toggling (main3d.F). Those closed forms reproduce ROMS only when nfast
    // is odd, since that is what makes indx1 flip parity exactly once per
    // baroclinic step. With an even nfast they desynchronise after the first
    // step and every subsequent substep reads the wrong time level, silently.
    // NDTFAST=44 gives nfast=61, but nothing else enforces this.
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nfast % 2 == 1,
        "barotropic time-index closed forms in advance_2d require an odd nfast; "
        "this ndtfast produces an even one and would desynchronise from ROMS");
}

// ---------------------------------------------------------------------------
// On the accumulator precision, because it looks like a defect and is not.
//
// set_weights.F declares wsum/shift/cff as `real(r16)` and mod_kinds.F labels
// r16 "128-bit", which reads as quad against the `Real` (double) used here. It
// is not quad on this build. mod_kinds.F only selects a 24-digit kind for
// SUN/AIX/NEC/SGI/CRAYX1/DEC; every other platform, Linux + gfortran included,
// gets
//     integer, parameter :: r16 = SELECTED_REAL_KIND(15,300)
// and 15 digits of precision with range 300 is satisfied by kind 8. Checked by
// compiling ROMS's own declaration and printing the kind: dp=8, r16=8,
// storage_size = 64 bits. ROMS accumulates these weights in DOUBLE, so double
// here is not a shortfall and `long double` would move AWAY from ROMS rather
// than towards it.
//
// What DID differ was the association order of the `scale` update above, fixed
// in place. Verified by transcribing set_weights.F verbatim into a standalone
// Fortran program and this routine into a standalone C++ one, compiling both
// under matched IEEE flags and comparing the final 61 weight pairs as raw bit
// patterns: 116/122 fields differed before, 0/122 after.
//
// One caveat that belongs with any bit-level claim about these weights. The
// production ROMS binary is built `-O3 -ffast-math` (Compilers/Linux-gfortran.mk
// line 68), and -ffast-math alone moves ROMS's OWN weights by up to 2.1e-14
// relative against the same source built with -ffp-contract=off. That is two
// orders larger than the round-off this routine can control, it is a build
// choice rather than an algorithm difference, and it is exactly the class of
// difference the campaign's ARM_KIND=compiler yardstick exists to measure. So
// the claim made here is the defensible one: this routine now matches ROMS's
// set_weights bit for bit when the two are compiled with equivalent semantics.
// ---------------------------------------------------------------------------
