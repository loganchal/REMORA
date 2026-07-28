// Standalone oracle harness for the REMORA Vstretching=5 port.
// Carries the EXACT expressions from Source/Utils/REMORA_DepthStretchTransform.H
// (branch moana/vstretching5). Prints s, Cs, and Vtransform=2 depths at full
// precision for comparison against the atlas Python reference.
#include <cmath>
#include <cstdio>

using Real = double;
static const Real zero = 0.0, one = 1.0, two = 2.0;

int main() {
    const int  N = 50;
    const Real theta_s = 6.0, theta_b = 2.0, hc = 250.0;
    const int  vstretching = 5;

    Real s_w[N+1], Cs_w[N+1], s_r[N], Cs_r[N];
    const Real ds = one / Real(N);

    for (int k = 0; k <= N; ++k) {
        Real Csur, Cbot;
        if (vstretching == 5) {
            const Real rN = Real(N);
            const Real rk = Real(k);
            if (k == N) {
                s_w[k] = zero;
            } else if (k == 0) {
                s_w[k] = -one;
            } else {
                s_w[k] = -(rk*rk - two*rk*rN + rk + rN*rN - rN)/(rN*rN - rN)
                         - Real(0.01)*(rk*rk - rk*rN)/(one - rN);
            }
        } else {
            s_w[k] = ds*(k-N);
        }
        if (theta_s > zero) {
            Csur = (one - std::cosh(theta_s*s_w[k]))/(std::cosh(theta_s) - one);
        } else {
            Csur = -s_w[k]*s_w[k];
        }
        if (theta_b > zero) {
            Cbot = (std::exp(theta_b*Csur) - one)/(one - std::exp(-theta_b));
            Cs_w[k] = Cbot;
        } else {
            Cs_w[k] = Csur;
        }
        if (vstretching == 5 && k == 0) Cs_w[k] = -one;

        if (k < N) {
            if (vstretching == 5) {
                const Real rN = Real(N);
                const Real rk = Real(k) + Real(0.5);
                s_r[k] = -(rk*rk - two*rk*rN + rk + rN*rN - rN)/(rN*rN - rN)
                         - Real(0.01)*(rk*rk - rk*rN)/(one - rN);
            } else {
                s_r[k] = ds*(k-N+Real(0.5));
            }
            if (theta_s > zero) {
                Csur = (one - std::cosh(theta_s*s_r[k]))/(std::cosh(theta_s) - one);
            } else {
                Csur = -s_r[k]*s_r[k];
            }
            if (theta_b > zero) {
                Cbot = (std::exp(theta_b*Csur) - one)/(one - std::exp(-theta_b));
                Cs_r[k] = Cbot;
            } else {
                Cs_r[k] = Csur;
            }
        }
    }

    // coefficient tables
    for (int k = 0; k <= N; ++k)
        printf("W %d %.17g %.17g\n", k, s_w[k], Cs_w[k]);
    for (int k = 0; k < N; ++k)
        printf("R %d %.17g %.17g\n", k, s_r[k], Cs_r[k]);

    // Vtransform=2 depths, as in stretch_transform: cff2=(hc*s+Cs*h)/(hc+h);
    // z = zeta + (zeta+h)*cff2
    const Real cases[5][2] = {{50,0},{1000,0},{4000,0},{1000,1.2},{25,-0.8}};
    for (int c = 0; c < 5; ++c) {
        const Real h = cases[c][0], zeta = cases[c][1];
        for (int k = 0; k < N; ++k) {
            Real cff2 = (hc*s_r[k] + Cs_r[k]*h)/(hc + h);
            printf("ZR %d %g %g %.17g\n", k, h, zeta, zeta + (zeta+h)*cff2);
        }
        for (int k = 0; k <= N; ++k) {
            Real cff2 = (hc*s_w[k] + Cs_w[k]*h)/(hc + h);
            printf("ZW %d %g %g %.17g\n", k, h, zeta, zeta + (zeta+h)*cff2);
        }
    }
    return 0;
}
