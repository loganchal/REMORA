#include <cstdio>
#include <string>
#include <cstdlib>
namespace amrex { using Real = double;
  inline void Print(){} }
#include "tu_body.h"
int main(){
    double dummy, ref2010, ref1950; bool ok;
    parse_cf_time_units("seconds since 2010-01-01", dummy, ref2010);
    parse_cf_time_units("seconds since 1950-01-01", dummy, ref1950);
    printf("epoch2010=%.0f expect 1262304000 %s\n", ref2010, ref2010==1262304000?"OK":"FAIL");
    printf("epoch1950=%.0f expect -631152000 %s\n", ref1950, ref1950==-631152000?"OK":"FAIL");
    double a=to_model_seconds(4018.0,"days since 2010-01-01",ref2010,true,ok);
    double b=to_model_seconds(25933.0,"days since 1950-01-01",ref2010,true,ok);
    double c=to_model_seconds(2240611200.0,"seconds since 1950-01-01",ref2010,true,ok);
    printf("meteo days/2010 -> %.0f %s\n", a, a==347155200?"OK":"FAIL");
    printf("clim  days/1950 -> %.0f %s (must equal meteo)\n", b, b==347155200?"OK":"FAIL");
    printf("bnd   secs/1950 -> %.0f %s (must equal meteo)\n", c, c==347155200?"OK":"FAIL");
    double l=to_model_seconds(4018.0,"days since 2010-01-01",0.0,false,ok);
    printf("legacy(no ref)  -> %.0f %s\n", l, l==347155200?"OK":"FAIL");
    double e; bool p=parse_cf_time_units("days since 2020-02-29 12:00:00", dummy, e);
    printf("leap parse ok=%d epoch=%.0f %s\n", p, e, (p&&e==1582977600)?"OK":"FAIL");
    double h=to_model_seconds(2.0,"hours since 2010-01-01",ref2010,true,ok);
    printf("hours unit      -> %.0f %s\n", h, h==7200?"OK":"FAIL");
    return 0;
}
// Build standalone (no AMReX needed -- the header's math is pure):
//   python3 - <<'PY'
//   src=open('Source/IO/REMORA_TimeUnits.H').read()
//   open('tu_body.h','w').write(src[src.index('inline long long days_from_civil'):src.index('} // namespace remora_time')].replace('namespace remora_time {',''))
//   PY
//   c++ -O0 -o tu Tests/MoanaParity/test_time_units.cpp && ./tu
//
// Verifies that the three time conventions in the Moana input set --
//   meteo  days    since 2010-01-01
//   clim   days    since 1950-01-01
//   bnd    seconds since 1950-01-01
// all resolve to the SAME model time (347155200 s = 2021-01-01 on a 2010
// reference). Before remora.time_ref existed they did not, and nothing
// complained.
