#ifndef _REMORA_TIME_UNITS_H_
#define _REMORA_TIME_UNITS_H_

#include <string>
#include <cstdlib>
#include <AMReX_REAL.H>
#include <AMReX_Print.H>

/**
 * \file REMORA_TimeUnits.H
 *
 * CF time-coordinate handling for NetCDF input series.
 *
 * REMORA historically required every time coordinate to be in "days" and
 * simply multiplied by 86400, ignoring the epoch in the units string. Real
 * ROMS input sets routinely mix both: the Moana nz5km files use
 *   ini    seconds since 2010-01-01
 *   meteo  days    since 2010-01-01
 *   bnd    seconds since 1950-01-01
 *   clim   days    since 1950-01-01
 *   rivers days    since 2010-01-01
 * Ignoring the epoch silently offsets one file from another by decades while
 * passing every existing check, so this converts each series to seconds since
 * a single model reference epoch.
 */
namespace remora_time {

/** Model reference epoch in seconds since 1970-01-01, set once from
 *  remora.time_ref during input parsing and read by every NetCDF time-series
 *  reader. Accessors rather than plain globals so this stays header-only. */
inline amrex::Real& model_ref_epoch () { static amrex::Real e = amrex::Real(0.0); return e; }
inline bool&        model_have_ref  () { static bool h = false; return h; }

/** Days from 1970-01-01 to the given civil date (proleptic Gregorian).
 *  Howard Hinnant's days_from_civil; exact for all dates in range. */
inline long long days_from_civil (long long y, unsigned m, unsigned d)
{
    y -= m <= 2;
    const long long era = (y >= 0 ? y : y - 399) / 400;
    const unsigned yoe = static_cast<unsigned>(y - era * 400);            // [0, 399]
    const unsigned doy = (153u * (m + (m > 2 ? -3 : 9)) + 2u) / 5u + d - 1; // [0, 365]
    const unsigned doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;           // [0, 146096]
    return era * 146097LL + static_cast<long long>(doe) - 719468LL;
}

/** Parse a CF units string of the form "<unit> since <YYYY-MM-DD[ HH:MM:SS]>".
 *
 * @param[in]  units      the units attribute
 * @param[out] scale      multiplier converting the stored values to seconds
 * @param[out] epoch_secs the epoch, in seconds since 1970-01-01T00:00:00
 * @return false if the string could not be parsed
 */
inline bool parse_cf_time_units (const std::string& units,
                                 amrex::Real& scale,
                                 amrex::Real& epoch_secs)
{
    const std::string u = units;
    if      (u.find("seconds") != std::string::npos) scale = amrex::Real(1.0);
    else if (u.find("minutes") != std::string::npos) scale = amrex::Real(60.0);
    else if (u.find("hours")   != std::string::npos) scale = amrex::Real(3600.0);
    else if (u.find("days")    != std::string::npos) scale = amrex::Real(86400.0);
    else return false;

    const std::size_t pos = u.find("since");
    if (pos == std::string::npos) return false;

    // Expect "... since YYYY-MM-DD[ HH:MM:SS]"
    std::string rest = u.substr(pos + 5);
    std::size_t b = rest.find_first_not_of(" \t");
    if (b == std::string::npos) return false;
    rest = rest.substr(b);

    long long Y = 0; unsigned M = 0, D = 0;
    int hh = 0, mm = 0; double ss = 0.0;
    // Date part is mandatory; the time part is optional.
    if (std::sscanf(rest.c_str(), "%lld-%u-%u", &Y, &M, &D) != 3) return false;
    const std::size_t tpos = rest.find_first_of(" T");
    if (tpos != std::string::npos) {
        std::sscanf(rest.c_str() + tpos + 1, "%d:%d:%lf", &hh, &mm, &ss);
    }
    if (M < 1 || M > 12 || D < 1 || D > 31) return false;

    const long long days = days_from_civil(Y, M, D);
    epoch_secs = static_cast<amrex::Real>(days) * amrex::Real(86400.0)
               + amrex::Real(hh) * amrex::Real(3600.0)
               + amrex::Real(mm) * amrex::Real(60.0)
               + static_cast<amrex::Real>(ss);
    return true;
}

/**
 * Convert a raw time value from a NetCDF file to seconds since the model's
 * reference epoch.
 *
 * @param raw         value as stored in the file
 * @param units       the variable's units attribute
 * @param ref_epoch   model reference epoch, seconds since 1970-01-01
 * @param have_ref    whether a model reference epoch was configured; if not,
 *                    the legacy behaviour (scale only, ignore epoch) is used
 *                    so existing single-epoch cases are unchanged
 * @param ok          set false if the units string could not be parsed
 */
inline amrex::Real to_model_seconds (amrex::Real raw,
                                     const std::string& units,
                                     amrex::Real ref_epoch,
                                     bool have_ref,
                                     bool& ok)
{
    amrex::Real scale = amrex::Real(86400.0);
    amrex::Real epoch = amrex::Real(0.0);
    ok = parse_cf_time_units(units, scale, epoch);
    if (!ok) {
        return raw * amrex::Real(86400.0);   // legacy assumption
    }
    if (!have_ref) {
        return raw * scale;                  // legacy: units respected, epoch ignored
    }
    return raw * scale + (epoch - ref_epoch);
}

} // namespace remora_time

#endif
