// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "WellConfigurationFile.h"
#include "Exception.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>

using namespace std;

namespace csmp {

namespace {

/// Everything before a '#', trimmed at both ends.
string StripComment(const string &line) {
    string s = line.substr(0, line.find('#'));
    const size_t b = s.find_first_not_of(" \t\r\n");
    if (b == string::npos) return "";
    const size_t e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
}

string ToLower(string s) {
    transform(s.begin(), s.end(), s.begin(),
              [](unsigned char c) { return static_cast<char>(tolower(c)); });
    return s;
}

/// Splits on whitespace and tabs.
vector<string> Tokenise(const string &s) {
    vector<string> out;
    istringstream iss(s);
    string t;
    while (iss >> t) out.push_back(t);
    return out;
}

[[noreturn]] void Fail(const string &path, size_t line_no, const string &what) {
    ostringstream msg;
    msg << path << ":" << line_no << ": " << what;
    throw csmp::Exception(ERROR, "WellConfigurationFile", msg.str());
}

double ToDouble(const string &path, size_t line_no, const string &token) {
    try {
        size_t used = 0;
        const double v = stod(token, &used);
        if (used != token.size())
            Fail(path, line_no, "'" + token + "' is not a number");
        return v;
    }
    catch (const csmp::Exception &) { throw; }
    catch (...) {
        Fail(path, line_no, "'" + token + "' is not a number");
    }
}

/** Keys are written with spaces in the file ("radius depths") for readability;
    match them case-insensitively and on the joined leading tokens. Returns the
    number of tokens the key consumed, or 0 if it does not match. */
size_t MatchKey(const vector<string> &tok, const string &key) {
    const vector<string> kt = Tokenise(key);
    if (tok.size() < kt.size()) return 0;
    for (size_t i = 0; i < kt.size(); ++i)
        if (ToLower(tok[i]) != kt[i]) return 0;
    return kt.size();
}

bool ToBool(const string &path, size_t line_no, const string &token) {
    const string t = ToLower(token);
    if (t == "true"  || t == "yes" || t == "on"  || t == "1") return true;
    if (t == "false" || t == "no"  || t == "off" || t == "0") return false;
    Fail(path, line_no, "'" + token + "' is not a boolean (use true/false)");
}

} // anonymous namespace

map<string, WellConfiguration> WellConfigurationFile::ReadAll(const string &path) {

    ifstream ifs(path);
    if (!ifs)
        throw csmp::Exception(ERROR, "WellConfigurationFile::ReadAll",
                              path, "cannot be opened for reading.");

    map<string, WellConfiguration> configs;
    string current_well;
    string line;
    size_t line_no = 0;

    while (getline(ifs, line)) {
        ++line_no;
        const string stripped = StripComment(line);
        if (stripped.empty()) continue;

        vector<string> tok = Tokenise(stripped);

        // ── section header ──────────────────────────────────────────────────
        if (ToLower(tok[0]) == "[well]") {
            if (tok.size() != 2)
                Fail(path, line_no, "[well] must be followed by exactly one well name");
            current_well = tok[1];
            if (configs.count(current_well))
                Fail(path, line_no, "well '" + current_well + "' is declared twice");
            configs[current_well];            // default-construct
            continue;
        }

        if (current_well.empty())
            Fail(path, line_no, "entry outside any [well] section");

        WellConfiguration &c = configs[current_well];

        // ── geometry ────────────────────────────────────────────────────────
        if (size_t k = MatchKey(tok, "radius depths")) {
            c.radius_depths.clear();
            for (size_t i = k; i < tok.size(); ++i)
                c.radius_depths.push_back(ToDouble(path, line_no, tok[i]));
            continue;
        }
        if (size_t k = MatchKey(tok, "radius values")) {
            c.radius_values.clear();
            for (size_t i = k; i < tok.size(); ++i)
                c.radius_values.push_back(ToDouble(path, line_no, tok[i]));
            if (c.radius_values.empty())
                Fail(path, line_no, "radius values needs at least one radius");
            continue;
        }
        if (size_t k = MatchKey(tok, "virtual top radius")) {
            c.virtual_top_radius = ToDouble(path, line_no, tok[k]);
            continue;
        }
        if (size_t k = MatchKey(tok, "pipe roughness")) {
            c.pipe_roughness = ToDouble(path, line_no, tok[k]);
            continue;
        }
        if (size_t k = MatchKey(tok, "default skin")) {
            c.default_skin = ToDouble(path, line_no, tok[k]);
            continue;
        }
        if (size_t k = MatchKey(tok, "angle fraction")) {
            c.angle_fraction = ToDouble(path, line_no, tok[k]);
            continue;
        }
        if (size_t k = MatchKey(tok, "init bulk smf")) {
            c.init_bulk_smf = ToDouble(path, line_no, tok[k]);
            continue;
        }

        // ── completions ─────────────────────────────────────────────────────
        if (size_t k = MatchKey(tok, "completion length")) {
            const string rule = ToLower(tok[k]) + (tok.size() > k + 1 ? " " + ToLower(tok[k + 1]) : "");
            if (rule == "element lengths")
                c.completion_length_rule = CompletionLengthRule::FromElementLengths;
            else if (rule == "host thickness")
                c.completion_length_rule = CompletionLengthRule::FromHostThickness;
            else
                Fail(path, line_no, "completion length must be 'element lengths' or "
                                    "'host thickness', not '" + rule + "'");
            continue;
        }
        if (size_t k = MatchKey(tok, "completion")) {
            const size_t n = tok.size() - k;
            if (n < 2 || n > 3)
                Fail(path, line_no, "completion needs a top and a bottom depth, "
                                    "optionally followed by a skin");
            CompletionInterval ci;
            ci.top_depth    = ToDouble(path, line_no, tok[k]);
            ci.bottom_depth = ToDouble(path, line_no, tok[k + 1]);
            if (n == 3) ci.skin = ToDouble(path, line_no, tok[k + 2]);
            c.completions.push_back(ci);
            continue;
        }

        // ── discretisation ──────────────────────────────────────────────────
        if (size_t k = MatchKey(tok, "target segment length")) {
            c.target_segment_length = ToDouble(path, line_no, tok[k]);
            continue;
        }
        if (size_t k = MatchKey(tok, "virtual top segments")) {
            const double v = ToDouble(path, line_no, tok[k]);
            if (v < 0.) Fail(path, line_no, "virtual top segments cannot be negative");
            c.number_of_virtual_top_segments = static_cast<unsigned>(v);
            continue;
        }

        // ── initial state ───────────────────────────────────────────────────
        if (size_t k = MatchKey(tok, "initial temperature")) {
            const string mode = ToLower(tok[k]) + (tok.size() > k + 1 ? " " + ToLower(tok[k + 1]) : "");
            if (mode == "follow reservoir") {
                c.initial_temperature.mode = InitialWellTemperatureMode::FollowReservoir;
            }
            else if (ToLower(tok[k]) == "uniform") {
                if (tok.size() < k + 2) Fail(path, line_no, "uniform needs a temperature");
                c.initial_temperature.mode  = InitialWellTemperatureMode::Uniform;
                c.initial_temperature.value = ToDouble(path, line_no, tok[k + 1]);
            }
            else if (mode == "reservoir minimum") {
                if (tok.size() < k + 3) Fail(path, line_no, "reservoir minimum needs a temperature");
                c.initial_temperature.mode  = InitialWellTemperatureMode::ReservoirWithMinimum;
                c.initial_temperature.value = ToDouble(path, line_no, tok[k + 2]);
            }
            else
                Fail(path, line_no, "initial temperature must be 'follow reservoir', "
                                    "'uniform <degC>' or 'reservoir minimum <degC>'");
            continue;
        }
        if (size_t k = MatchKey(tok, "top pressure")) {
            c.top_pressure = ToDouble(path, line_no, tok[k]);
            continue;
        }
        if (size_t k = MatchKey(tok, "top temperature")) {
            c.top_temperature = ToDouble(path, line_no, tok[k]);
            continue;
        }

        if (size_t k = MatchKey(tok, "radial heat model")) {
            const string m = ToLower(tok[k]);
            if      (m == "shapefactor" || m == "shape") c.radial_heat_model = RadialHeatModel::ShapeFactor;
            else if (m == "ramey")                        c.radial_heat_model = RadialHeatModel::Ramey;
            else Fail(path, line_no, "radial heat model must be 'shape factor' or 'ramey'");
            continue;
        }
        if (size_t k = MatchKey(tok, "formation thermal conductivity")) {
            c.formation_thermal_conductivity = ToDouble(path, line_no, tok[k]);
            continue;
        }
        if (size_t k = MatchKey(tok, "formation thermal diffusivity")) {
            c.formation_thermal_diffusivity = ToDouble(path, line_no, tok[k]);
            continue;
        }

        // ── physics ─────────────────────────────────────────────────────────
        if (size_t k = MatchKey(tok, "with friction")) {
            c.with_friction = ToBool(path, line_no, tok[k]); continue; }
        if (size_t k = MatchKey(tok, "with gravitational potential energy")) {
            c.with_grav_pot_energy = ToBool(path, line_no, tok[k]); continue; }
        if (size_t k = MatchKey(tok, "with kinetic energy")) {
            c.with_kinetic_energy = ToBool(path, line_no, tok[k]); continue; }
        if (size_t k = MatchKey(tok, "with inertial terms")) {
            c.with_inertial_terms_in_momentum = ToBool(path, line_no, tok[k]); continue; }
        if (size_t k = MatchKey(tok, "with radial heat")) {
            c.with_radial_heat = ToBool(path, line_no, tok[k]); continue; }
        if (size_t k = MatchKey(tok, "spread source")) {
            c.spread_source = ToBool(path, line_no, tok[k]); continue; }

        Fail(path, line_no, "unknown key '" + stripped + "'");
    }

    if (configs.empty())
        throw csmp::Exception(ERROR, "WellConfigurationFile::ReadAll",
                              path, "contains no [well] section.");

    // Validate each well here rather than leaving it to Initialize_Well, so a
    // bad file is reported before the model is built.
    for (const auto &kv : configs) {
        try { kv.second.Validate(); }
        catch (const csmp::Exception &e) {
            throw csmp::Exception(ERROR, "WellConfigurationFile::ReadAll",
                                  "well '" + kv.first + "' in " + path, e.What());
        }
    }

    cerr << endl << "Read " << configs.size() << " well configuration(s) from " << path << ":";
    for (const auto &kv : configs)
        cerr << endl << "  " << kv.first << ": " << kv.second.completions.size()
             << " completion interval(s), " << kv.second.radius_values.size()
             << " radius section(s)";

    return configs;
}

void WellConfigurationFile::WriteTemplate(const string &path) {

    ofstream out(path);
    if (!out)
        throw csmp::Exception(ERROR, "WellConfigurationFile::WriteTemplate",
                              path, "cannot be opened for writing.");

    const WellConfiguration d;   // defaults, so this file cannot drift from the code

    out <<
        "# ============================================================================\n"
        "# Well configuration\n"
        "#\n"
        "# One [well] section per well, named after its region in the mesh. Entries may\n"
        "# appear in any order; anything you leave out keeps its default, so a section\n"
        "# can be as short as a single 'completion' line.\n"
        "#\n"
        "# '#' starts a comment. Tokens are separated by spaces or tabs. Depths are\n"
        "# ABSOLUTE, in the same datum as the mesh, negative downward.\n"
        "#\n"
        "# NOT set here: wellhead pressure and temperature, target rates and injection\n"
        "# rates. Those are operating conditions that change during a run and are set\n"
        "# through the scheme's Set_well_* methods.\n"
        "# ============================================================================\n"
        "\n"
        "[well]	INJECTOR\n"
        "\n"
        "# ---- Geometry --------------------------------------------------------------\n"
        "# Radius as a function of depth, as two parallel lists. The first radius applies\n"
        "# from the well top down to the first depth, the next from there to the second,\n"
        "# and the LAST radius continues to the toe — so 'radius values' has exactly one\n"
        "# more entry than 'radius depths'. Depths must be strictly decreasing.\n"
        "# For a uniform well, omit 'radius depths' and give a single radius.\n"
        "#\n"
        "# The radius enters the well index only through log(re/rw), so it is a weak\n"
        "# parameter; it matters more for the flow equations through the cross-section.\n"
        "radius depths		-1100.	-1900.\n"
        "radius values		0.17	0.12	0.11\n"
        "\n"
        "# Pipe roughness [m], for the friction factor.\n"
        "pipe roughness		" << d.pipe_roughness << "\n"
                               "\n"
                               "# ---- Completions -----------------------------------------------------------\n"
                               "# 'completion  <top depth>  <bottom depth>  [skin]'\n"
                               "# Repeat the line for a multi-completion well. A node is perforated if it falls\n"
                               "# inside ANY interval; overlaps are allowed and the first match wins. ONLY\n"
                               "# perforated nodes exchange with the reservoir, so this list is the well's\n"
                               "# connection to the model, not a detail.\n"
                               "#\n"
                               "# The optional third value is the skin for that interval. Omit it and the\n"
                               "# interval uses 'default skin' below.\n"
                               "completion		-1900.	-10000.\n"
                               "\n"
                               "# Skin for this well. Applies to every perforated node except those in an\n"
                               "# interval that gave its own skin (the optional third value on its 'completion'\n"
                               "# line), so there is no need to set one per completion. It is a STARTING value:\n"
                               "# it seeds the 'well skin' node property, which the model re-reads every step, so\n"
                               "# the skin can be changed during a run by writing that property.\n"
                               "default skin		" << d.default_skin << "\n"
                             "\n"
                             "# How the completion length of each perforated node is obtained. It scales the\n"
                             "# well index LINEARLY, so it matters considerably more than the radius.\n"
                             "#   element lengths  half of each adjoining well element (default)\n"
                             "#   host thickness   per node, where the host reservoir element is lower-\n"
                             "#                    dimensional (a surface element, i.e. a fault mid-region),\n"
                             "#                    use that element's thickness — the aperture. Nodes whose\n"
                             "#                    host is a volume element fall back to element lengths.\n"
                             "completion length	element lengths\n"
                             "\n"
                             "# ---- Discretisation --------------------------------------------------------\n"
                             "# Well elements are subdivided into segments of about this length [m].\n"
                             "target segment length	" << d.target_segment_length << "\n"
                                      "\n"
                                      "# Segments modelled ABOVE the top of the mesh, for a domain that does not reach\n"
                                      "# the surface. UNTESTED: this has been 0 for a long time and nothing exercises\n"
                                      "# those branches. Leave at 0 unless you are prepared to verify the result.\n"
                                      "virtual top segments	" << d.number_of_virtual_top_segments << "\n"
                                               "\n"
                                               "# Radius of those virtual segments [m]; negative means 'use the topmost radius'.\n"
                                               "virtual top radius	" << d.virtual_top_radius << "\n"
                                   "\n"
                                   "# ---- Initial state ---------------------------------------------------------\n"
                                   "# How the well temperature is set up at t = 0. One of:\n"
                                   "#   follow reservoir            Twell = Treservoir\n"
                                   "#   uniform <degC>              Twell = value everywhere\n"
                                   "#   reservoir minimum <degC>    Twell = max(Treservoir, value)\n"
                                   "initial temperature	follow reservoir\n"
                                   "\n"
                                   "# Initial bulk SALT mass fraction in the well [kg NaCl / kg bulk] — one of the\n"
                                   "# model's primary variables, not a steam fraction. Set it to the salinity of the\n"
                                   "# fluid the well holds at t = 0.\n"
                                   "init bulk smf		" << d.init_bulk_smf << "\n"
                              "\n"
                              "# Starting values at the top of the well. The run-time wellhead conditions are\n"
                              "# set separately, through the scheme.\n"
                              "top pressure		" << d.top_pressure << "\n"
                             "top temperature		" << d.top_temperature << "\n"
                                "\n"
                                "# ---- Physics ---------------------------------------------------------------\n"
                                "# Each of these has a working alternative; the defaults are what every result to\n"
                                "# date was produced with. The comments give the usual magnitude of each term.\n"
                                "with friction				true	# wall friction\n"
                                "with gravitational potential energy	true	# generally important\n"
                                "with kinetic energy			true	# generally negligible\n"
                                "with inertial terms			true	# generally very negligible\n"
                                "with radial heat			true	# conduction to the formation\n"
                                "\n"
                                "# How that conduction is computed:\n"
                                "#   shape factor  steady shape factor to the first reservoir control volume, the\n"
                                "#                 thermal analogue of the Peaceman well index (default). Needs no\n"
                                "#                 time variable and tightens as the mesh is refined.\n"
                                "#   ramey         Ramey (1962) transient line source. For meshes that do NOT\n"
                                "#                 resolve the near-well formation. Long-time approximation: it\n"
                                "#                 returns zero below roughly alpha*t/r_w^2 = 1.\n"
                                "radial heat model	shape factor\n"
                                "\n"
                                "# Formation thermal conductivity [W/(m.K)] for the heat-loss model. Set it to\n"
                                "# match the reservoir's own value. It is not read from the reservoir because\n"
                                "# 'thermal conductivity' is an element property and the well model only has\n"
                                "# nodal values available.\n"
                                "formation thermal conductivity	2.0\n"
                                "\n"
                                "# Formation thermal diffusivity [m2/s], alpha = kappa/(rho*c). Used ONLY by the\n"
                                "# ramey model. ~9.3e-7 for kappa 2, rho 2500, c 860.\n"
                                "formation thermal diffusivity	9.3e-7\n"
                                "\n"
                                "# NOTE: water-table displacement is NOT set here. It is switched on together\n"
                                "# with the depth it starts from, through the scheme's Set_well_water_table(),\n"
                                "# which is also how it changes during a run. It defaults to off.\n"
                                "\n"
                                "# Advanced. Where a node's reservoir exchange is placed among the sub-segments\n"
                                "# of its well element. false: all of it on the last sub-segment (what every\n"
                                "# result to date used). true: spread evenly over the element's sub-segments,\n"
                                "# which represents a distributed completion better. Re-check global mass balance\n"
                                "# if you change it.\n"
                                "spread source				false\n"
                                "\n"
                                "# ============================================================================\n"
                                "[well]	PRODUCER\n"
                                "\n"
                                "radius depths		-1100.	-1900.\n"
                                "radius values		0.17	0.12	0.11\n"
                                "completion		-1900.	-10000.\n"
                                "default skin		0.\n"
                                "init bulk smf		0.\n"
                                "completion length	host thickness\n"
                                "initial temperature	follow reservoir\n";

    cerr << endl << "Wrote a commented well-configuration template to " << path;
}

} // namespace csmp