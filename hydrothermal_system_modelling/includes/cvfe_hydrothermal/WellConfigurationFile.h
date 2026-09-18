// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CVFEM_WELL_CONFIGURATION_FILE_H
#define CVFEM_WELL_CONFIGURATION_FILE_H

/*  Reading WellConfiguration objects from a text file.

    Deliberately NOT part of InputDataManager: that reader's blocks are
    positional (read in order, separated by blank lines, with the `# (blockN)`
    lines being comments the parser never sees), and its existing well block is a
    flat map<string, vector<double>> — neither fits a structure with repeated
    completion intervals, optional fields and enums. This reader is keyed rather
    than positional, so entries can appear in any order and an absent key simply
    leaves the WellConfiguration default in place.

    Conventions follow the CSMP configuration files: '#' starts a comment,
    tokens are whitespace- or tab-separated, blank lines are ignored.

    The authoritative description of every option is WriteTemplate(), which emits
    a fully commented file. Keeping the documentation in the same place as the
    parser is the point — a template written by hand drifts from the code.

    Benoit DD/MM/YYYY
*/

#include "WellModelPrototype.h"

#include <map>
#include <string>

namespace csmp {

class WellConfigurationFile {
public:
    /** Reads every [well] section of `path`.

        @return one WellConfiguration per well name. Keys not present keep their
                WellConfiguration default, so a section may be as short as a
                single `completion` line.
        @throws csmp::Exception on a missing file, an unknown key, a malformed
                value or a value count that cannot be right — always naming the
                file and line, because a silently ignored typo in a config file
                is worse than no config file at all. */
    static std::map<std::string, WellConfiguration> ReadAll(const std::string &path);

    /** Writes a fully commented template with every key at its default value.
        This IS the documentation of the format. */
    static void WriteTemplate(const std::string &path);
};

} // namespace csmp

#endif /* CVFEM_WELL_CONFIGURATION_FILE_H */
