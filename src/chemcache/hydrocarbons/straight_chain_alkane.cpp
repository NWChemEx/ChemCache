/*
 * Copyright 2024 GhostFragment
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "hydrocarbons.hpp"
#include <simde/simde.hpp>

namespace chemcache::hydrocarbons {
namespace {

std::array<float, 3> position_carbon(const std::vector<float>& source_coords,
                                     float carbon_bond, int num,
                                     float angle_deg) {
    // Creates the pointer for positions
    std::array<float, 3> coords{0.0, 0.0, 0.0};
    // Convert degrees to radians
    float angle_rad = angle_deg * (3.14159265358979 / 180);

    // Sets the distance in each axis from the source coordinates
    // X: r*sin(theta)
    // Y: r*cos(theta)
    // Z: 0
    coords[0] = carbon_bond * sin(angle_rad / 2);
    coords[1] = carbon_bond * cos(angle_rad / 2);
    coords[2] = 0;

    // If the source atom is even, the new one must go up, otherwise,
    // it must go down
    coords[1] *= (num % 2 == 0) ? 1 : -1;

    for(int i = 0; i < 3; i++) { coords[i] += source_coords[(3 * num) + i]; }

    return coords;
}

std::array<float, 3> position_hydrogen(const std::vector<float>& source_coords,
                                       int flag, int num, float hydrogen_bond,
                                       float angle_deg) {
    // Creates the pointer for positions
    std::array<float, 3> coords{0.0, 0.0, 0.0};
    // Convert degrees to radians
    float angle_rad = angle_deg * (3.14159 / 180);

    // Because multiple hydrogen atoms get the same carbon atom, they need to be
    // defined by which place they are taking
    switch(flag) {
        // Goes behind the carbon (z is negative)
        case 0:
            // X: 0
            // Y: r*cos(theta)
            // Z: r*sin(theta)
            coords[0] = 0;
            coords[1] = hydrogen_bond * cos(angle_rad / 2);
            coords[2] = (-1) * hydrogen_bond * sin(angle_rad / 2);

            // If the carbon atom is odd, the hydrogen must be above, otherwise,
            // it must be below
            coords[1] *= (num % 2 == 0) ? -1 : 1;
            break;

        // Goes in front of the carbon (z is positive)
        case 1:
            // X: 0
            // Y: r*cos(theta)
            // Z: r*sin(theta)
            coords[0] = 0;
            coords[1] = hydrogen_bond * cos(angle_rad / 2);
            coords[2] = hydrogen_bond * sin(angle_rad / 2);

            // If the carbon atom is odd, the hydrogen must be above, otherwise,
            // it must be below
            coords[1] *= (num % 2 == 0) ? -1 : 1;
            break;

        // Goes to the left of the first carbon
        case 2:
            // X: r*sin(theta)
            // Y: r*cos(theta)
            // Z: 0
            coords[0] = (-1) * hydrogen_bond * sin(angle_rad / 2);
            coords[1] = hydrogen_bond * cos(angle_rad / 2);
            coords[2] = 0;
            break;

        // Goes to the right of the last carbon
        case 3:
            // X: r*sin(theta)
            // Y: r*cos(theta)
            // Z: 0
            coords[0] = hydrogen_bond * sin(angle_rad / 2);
            coords[1] = hydrogen_bond * cos(angle_rad / 2);
            coords[2] = 0;

            // If the last carbon is even, the hydrogen must be above,
            // otherwise, it must be below
            coords[1] *= (num % 2 == 0) ? 1 : -1;
            break;
    }

    for(int i = 0; i < 3; i++) { coords[i] += source_coords[(3 * num) + i]; }

    return coords;
}

#define C_C_BOND 2.89
#define H_C_BOND 2.06
#define ANGLE 109.5

#define C_MASS 21874.662
#define H_MASS 1837.289
#define C_NPROTON 6
#define H_NPROTON 1

chemist::Molecule straight_chain_alkane(int num_carbon) {
    chemist::Molecule m;

    // A log of all carbon atom positions.
    // Used to determine the position of the next carbon and connected hydrogens
    std::vector<float> source_coords{0.0, 0.0, 0.0};
    source_coords.reserve(num_carbon * 3);

    // The first carbon atom positioned at {0, 0, 0}
    m.push_back(chemist::Atom("C", C_NPROTON, C_MASS, 0, 0, 0));

    // Positions each of the following carbon atoms
    for(int i = 1; i < num_carbon; i++) {
        // Finds the new coordinates using the previous ones
        std::array<float, 3> coords =
          position_carbon(source_coords, C_C_BOND, i - 1, ANGLE);

        // Adds the new atom
        m.push_back(chemist::Atom("C", C_NPROTON, C_MASS, coords[0], coords[1],
                                  coords[2]));

        // Logs the new atoms coordinates
        for(int j = 0; j < 3; j++) { source_coords[(3 * i) + j] = coords[j]; }
    }

    // Positions each connected hydrogen atom
    for(int i = 0; i < num_carbon; i++) {
        // Creates the first hydrogen in the chain
        if(i == 0) {
            // Finds the coordinates using the carbon
            std::array<float, 3> coords =
              position_hydrogen(source_coords, 2, i, H_C_BOND, ANGLE);

            // Adds the new atom
            m.push_back(chemist::Atom("H", H_NPROTON, H_MASS, coords[0],
                                      coords[1], coords[2]));
        }

        // Creates the side hydrogens
        for(int j = 0; j <= 1; j++) {
            // Finds the coordinates using the carbon
            std::array<float, 3> coords =
              position_hydrogen(source_coords, j, i, H_C_BOND, ANGLE);

            // Adds the new atom
            m.push_back(chemist::Atom("H", H_NPROTON, H_MASS, coords[0],
                                      coords[1], coords[2]));
        }

        if(i == num_carbon - 1) {
            // Finds the coordinates using the carbon
            std::array<float, 3> coords =
              position_hydrogen(source_coords, 3, i, H_C_BOND, ANGLE);

            // Adds the new atom
            m.push_back(chemist::Atom("H", H_NPROTON, H_MASS, coords[0],
                                      coords[1], coords[2]));
        }
    }

    return m;
}

} // namespace

const auto mod_desc = R"(
Straight-Chain Alkane Generator
-------------------------------

This module performs basic geometry to compute the positions of the carbons and
hydrogens in an arbitrary length straight-chain alkane. The resulting alkane
will have carbon-carbon bonds determined by ``c_c_bond`` parameter and carbon-
hydrogen bonds determined by the ``h_c_bond`` parameter. All angles will be set
to 109.5 (idealized tetrahedron).

TODO: Expose the bond length and angle parameters.
TODO: Use non-idealized angles.
TODO: Make masses parameters.
)";

MODULE_CTOR(StraightChainAlkane) {
    description(mod_desc);
    satisfies_property_type<simde::MoleculeFromString>();
}

MODULE_RUN(StraightChainAlkane) {
    const auto& [num_carbon_string] =
      simde::MoleculeFromString::unwrap_inputs(inputs);

    auto hc = straight_chain_alkane(std::stoi(num_carbon_string));

    auto rv = results();

    return simde::MoleculeFromString::wrap_results(rv, hc);
}

} // namespace chemcache::hydrocarbons
