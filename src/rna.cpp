/* rna.cpp
 *
 * This file is part of RNA.
 *
 * Copyright 2014 David B. Knoester.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <boost/algorithm/string/predicate.hpp>
#include <mkv/markov_network_evolution.h>
#include <ea/generational_models/moran_process.h>
#include <ea/selection/rank.h>
#include <ea/selection/random.h>
#include <ea/datafiles/fitness.h>
#include <ea/cmdline_interface.h>
#include <mkv/analysis.h>
using namespace ealib;

#include "rna.h"

// Evolutionary algorithm definition.
typedef mkv::markov_network_evolution
< rna_fitness
, recombination::asexual
, generational_models::moran_process<selection::proportionate< >, selection::rank< > >
> ea_type;


/*! Define the EA's command-line interface.
 */
template <typename EA>
class cli : public cmdline_interface<EA> {
public:
    virtual void gather_options() {
        mkv::add_options(this);
        
        add_option<MORAN_REPLACEMENT_RATE_P>(this);
        add_option<POPULATION_SIZE>(this);
        add_option<RUN_UPDATES>(this);
        add_option<RUN_EPOCHS>(this);
        add_option<CHECKPOINT_PREFIX>(this);
        add_option<RNG_SEED>(this);
        add_option<RECORDING_PERIOD>(this);
        
        add_option<RNA_INPUT>(this);
        add_option<CA_RADIUS>(this);
    }
    
    virtual void gather_tools() {
        add_tool<mkv::dominant_reduced_graph>(this);
        add_tool<rna_confusion>(this);
    }
    
    virtual void gather_events(EA& ea) {
        add_event<datafiles::fitness_dat>(ea);
        //        add_event<rna_dat>(ea);
    }
    
    //! Called before initialization (good place to calculate config options).
    virtual void before_initialization(EA& ea) {
        int nin=(get<CA_RADIUS>(ea)*2+1)*20;
        put<mkv::MKV_INPUT_N>(nin, ea);
        put<mkv::MKV_OUTPUT_N>(10, ea);
    }
};
LIBEA_CMDLINE_INSTANCE(ea_type, cli);
