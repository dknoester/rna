/* rna_nsga2.cpp
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
#include <mkv/markov_network_evolution.h>
#include <ea/nsga2.h>
#include <ea/datafiles/multiobjective_fitness.h>
#include <ea/cmdline_interface.h>
#include <ea/math/roc.h>
#include <mkv/analysis.h>
using namespace ealib;

#include "rna_nsga2.h"

template <typename T>
struct rna_traits : nsga2_traits<T> {
    typedef nsga2_traits<T> parent;
    
    //! Constructor.
    rna_traits() {
    }
    
    //! Serialization.
    template <class Archive>
    void serialize(Archive& ar, const unsigned int version) {
        ar & boost::serialization::make_nvp("nsga2_trait", boost::serialization::base_object<parent>(*this));
        ar & boost::serialization::make_nvp("train", train);
        ar & boost::serialization::make_nvp("test", test);
    }
    
    std::vector<roc> train;
    std::vector<roc> test;
};

// Evolutionary algorithm definition.
typedef mkv::markov_network_evolution
< multi_rna_fitness
, recombination::two_point_crossover
, generational_models::nsga2
, dont_stop
, fill_population
, mkv::markov_network_lifecycle
, rna_traits
> ea_type;

template <typename EA>
class cli : public cmdline_interface<EA> {
public:
    virtual void gather_options() {
        mkv::add_options(this);
        
        add_option<POPULATION_SIZE>(this);
        add_option<RUN_UPDATES>(this);
        add_option<RUN_EPOCHS>(this);
        add_option<CHECKPOINT_PREFIX>(this);
        add_option<RNG_SEED>(this);
        add_option<RECORDING_PERIOD>(this);
        
        add_option<RNA_INPUT>(this);
        add_option<RNA_LIMIT>(this);
        add_option<RNA_CROSSVAL>(this);
        add_option<CA_RADIUS>(this);
    }
    
    virtual void gather_tools() {
        add_tool<rna_multiobjective_detail>(this);
        add_tool<rna_train>(this);
        add_tool<rna_test>(this);
        add_tool<rna_multi_confusion>(this);
        add_tool<mkv::multi_reduced_graph>(this);
    }
    
    virtual void gather_events(EA& ea) {
        add_event<datafiles::multiobjective_fitness_dat>(ea);
    }
    
    //! Called before initialization (good place to calculate config options).
    virtual void before_initialization(EA& ea) {
        int nin=(get<CA_RADIUS>(ea)*2+1)*20;
        put<mkv::MKV_INPUT_N>(nin, ea);
        put<mkv::MKV_OUTPUT_N>(10, ea);
    }
};
LIBEA_CMDLINE_INSTANCE(ea_type, cli);
