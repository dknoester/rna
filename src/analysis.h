/* analysis.h
 *
 * This file is part of Self-Assembly.
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
#ifndef _ANALYSIS_H_
#define _ANALYSIS_H_

#include <boost/lexical_cast.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/sum.hpp>

#include <ea/algorithm.h>
#include <ea/analysis.h>
#include <ea/analysis/dominant.h>
#include "rna_record.h"
#include "rna.h"

//! Callback used to record an entry in the confusion matrix.
template <typename FitnessFunction>
struct confusion_callback : public FitnessFunction::callback {
    confusion_callback() : f(0.0), fl(0.0) {
        _C = scalar_matrix_type(5, 6, 0); // 6 + 9 bits for sequence+structure
    }

    virtual void classification(rna_record& r, matrix_type& s) {
        for(std::size_t i=0; i<r.output_size2(); ++i) { // columns
            int lbl = r.NL[i];
            bool none=true;
            for(std::size_t j=0; j<r.output_size1(); ++j, ++fl) { // rows
                if(s(j,i)) {
                    ++_C(lbl,j);
                    none=false;
                }
                f += (s(j,i) == r.N(j,i));
            }
            if(none) {
                ++_C(lbl,5);
            }
        }
    }

    matrix_type _C;
    double f,fl;
};

LIBEA_ANALYSIS_TOOL(rna_confusion) {
    typename EA::iterator i=analysis::dominant(ea);
    
    datafile df("rna_confusion.dat");
    // M == 5x5 matrix
    //  label   classification
    //         0  1  2  3  4  5
    //    0    A  B  C  ...
    //    1
    //    ...
    //
    // A is true positives; correct classification
    // B is 0s that were misclassified as 1, etc.
    
    confusion_callback<typename EA::fitness_function_type> cc;
    ea.fitness_function().reset_callback(&cc);
    recalculate_fitness(*i,ea);
    
    for(std::size_t i=0; i<cc._C.size1(); ++i) {
        df.write(i);
        for(std::size_t j=0; j<cc._C.size2(); ++j) {
            df.write(cc._C(i,j));
        }
        df.endl();
    }
    df.write(cc.f).write(cc.fl).write(cc.f/cc.fl).endl();
}


LIBEA_ANALYSIS_TOOL(rna_multiobjective_detail) {
    datafile df("rna_multiobjective_detail.dat");
    df.add_field("individual")
    .add_field("objective")
    .add_field("fitness");

    for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
        for(std::size_t j=0; j<ea.fitness_function().size(); ++j) {
            df.write(get<IND_NAME>(*i)).write(j).write(static_cast<double>(ealib::fitness(*i,ea)[j])).endl();
        }
    }
}


LIBEA_ANALYSIS_TOOL(rna_test) {
    // build up a classifier where each index i is the set of most-fit individuals on
    // objective i:
    std::vector<typename EA::population_type> classifiers(ea.fitness_function().size());
    
    for(std::size_t i=0; i<ea.fitness_function().size(); ++i) {
        std::sort(ea.population().begin(), ea.population().end(), comparators::objective<EA>(i,ea));
        // take all of the most fit on objective i, stuff them in classifiers[i].
        typename EA::population_type::reverse_iterator f=ea.population().rbegin();
        typename EA::population_type::reverse_iterator n=f+1;
        classifiers[i].push_back(*f);
        while((n!=ea.population().rend()) && (ealib::fitness(**f,ea)[i]==ealib::fitness(**n,ea)[i])) {
            classifiers[i].push_back(*n);
            ++f; ++n;
        }
    }
    
    // eval the performance of each population of individuals on their objective:
    datafile df("rna_test.dat");
    df.add_field("objective")
    .add_field("mean_p")
    .add_field("mean_tp")
    .add_field("mean_n")
    .add_field("mean_tn")
    .add_field("mean_mcc")
    .add_field("mean_acc")
    .add_field("mean_err");

    
    for(std::size_t i=0; i<classifiers.size(); ++i) {
        using namespace boost::accumulators;
        accumulator_set<double, stats<tag::mean> > mcc, acc, err;
        accumulator_set<double, stats<tag::mean> > p, tp, n, tn;
        for(typename EA::population_type::iterator j=classifiers[i].begin(); j!=classifiers[i].end(); ++j) {
            p((*j)->traits().test[i].p());
            tp((*j)->traits().test[i].tp());
            n((*j)->traits().test[i].n());
            tn((*j)->traits().test[i].tn());
            mcc((*j)->traits().test[i].mcc());
            acc((*j)->traits().test[i].acc());
            err((*j)->traits().test[i].err());
        }
        
        df.write(i).write(mean(p)).write(mean(tp)).write(mean(n)).write(mean(tn)).write(mean(mcc)).write(mean(acc)).write(mean(err)).endl();
    }
}

LIBEA_ANALYSIS_TOOL(rna_train) {
    // build up a classifier where each index i is the set of most-fit individuals on
    // objective i:
    std::vector<typename EA::population_type> classifiers(ea.fitness_function().size());
    
    for(std::size_t i=0; i<ea.fitness_function().size(); ++i) {
        std::sort(ea.population().begin(), ea.population().end(), comparators::objective<EA>(i,ea));
        // take all of the most fit on objective i, stuff them in classifiers[i].
        typename EA::population_type::reverse_iterator f=ea.population().rbegin();
        typename EA::population_type::reverse_iterator n=f+1;
        classifiers[i].push_back(*f);
        while((n!=ea.population().rend()) && (ealib::fitness(**f,ea)[i]==ealib::fitness(**n,ea)[i])) {
            classifiers[i].push_back(*n);
            ++f; ++n;
        }
    }
    
    // eval the performance of each population of individuals on their objective:
    datafile df("rna_train.dat");
    df.add_field("objective")
    .add_field("mean_p")
    .add_field("mean_tp")
    .add_field("mean_n")
    .add_field("mean_tn")
    .add_field("mean_mcc")
    .add_field("mean_acc")
    .add_field("mean_err");
    
    
    for(std::size_t i=0; i<classifiers.size(); ++i) {
        using namespace boost::accumulators;
        accumulator_set<double, stats<tag::mean> > mcc, acc, err;
        accumulator_set<double, stats<tag::mean> > p, tp, n, tn;
        for(typename EA::population_type::iterator j=classifiers[i].begin(); j!=classifiers[i].end(); ++j) {
            p((*j)->traits().train[i].p());
            tp((*j)->traits().train[i].tp());
            n((*j)->traits().train[i].n());
            tn((*j)->traits().train[i].tn());
            mcc((*j)->traits().train[i].mcc());
            acc((*j)->traits().train[i].acc());
            err((*j)->traits().train[i].err());
        }
        
        df.write(i).write(mean(p)).write(mean(tp)).write(mean(n)).write(mean(tn)).write(mean(mcc)).write(mean(acc)).write(mean(err)).endl();
    }
}


LIBEA_ANALYSIS_TOOL(rna_multi_confusion) {
    for(std::size_t k=0; k<ea.fitness_function().size(); ++k) {
        std::sort(ea.population().begin(), ea.population().end(), comparators::objective<EA>(k,ea));
        typename EA::reverse_iterator ind=ea.rbegin();
        
        std::ostringstream fname;
        fname << "rna_confusion_ojb" << k << ".dat";
        datafile df(fname.str());

        // M == 5x5 matrix
        //  label   classification
        //         0  1  2  3  4  5
        //    0    A  B  C  ...
        //    1
        //    ...
        //
        // A is true positives; correct classification
        // B is 0s that were misclassified as 1, etc.
        
        confusion_callback<typename EA::fitness_function_type> cc;
        ea.fitness_function().reset_callback(&cc);
        recalculate_fitness(*ind,ea);
        
        for(std::size_t i=0; i<cc._C.size1(); ++i) {
            df.write(i);
            for(std::size_t j=0; j<cc._C.size2(); ++j) {
                df.write(cc._C(i,j));
            }
            df.endl();
        }
        df.write(cc.f).write(cc.fl).write(cc.f/cc.fl).endl();
    }
}


#endif
