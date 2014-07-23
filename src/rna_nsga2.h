/* rna.h
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
#ifndef _RNA_NSGA2_H_
#define _RNA_NSGA2_H_

#include <ea/fitness_function.h>

#include "rna.h"

struct multi_rna_fitness : public fitness_function<multivalued_fitness<double>, constantS, stochasticS> {
    typedef std::vector<rna_record> rna_vector;
    
    rna_vector _rna; //!< List of all available RNA protein binding records.
    std::size_t _length; //!< Total length of all sequences.
    std::size_t _max_length; //!< Maximum length of a sequence
    
    struct callback {
        virtual void classification(rna_record& r, matrix_type& s) = 0;
    };
    
    callback* _cb;
    
    //! Reset the callback pointer.
    void reset_callback(callback* c=0) {
        _cb = c;
    }
    
    //! Initialize this fitness function.
    template <typename RNG, typename EA>
    void initialize(RNG& rng, EA& ea) {
        _cb = 0;
        idatafile df(get<RNA_INPUT>(ea), ",");
        df.translate(_rna,ea);
        
        _max_length = 0;
        _length = 0;
        for(rna_vector::iterator i=_rna.begin(); i!=_rna.end(); ++i) {
            _length += i->S.size();
            _max_length = std::max(i->S.size(), _max_length);
        }
        std::random_shuffle(_rna.begin(), _rna.end(), rng);
        
        if(exists<RNA_LIMIT>(ea)) {
            assert(get<RNA_LIMIT>(ea) < _rna.size());
            _rna.erase(_rna.begin()+get<RNA_LIMIT>(ea), _rna.end());
        }
    }
    
    //! Number of objectives associated with this fitness function.
    std::size_t size() const {
        return 5;
    }

    //! Returns the range of the m'th objective.
    double range(std::size_t m) {
        return 1.0;
    }
    
    //! Calculate fitness of an individual.
    template <typename Individual, typename RNG, typename EA>
    value_type operator()(Individual& ind, RNG& rng, EA& ea) {
        // state containers:
        typedef boost::numeric::ublas::matrix<int> matrix_type;
        const int nrow=5;
        matrix_type S1(nrow,_max_length), S2(nrow,_max_length);
        matrix_type* pt=&S1;
        matrix_type* ptp1=&S2;
        
        // ca:
        const int radius = get<CA_RADIUS>(ea);
        std::vector<typename EA::phenotype_type> ca(_max_length, ealib::phenotype(ind, ea));
        for(std::size_t i=0; i<ca.size(); ++i) {
            ca[i].reset(rng.seed());
        }
        
        ind.traits().train.resize(5);
        ind.traits().test.resize(5);
        
        // for each RNA record:
        for(std::size_t i=0; i<_rna.size(); ++i) {
            rna_record& r = _rna[i];
            int N = r.input_size2();
            
            // reset the ca:
            for(int j=0; j<N; ++j) {
                ca[j].clear();
            }
            S1 = scalar_matrix_type(nrow,_max_length,0);
            S2 = scalar_matrix_type(nrow,_max_length,0);
            
            /* each column is:
             15 bits (sequence + structure) WARNING: bug here; these are overlaid; only 0-8 used
             5 bits (label produced by the MN at this column)
             
             n-4          *
             0  20 40 60 80 100 120 140 160
             1
             2
             3
             4
             5
             6                  126
             7
             8
             9                      149
             10                     150
             11
             12    52
             13       73
             14          94         154
             15          95
             16
             17
             18    58
             19 39 59 79 99 119 139 159 179
             */
            
            // for each update:
            for(int u=0; u<(2*N); ++u) {
                bool changed=false;
                neighborhood n(&r, pt); // point the neighborhood at the right rna and state matrix
                
                for(int a=0; a<N; ++a) {
                    n.reset(a-radius);
                    ca[a].update(n);
                    col_type c1(*ptp1,a), c2(*pt,a);
                    std::size_t oc=0;
                    for(typename EA::phenotype_type::iterator j=ca[a].begin_output(); j!=ca[a].end_output(); j+=2, ++oc) {
                        c1[oc] = (*j) & (!(*(j+1)));
                        changed = changed || (c1[oc] != c2[oc]);
                    }
                }
                
                // rotate the state vector; BE SURE to use pt from here on!
                std::swap(pt, ptp1);
                
                // early stopping:
                if(!changed) {
                    break;
                }
            }
            
            // record the state:
            if(_cb != 0) {
                _cb->classification(r, *pt);
            }

            // train or test?
            std::vector<roc>* A=0;
            if(i <= (get<RNA_CROSSVAL>(ea,0.8) * _rna.size())) {
                A = &ind.traits().train;
            } else {
                A = &ind.traits().test;
            }
            
            // fitness contribution:
            for(std::size_t j=0; j<r.output_size2(); ++j) { // columns
                // check to see if we do better w/ restricted outputs or not:
                int decision=0;
                for(std::size_t k=0; k<3; ++k) {
                    decision |= (*pt)(k,j) << k;
                }
                for(std::size_t k=0; k<r.output_size1(); ++k) { // rows
                    int c=r.N(k,j); // condition
                    int t=(decision == k); // test outcome
                    (*A)[k](c, t);
                }

//                this is the multi-output setup:
//                for(std::size_t k=0; k<r.output_size1(); ++k) { // rows
//                    int c=r.N(k,j); // condition
//                    int t=(*pt)(k,j); // test outcome
//                    (*A)[k](c, t);
//                }
            }
        }
        
        std::vector<roc>* A=&ind.traits().train;
        value_type f;
        for(std::size_t i=0; i<A->size(); ++i) {
            double mcc = (*A)[i].mcc();
            f.push_back(mcc+1.0);
        }
        
        return f;
    }
};

#endif
