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
#ifndef _RNA_H_
#define _RNA_H_

#include <ea/fitness_function.h>

#include "rna_record.h"
#include "analysis.h"

LIBEA_MD_DECL(RNA_INPUT, "rna.input_filename", std::string); // name of input datafile
LIBEA_MD_DECL(RNA_CROSSVAL, "rna.cross_validation", double); // cross validation fraction
LIBEA_MD_DECL(RNA_LIMIT, "rna.limit", int); // number to stop eval at
LIBEA_MD_DECL(CA_RADIUS, "cellular_automata.radius", int);

/*! Constructs a CA neighborhood for a given RNA record and state matrix.
 */
struct neighborhood {
    neighborhood(rna_record* r, matrix_type* m) : _r(r), _m(m), _n(0) {
        _nrowr = _r->input_size1();
        _nrow = _r->input_size1() + _m->size1();
        _ncol = _r->input_size2();
        _size = _nrow * _ncol;
    }
    
    void reset(int i) { // i is number of columns away from 0 (can be negative); must * by nrows
        _n = i * _nrow;
    }
    
    inline int operator[](int i) {
        i += _n;
        if((i < 0) || (i >= _size)) {
            return 0;
        }
        
        // in bounds:
        int c = i / _nrow;
        int r = i % _nrow;
        
        if(r < _nrowr) {
            return _r->M(r,c);
        } else {
            return (*_m)(r-_nrowr,c);
        }
    }
    
    rna_record* _r; // rna record
    matrix_type* _m; // matrix for state information
    int _n; // base index for this neighborhood
    int _nrowr, _nrow, _ncol, _size;
};


/*! Fitness function for RNA labeling project.
 */
struct rna_fitness : fitness_function<unary_fitness<double>, constantS, stochasticS> {
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
    }
    
    //! Calculate fitness of an individual.
    template <typename Individual, typename RNG, typename EA>
    double operator()(Individual& ind, RNG& rng, EA& ea) {
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
        
        double f=0.0;
        matrix_type F=scalar_matrix_type(5,2,0);
        
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
            
            // fitness contribution:
            for(std::size_t j=0; j<r.output_size2(); ++j) { // num columns
                bool correct=true;
                int nl=r.NL[j];
                for(std::size_t k=0; k<r.output_size1(); ++k) { // num rows
                    int c=(*pt)(k,j);
                    int l=r.N(k,j);
                    correct = correct && (c==l);
                }
                if(correct) {
                    f += exp((-1.0/25.0) * static_cast<double>(++F(nl,0)));
                }
            }
        }
        
        return f;
    }
};

#endif
