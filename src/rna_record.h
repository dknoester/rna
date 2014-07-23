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
#ifndef _RNA_RECORD_H_
#define _RNA_RECORD_H_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/regex.hpp>

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

typedef boost::numeric::ublas::matrix<int> matrix_type; //!< Type for matrix that will store raw input.
typedef boost::numeric::ublas::matrix_row<matrix_type> row_type; //!< Row type for the matrix.
typedef boost::numeric::ublas::matrix_column<matrix_type> col_type; //!< Column type for the matrix.
typedef boost::numeric::ublas::scalar_matrix<int> scalar_matrix_type; //!< Scalar matrix type to aid in initialization.
typedef boost::numeric::ublas::vector<int> vector_type; //!< Vector type.
typedef boost::numeric::ublas::scalar_vector<int> scalar_vector_type; //!< Vector type.

//! Single RNA protein binding record.
struct rna_record {
    std::string S; //!< Sequence (A,C,G,U,I,X): 6
    std::string T; //!< Structure (a-i): 9
    std::string L; //!< Label (S,N,A,X,B): 5
    
    // the matrices look like this:
    //
    // sequence:   a c g t c g t g c a i x ...
    // structure:  a e i c d b c b a i h f ...
    // label:      s n a x b n x b s a n x ...
    //
    // input matrix: (#sequence+#structure, #bases)
    //   seq: a    1 0 0 0 0 0 0 0 0 1 0 0 ...
    //        c    0 1 0 0 1 0 0 0 1 0 0 0 ...
    //      ...
    //        x    0 0 0 0 0 0 0 0 0 0 0 1 ...
    //   str: a    1 0 0 0 0 0 0 0 1 0 0 0 ...
    //        b    0 0 0 0 0 1 0 1 0 0 0 0 ...
    //      ...
    //        i    0 0 1 0 0 0 0 0 0 1 0 0 ...
    //
    // label matrix (#labels, #bases)
    //   lbl: s    1 0 0 0 0 0 0 0 1 0 0 0 ...
    //        n    0 1 0 0 0 1 0 0 0 0 1 0 ...
    //      ...
    //        b    0 0 0 0 1 0 0 1 0 0 0 0 ...
    
    matrix_type M; //!< Pre-processed input data for this record
    matrix_type N; //!< Pre-processed label data for this record
    std::vector<int> NL; //!< Numeric labels
    
    //! Constructor for RNA records.
    template <typename RandomAccess, typename EA>
    rna_record(RandomAccess ri, EA& ea) : S(ri[0]), T(ri[1]), L(ri[2]) {
        using namespace boost;
        assert(S.size() == T.size());
        assert(S.size() == L.size());
        
        M = scalar_matrix_type(15, S.size(), 0); // 6 + 9 bits for sequence+structure
        N = scalar_matrix_type(5, S.size(), 0); // 5 bits for label
        
        for(std::size_t i=0; i<M.size2(); ++i) {
            int seq=seq2int(std::string("") + S[i]);
            int str=str2int(std::string("") + T[i]);
            int lbl=lbl2int(std::string("") + L[i]);
            M(seq,i) = 1;
            M(6 + str,i) = 1; // +6 for the seq
            N(lbl,i) = 1;
            NL.push_back(lbl);
        }
    }
    
    std::size_t input_size1() const { return M.size1(); }
    std::size_t input_size2() const { return M.size2(); }
    std::size_t output_size1() const { return N.size1(); }
    std::size_t output_size2() const { return N.size2(); }
    
    
    int seq2int(const std::string& s) {
        using namespace boost;
        static const regex a("a", regex::icase);
        static const regex c("c", regex::icase);
        static const regex g("g", regex::icase);
        static const regex u("u", regex::icase);
        static const regex i("i", regex::icase);
        static const regex x("x", regex::icase);
        
        if(regex_match(s,i)) {
            return 0x00;
        } else if(regex_match(s,x)) {
            return 0x01;
        } else if(regex_match(s,a)) {
            return 0x02;
        } else if(regex_match(s,c)) {
            return 0x03;
        } else if(regex_match(s,g)) {
            return 0x04;
        } else if(regex_match(s,u)) {
            return 0x05;
        } else {
            throw bad_argument_exception("rna.h::seq2int: invalid base " + s);
        }
    }
    
    int str2int(const std::string& l) {
        using namespace boost;
        
        static const regex W("a", regex::icase);
        static const regex WS("b", regex::icase);
        static const regex Ws("c", regex::icase);
        static const regex WSs("d", regex::icase);
        static const regex S("e", regex::icase);
        static const regex Ss("f", regex::icase);
        static const regex s("g", regex::icase);
        static const regex u("h", regex::icase);
        static const regex i("i", regex::icase);
        
        if(regex_match(l,i)) {
            return 0x00;
        } else if(regex_match(l,W)) {
            return 0x01;
        } else if(regex_match(l,WS)) {
            return 0x02;
        } else if(regex_match(l,Ws)) {
            return 0x03;
        } else if(regex_match(l,WSs)) {
            return 0x04;
        } else if(regex_match(l,S)) {
            return 0x05;
        } else if(regex_match(l,Ss)) {
            return 0x06;
        } else if(regex_match(l,s)) {
            return 0x07;
        } else if(regex_match(l,u)) {
            return 0x08;
        } else {
            throw bad_argument_exception("rna.h::str2int: invalid structure " + l);
        }
    }
    
    int lbl2int(const std::string& l) {
        using namespace boost;
        
        static const regex s("s", regex::icase);
        static const regex n("n", regex::icase);
        static const regex a("a", regex::icase);
        static const regex x("x", regex::icase);
        static const regex b("b", regex::icase);
        
        if(regex_match(l,x)) {
            return 0x00;
        } else if(regex_match(l,s)) {
            return 0x01;
        } else if(regex_match(l,n)) {
            return 0x02;
        } else if(regex_match(l,b)) {
            return 0x03;
        } else if(regex_match(l,a)) {
            return 0x04;
        } else {
            throw bad_argument_exception("rna.h::lbl2int: invalid binding " + l);
        }
    }
};

#endif
