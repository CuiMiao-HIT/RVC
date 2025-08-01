/*
 * SNV_code.hpp
 *
 *  Created on: 2023年02月15日
 *      Author: miaocui
 */

#ifndef CALLERSTEP1A3_HPP_
#define CALLERSTEP1A3_HPP_
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <map>
#include <bitset>
#include <cmath>
#include <regex>
#include <cctype>
#include <numeric>
#include <ctime>
#include <cstdio>

#include "../cmtools/math_func.hpp"
#include "../cmtools/form.hpp"

extern "C"
{
    #include "../clib/utils.h"
}

using namespace std;

#define ERR_snp 0.08;//0.08
#define ERR_indel 0.12;//0.12


typedef struct VRT_info {
    // int chrID;
    uint32_t ref_pos;
    int vrttype;//0:X 1:I|D
    std :: string REF = "";
    std :: string ALT = "";
    int length_in_ref;
    int length_in_read;
    int vrt_num = 0;//reads support the vrt
    int vrt_num0 = 0;
    int vrt_num16 = 0;
    int ref_num = 0;//reads support the vrt and reference 
    int all_num = 0;//all reads support this pos

    bool operator== (const VRT_info &a) const{
        return (ref_pos==a.ref_pos) && (length_in_ref == a.length_in_ref) && (ALT == a.ALT);
    }
    bool operator< (const VRT_info &a) const{
        return ref_pos < a.ref_pos;
    }
};

typedef struct result_VRT_info {
    uint32_t ref_pos;
    std :: string REF = "";
    std :: string ALT = "";
    FormatInfo outinfo;

    bool operator== (const result_VRT_info &a) const{
        return (ref_pos==a.ref_pos) && (REF == a.REF) && (ALT == a.ALT);
    }
    bool operator< (const result_VRT_info &a) const{
        return ref_pos < a.ref_pos;
    }
};


int run_callerStep3(int argc, char *argv[]);
int printALLsnv(std :: vector <VRT_info> VRT_info_p3);
#endif /* CALLERSTEP1A3_HPP_ */