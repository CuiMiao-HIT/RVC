/*
 * callerstep1.cpp
 *
 *  Created on: 2023年02月15日
 *      Author: miaocui
 */

#include "callerstep1A3.hpp"

opt_t *opt;

int run_callerStep3(int argc, char *argv[]){
    double cpu_time_all = cputime_hp();
    struct timeval start_all;
    gettimeofday(&start_all, NULL);

    opt_t opt_ent;   
    opt = &opt_ent; //global opt
    if (opt_parse(argc, argv, opt) != 0)  return 2;
    std :: vector <VRT_info> VRT_info_p3;
    //load vrt info from step1 
    ifstream ifvrtinfo;
    ifvrtinfo.open(opt->snvinfo_path, ios::in);
    std :: string vrtline;
    while (getline(ifvrtinfo , vrtline)){
        stringstream bed_ss(vrtline);
        string sp;
        while(bed_ss >> sp){
            VRT_info add_one;
            std :: vector<string> SP = split2(sp,"_");
            // add_one.chrID = stoi(SP[0]);
            add_one.ref_pos = stoi(SP[1]);
            add_one.vrttype = stoi(SP[2]);
            add_one.REF = SP[3];
            add_one.length_in_ref = SP[3].size();
            add_one.ALT = SP[4];
            add_one.length_in_read = SP[4].size();
            add_one.all_num = stoi(SP[5]);
            add_one.ref_num = stoi(SP[6]);
            add_one.vrt_num = stoi(SP[7]);
            add_one.vrt_num0 = stoi(SP[8]);
            add_one.vrt_num16 = stoi(SP[9]);
            VRT_info_p3.push_back(add_one);
        }
    }
    ifvrtinfo.close();
    int Pvcf = printALLsnv(VRT_info_p3);

    fprintf(stderr, "\nClassify CPU STEPcm of [%s] : %.3f sec\n\n", opt->chrID, cputime_hp() - cpu_time_all);
    return 0;
}

int printALLsnv(std :: vector <VRT_info> VRT_info_p3){
    ofstream bedfile;
    bedfile.open(opt->bed_path);
    std :: vector <result_VRT_info> Result_vcfs;
    for (int i = 0; i < VRT_info_p3.size(); i++){
        int VaRnum = VRT_info_p3[i].vrt_num + VRT_info_p3[i].ref_num;
        int n_alts = VRT_info_p3[i].vrt_num, n_total = VaRnum;
        float error_rate;
        if(VRT_info_p3[i].vrttype == 0){
            error_rate = ERR_snp;
            std ::vector<float> GL_P;
            rescale_read_counts(n_alts, n_total);
            GL_P = cal_GL(n_alts, n_total, error_rate);
            FormatInfo outinfo = Likelihood_Estimation(GL_P);
            // std::vector<float> priors = {0.0f, 0.4f, 0.6f};
            // FormatInfo outinfo = genotype_by_majority_vote_with_quality(GL_P, priors, VRT_info_p3[i].ref_num, n_alts, error_rate);
            // if(outinfo.GT == "0/0"||outinfo.QUAL == 0 || outinfo.GQ == 0) {
            if(outinfo.GT == "0/0"||outinfo.QUAL < opt->Qual || outinfo.GQ == 0) {
                bedfile << opt->chrID << "\t" <<  VRT_info_p3[i].ref_pos << "\t" <<  VRT_info_p3[i].ref_pos <<endl;
                continue;
            }
            else {
                result_VRT_info addvrt;
                addvrt.ref_pos = VRT_info_p3[i].ref_pos;
                addvrt.REF = VRT_info_p3[i].REF;
                addvrt.ALT = VRT_info_p3[i].ALT;
                addvrt.outinfo = outinfo;
                Result_vcfs.push_back(addvrt);
            }
        }else{
            error_rate = ERR_indel;
            std ::vector<float> GL_P;
            rescale_read_counts(n_alts, n_total);
            GL_P = cal_GL(n_alts, n_total, error_rate);
            std::vector<float> priors = {0.0f, 0.5f, 0.5f};
            FormatInfo outinfo = genotype_by_majority_vote_with_quality(GL_P, priors, VRT_info_p3[i].ref_num, n_alts, error_rate);
            if(outinfo.GT == "0/0"||outinfo.QUAL < opt->Qual || outinfo.GQ == 0) {
                bedfile << opt->chrID << "\t" <<  VRT_info_p3[i].ref_pos << "\t" <<  VRT_info_p3[i].ref_pos <<endl;
                continue;
            }
            else {
                result_VRT_info addvrt;
                addvrt.ref_pos = VRT_info_p3[i].ref_pos;
                addvrt.REF = VRT_info_p3[i].REF;
                addvrt.ALT = VRT_info_p3[i].ALT;
                addvrt.outinfo = outinfo;
                Result_vcfs.push_back(addvrt);
            }
        }
    }

    std :: sort(Result_vcfs.begin() , Result_vcfs.end());
    std :: vector <uint32_t> ResultSnv_refpos;
    for (auto cc : Result_vcfs) ResultSnv_refpos.push_back(cc.ref_pos);
    ResultSnv_refpos.push_back(4294967294);
    sort(ResultSnv_refpos.begin() , ResultSnv_refpos.end());

    ofstream ALLsnvfile;
    ALLsnvfile.open(opt->snvvcf_path);
    // body
    int snvpos_next = 1,pass = 0;
    for(auto onelow : Result_vcfs){
        if (pass == 1) {
            snvpos_next = snvpos_next + 1;
            pass = 0;
            continue;
        }
        uint32_t nextpos = ResultSnv_refpos[snvpos_next];
        std::string Filter = "PASS";
        if (onelow.ref_pos == nextpos) {//the same pos has more than one vrt
            ALLsnvfile << opt->chrID << "\t" << onelow.ref_pos << "\t"
                << ".\t" << onelow.REF << "\t" << onelow.ALT << "\t" << onelow.outinfo.QUAL << "\t"
                << Filter
                << "\t"
                << "."
                << "\t"
                << "GT:GQ"
                << "\t" << "0/1" << ":" << onelow.outinfo.GQ << endl;
            ALLsnvfile << opt->chrID << "\t" << onelow.ref_pos << "\t"
                << ".\t" << Result_vcfs[snvpos_next].REF << "\t" << Result_vcfs[snvpos_next].ALT << "\t" << Result_vcfs[snvpos_next].outinfo.QUAL << "\t"
                << Filter
                << "\t"
                << "."
                << "\t"
                << "GT:GQ"
                << "\t" << "1/0" << ":" << Result_vcfs[snvpos_next].outinfo.GQ << endl;
            pass = 1;
            snvpos_next = snvpos_next + 1;
        }else{// 0/1 or 1/1
            ALLsnvfile << opt->chrID << "\t" << onelow.ref_pos << "\t"
                << ".\t" << onelow.REF << "\t" << onelow.ALT << "\t" << onelow.outinfo.QUAL << "\t"
                << Filter
                << "\t"
                << "."
                << "\t"
                << "GT:GQ"
                << "\t" << onelow.outinfo.GT << ":" << onelow.outinfo.GQ << endl;
            snvpos_next = snvpos_next + 1;
        }
    }
    ResultSnv_refpos.clear();
    ALLsnvfile.close();
    bedfile.close();
    return 0;
}
