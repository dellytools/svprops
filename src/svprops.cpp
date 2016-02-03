/*
============================================================================
SV VCF properties
============================================================================
Copyright (C) 2015 Tobias Rausch

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
============================================================================
Contact: Tobias Rausch (rausch@embl.de)
============================================================================
*/

#include <iostream>
#include <set>
#include <vector>
#include <algorithm>
#include "htslib/vcf.h"


template<typename TVector>
inline void
_getMedian(TVector& v, typename TVector::value_type& med) {
  med = 0;
  if (v.size()) {
    typename TVector::iterator begin = v.begin();
    typename TVector::iterator end = v.end();
    std::nth_element(begin, begin + (end - begin) / 2, end);
    med = *(begin + (end - begin) / 2);
  }
}


int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <in.vcf.gz> " << std::endl;
    return 1; 
  }

  // Load bcf file
  htsFile* ifile = bcf_open(argv[1], "r");
  if (!ifile) {
    std::cerr << "Fail to load " << argv[1] << "!" << std::endl;
    return 1;
  }
  
  bcf_hdr_t* hdr = bcf_hdr_read(ifile);
  bcf1_t* rec = bcf_init();


  int32_t nsvend = 0;
  int32_t* svend = NULL;
  int32_t ninslen = 0;
  int32_t* inslen = NULL;
  int32_t ncipos = 0;
  int32_t* cipos = NULL;
  int32_t nsvt = 0;
  char* svt = NULL;
  int ngt = 0;
  int32_t* gt = NULL;
  int nrc = 0;
  int32_t* rc = NULL;
  int nrcl = 0;
  int32_t* rcl = NULL;
  int nrcr = 0;
  int32_t* rcr = NULL;
  int ndv = 0;
  int32_t* dv = NULL;
  int ndr = 0;
  int32_t* dr = NULL;
  int nrv = 0;
  int32_t* rv = NULL;
  int nrr = 0;
  int32_t* rr = NULL;
  int ngq = 0;
  int32_t* gq = NULL;
  int nft = 0;
  char** ft = NULL;
  std::cout << "chr\tstart\tend\tid\tsvType\tprecise\tvac\tvaf\tsize\tci\trefratio\taltratio\trefgq\taltgq\trefpass\taltpass\trdratio\tmedianrc\tsingleton\tmissingrate" << std::endl;
  while (bcf_read(ifile, hdr, rec) == 0) {
    int ac[2];
    ac[0] = 0;
    ac[1] = 0;
    int uncalled = 0;
    bcf_unpack(rec, BCF_UN_ALL);
    bcf_get_info_int32(hdr, rec, "END", &svend, &nsvend);
    bcf_get_info_int32(hdr, rec, "INSLEN", &inslen, &ninslen);
    bcf_get_info_int32(hdr, rec, "CIPOS", &cipos, &ncipos);
    bcf_get_info_string(hdr, rec, "SVTYPE", &svt, &nsvt);
    bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt);
    bcf_get_format_int32(hdr, rec, "RC", &rc, &nrc);
    bcf_get_format_int32(hdr, rec, "RCL", &rcl, &nrcl);
    bcf_get_format_int32(hdr, rec, "RCR", &rcr, &nrcr);
    bcf_get_format_int32(hdr, rec, "DV", &dv, &ndv);
    bcf_get_format_int32(hdr, rec, "DR", &dr, &ndr);
    bcf_get_format_int32(hdr, rec, "RV", &rv, &nrv);
    bcf_get_format_int32(hdr, rec, "RR", &rr, &nrr);
    bcf_get_format_int32(hdr, rec, "GQ", &gq, &ngq);
    bcf_get_format_string(hdr, rec, "FT", &ft, &nft);
    bool precise = false;
    if (bcf_get_info_flag(hdr, rec, "PRECISE", 0, 0) > 0) precise = true;
    std::string rareCarrier;
    typedef std::vector<double> TDVector;
    typedef std::vector<int32_t> TIVector;
    TIVector gqRef;
    TIVector gqAlt;
    TDVector ratioRef;
    TDVector ratioAlt;
    TDVector rcRefRatio;
    TDVector rcAltRatio;
    TIVector rcRef;
    int32_t passRef = 0;
    int32_t passAlt = 0;
    for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i) {
      if ((bcf_gt_allele(gt[i*2]) != -1) && (bcf_gt_allele(gt[i*2 + 1]) != -1)) {
        int gt_type = bcf_gt_allele(gt[i*2]) + bcf_gt_allele(gt[i*2 + 1]);
	++ac[bcf_gt_allele(gt[i*2])];
	++ac[bcf_gt_allele(gt[i*2 + 1])];
	std::string pStr("PASS");
	bool pass = true;
	if (pStr.compare(ft[i]) != 0) pass = false;
	if (gt_type == 0) {
	  // Non-carrier
	  gqRef.push_back( gq[i] );
	  rcRef.push_back( rc[i] );
	  rcRefRatio.push_back( (double) rc[i] / (double) (rcl[i] + rcr[i]) );
	  if (pass) ++passRef;
	  if (precise) ratioRef.push_back( (double) rv[i] / (double) (rr[i] + rv[i]) );
	  else ratioRef.push_back( (double) dv[i] / (double) (dr[i] + dv[i]) );
	} else {
	  // Carrier
	  if (ac[1] == 1) rareCarrier = hdr->samples[i];
	  gqAlt.push_back( gq[i] );
	  if (pass) ++passAlt;
	  if (gt_type == 1) {
	    if (precise) ratioAlt.push_back( (double) rv[i] / (double) (rr[i] + rv[i]) );
	    else ratioAlt.push_back( (double) dv[i] / (double) (dr[i] + dv[i]) );
	    rcAltRatio.push_back( (double) rc[i] / (double) (rcl[i] + rcr[i]) );
	  }
	}
      } else {
        ++uncalled;
      }
    }
    double af= (double) ac[1] / (double) (ac[0] + ac[1]);
    int svlen = *svend - rec->pos;
    if (std::string(svt) == "INS") svlen = *inslen;
    double missingRate = (double) uncalled / (double) bcf_hdr_nsamples(hdr);
    if (ac[1] != 1) rareCarrier = "NA";
    TDVector::value_type refratio = 0;
    _getMedian(ratioRef, refratio);
    TDVector::value_type altratio = 0;
    _getMedian(ratioAlt, altratio);
    TIVector::value_type refgq = 0;
    _getMedian(gqRef, refgq);
    TIVector::value_type altgq = 0;
    _getMedian(gqAlt, altgq);
    TDVector::value_type altRC = 0;
    _getMedian(rcAltRatio, altRC);
    TDVector::value_type refRC = 0;
    _getMedian(rcRefRatio, refRC);
    TDVector::value_type rdRatio = altRC/refRC;
    TIVector::value_type rcMed = 0;
    _getMedian(rcRef, rcMed);
    std::cout << bcf_hdr_id2name(hdr, rec->rid) << "\t" << (rec->pos + 1) << "\t" << *svend << "\t" << rec->d.id << "\t" << svt << "\t" << precise << "\t" << ac[1] << "\t" << af  << "\t" << svlen << "\t" << cipos[1] << "\t" << refratio << "\t" << altratio << "\t" << refgq << "\t" << altgq << "\t" << (double) passRef / (double) gqRef.size() << "\t" << (double) passAlt / (double) gqAlt.size() << "\t" << rdRatio << "\t" << rcMed << "\t" << rareCarrier << "\t" << missingRate << std::endl;
  }

  // Clean-up
  free(svend);
  free(inslen);
  free(cipos);
  free(svt);
  free(gt);
  free(rc);
  free(rcl);
  free(rcr);
  free(dv);
  free(dr);
  free(rv);
  free(rr);
  free(gq);
  if (ft != NULL) {
    free(ft[0]);
    free(ft);
  }
  bcf_hdr_destroy(hdr);
  bcf_close(ifile);
  bcf_destroy(rec);

  return 0;
}
