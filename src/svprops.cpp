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
#include "htslib/tbx.h"


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
  
  // Load tabix index file
  tbx_t* tbx = tbx_index_load(argv[1]);
  if (!tbx) {
    std::cerr << "Fail to load tabix index file for " << argv[1] << "!" << std::endl;
    return 1;
  }

  bcf_hdr_t* hdr = bcf_hdr_read(ifile);
  bcf1_t* rec = bcf_init();


  int32_t nsvend = 0;
  int32_t* svend = NULL;
  int32_t nsvt = 0;
  char* svt = NULL;
  int ngt = 0;
  int32_t* gt = NULL;
  int ndv = 0;
  int32_t* dv = NULL;
  int ndr = 0;
  int32_t* dr = NULL;
  std::cout << "chr\tstart\tend\tid\tsvType\tvac\taf\tsize\thetperatio\tsingleton\tmissingrate" << std::endl;
  while (bcf_read(ifile, hdr, rec) == 0) {
    int ac[3];
    ac[0] = 0;
    ac[1] = 0;
    ac[2] = 0;
    int uncalled = 0;
    bcf_unpack(rec, BCF_UN_ALL);
    bcf_get_info_int32(hdr, rec, "END", &svend, &nsvend);
    bcf_get_info_string(hdr, rec, "SVTYPE", &svt, &nsvt);
    int ngts = bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt);
    bcf_get_format_int32(hdr, rec, "DV", &dv, &ndv);
    bcf_get_format_int32(hdr, rec, "DR", &dr, &ndr);
    int gtPerSample = ngts / bcf_hdr_nsamples(hdr);
    std::string rareCarrier;
    typedef std::vector<double> TDVector;
    TDVector hetPERatio;
    for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i) {
      for (int j = 0; j < gtPerSample; ++j) {
	if (bcf_gt_allele(gt[i*2 + j]) != -1) {
	  if (bcf_gt_allele(gt[i*2 + j]) != 0) rareCarrier = hdr->samples[i];
	  if ((bcf_gt_allele(gt[i*2 + j]) == 1) && (dv[i]>0) && (dr[i]>0)) hetPERatio.push_back( (double) dv[i] / (double) (dr[i] + dv[i]) );
	  ++ac[bcf_gt_allele(gt[i*2 + j])];
	} else {
	  ++uncalled;
	}
      }
    }
    int vac=ac[1] + ac[2];
    double af= (double) vac / (double) (ac[0] + ac[1] + ac[2]);
    int svlen = *svend - rec->pos;
    double missingRate = (double) uncalled / (double) ngts;
    if (vac > 1) rareCarrier = "NA";
    double hper = 0;
    if (hetPERatio.size()) {
      TDVector::iterator begin = hetPERatio.begin();
      TDVector::iterator end = hetPERatio.end();
      std::nth_element(begin, begin + (end - begin) / 2, end);
      hper = *(begin + (end - begin) / 2);
    }
    std::cout << bcf_hdr_id2name(hdr, rec->rid) << "\t" << (rec->pos + 1) << "\t" << *svend << "\t" << rec->d.id << "\t" << svt << "\t" << vac << "\t" << af  << "\t" << svlen << "\t" << hper << "\t" << rareCarrier << "\t" << missingRate << std::endl;
  }

  // Clean-up
  free(svend);
  free(svt);
  free(gt);
  free(dv);
  free(dr);
  tbx_destroy(tbx);
  bcf_hdr_destroy(hdr);
  bcf_close(ifile);
  bcf_destroy(rec);

  return 0;
}
