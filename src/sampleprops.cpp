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
#include <vector>
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

  typedef std::vector<unsigned int> TCountMap;
  TCountMap missing;
  TCountMap ref;
  TCountMap het;
  TCountMap hom;
  // Initialization
  missing.resize(bcf_hdr_nsamples(hdr), 0);
  ref.resize(bcf_hdr_nsamples(hdr), 0);
  het.resize(bcf_hdr_nsamples(hdr), 0);
  hom.resize(bcf_hdr_nsamples(hdr), 0);

  // Counting
  int ngt = 0;
  int32_t* gt = NULL;
  std::cout << "sample\tmissing\thomref\thet\thomalt" << std::endl;
  while (bcf_read(ifile, hdr, rec) == 0) {
    bcf_unpack(rec, BCF_UN_ALL);
    bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt);
    for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i) {
      if ((bcf_gt_allele(gt[i*2]) != -1) && (bcf_gt_allele(gt[i*2 + 1]) != -1)) {
	int gt_type = bcf_gt_allele(gt[i*2]) + bcf_gt_allele(gt[i*2 + 1]);
	if (gt_type == 0) ++ref[i];
	else if (gt_type == 1) ++het[i];
	else ++hom[i];
      } else {
	++missing[i];
      }
    }
  }

  // Output
  for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i) std::cout << hdr->samples[i] << '\t' << missing[i] << '\t' << ref[i] << '\t' << het[i] << '\t' << hom[i] << std::endl;

  // Clean-up
  free(gt);
  tbx_destroy(tbx);
  bcf_hdr_destroy(hdr);
  bcf_close(ifile);
  bcf_destroy(rec);

  return 0;
}
