#include <iostream>
#include <set>
#include <vector>
#include <algorithm>
#include <map>
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

template<typename TVector>
inline void
_getMax(TVector& v, typename TVector::value_type& max) {
  if (v.size()) {
    max = v[0];
    for(typename TVector::iterator itv = v.begin();  itv != v.end(); ++itv)
      if (*itv > max) max = *itv;
  }
}

inline bool
_isKeyPresent(bcf_hdr_t const* hdr, std::string const& key) {
  return (bcf_hdr_id2int(hdr, BCF_DT_ID, key.c_str())>=0);
}

inline int
_getInfoType(bcf_hdr_t const* hdr, std::string const& key) {
  return bcf_hdr_id2type(hdr, BCF_HL_INFO, bcf_hdr_id2int(hdr, BCF_DT_ID, key.c_str()));
}

inline int
_getFormatType(bcf_hdr_t const* hdr, std::string const& key) {
  return bcf_hdr_id2type(hdr, BCF_HL_FMT, bcf_hdr_id2int(hdr, BCF_DT_ID, key.c_str()));
}

inline bool _missing(bool const value) {
  return !value;
}

inline bool _missing(float const value) {
  return bcf_float_is_missing(value);
}

inline bool _missing(int8_t const value) {
  return (value == bcf_int8_missing);
}

inline bool _missing(int16_t const value) {
  return (value == bcf_int16_missing);
}

inline bool _missing(int32_t const value) {
  return (value == bcf_int32_missing);
}

inline bool _missing(std::string const& value) {
  return ((value.empty()) || (value == "."));
}


int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <in.vcf.gz> " << std::endl;
    return 1; 
  }

  // Load bcf file
  htsFile* ifile = bcf_open(argv[1], "r");
  if (ifile == NULL) {
    std::cerr << "Fail to load " << argv[1] << std::endl;
    return 1;
  }
  bcf_hdr_t* hdr = bcf_hdr_read(ifile);

  // Read SV information
  int32_t nsvend = 0;
  int32_t* svend = NULL;
  int32_t ninslen = 0;
  int32_t* inslen = NULL;
  int32_t npos2 = 0;
  int32_t* pos2 = NULL;
  int32_t nhomlen = 0;
  int32_t* homlen = NULL;
  int32_t ncipos = 0;
  int32_t* cipos = NULL;
  int32_t nsvt = 0;
  char* svt = NULL;
  int32_t nfic = 0;
  float* fic = NULL;
  int32_t nce = 0;
  float* ce = NULL;
  int32_t nexchet = 0;
  float* exchet = NULL;
  int32_t nhwepval = 0;
  float* hwepval = NULL;
  int32_t nchr2 = 0;
  char* chr2 = NULL;
  int32_t nct = 0;
  char* ct = NULL;
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
  int32_t* gqInt = NULL;
  float* gqFloat = NULL;

  // Get the valid columns
  uint32_t fieldIndex = 0;
  typedef std::map<std::string, uint32_t> TColumnMap;
  TColumnMap cMap;
  cMap["chr"] = fieldIndex++;
  cMap["start"] = fieldIndex++;
  cMap["chr2"] = fieldIndex++;
  cMap["end"] = fieldIndex++;
  cMap["id"] = fieldIndex++;
  cMap["size"] = fieldIndex++;
  cMap["qual"] = fieldIndex++;
  cMap["vac"] = fieldIndex++;
  cMap["vaf"] = fieldIndex++;
  cMap["pass"] = fieldIndex++;
  cMap["singleton"] = fieldIndex++;
  cMap["missingrate"] = fieldIndex++;
  if (_isKeyPresent(hdr, "SVTYPE")) cMap["svtype"] = fieldIndex++;
  if (_isKeyPresent(hdr, "CT")) cMap["ct"] = fieldIndex++;
  if (_isKeyPresent(hdr, "IMPRECISE")) cMap["precise"] = fieldIndex++;
  if (_isKeyPresent(hdr, "CIPOS")) cMap["ci"] = fieldIndex++;
  if (_isKeyPresent(hdr, "INSLEN")) cMap["inslen"] = fieldIndex++;
  if (_isKeyPresent(hdr, "HOMLEN")) cMap["homlen"] = fieldIndex++;
  if (_isKeyPresent(hdr, "FIC")) cMap["fic"] = fieldIndex++;
  if (_isKeyPresent(hdr, "CE")) cMap["ce"] = fieldIndex++;
  if (_isKeyPresent(hdr, "ExcHet")) cMap["exchet"] = fieldIndex++;
  if (_isKeyPresent(hdr, "HWE")) cMap["hwepval"] = fieldIndex++;
  if (_isKeyPresent(hdr, "GQ")) {
    cMap["refgq"] = fieldIndex++;
    cMap["altgq"] = fieldIndex++;
    cMap["gqsum"] = fieldIndex++;
    cMap["altgqsum"] = fieldIndex++;
  }
  if (_isKeyPresent(hdr, "RC")) {
    cMap["rdratio"] = fieldIndex++;
    cMap["medianrc"] = fieldIndex++;
  }
  if (_isKeyPresent(hdr, "DV")) {
    cMap["refratio"] = fieldIndex++;
    cMap["altratio"] = fieldIndex++;
    cMap["maxaltratio"] = fieldIndex++;
    cMap["PEsupport"] = fieldIndex++;
    cMap["SRsupport"] = fieldIndex++;
    cMap["supportsum"] = fieldIndex++;
  }

  typedef std::vector<std::string> TColumnHeader;
  TColumnHeader cHeader(cMap.size());
  for(TColumnMap::const_iterator cIt = cMap.begin(); cIt != cMap.end(); ++cIt) cHeader[cIt->second] = cIt->first;

  // Write header
  for(TColumnHeader::const_iterator cHead = cHeader.begin(); cHead != cHeader.end(); ++cHead) {
    if (cHead != cHeader.begin()) std::cout << "\t";
    std::cout << *cHead;
  }
  std::cout << std::endl;

  // Parse VCF records
  bcf1_t* rec = bcf_init();
  uint32_t siteCount = 0;
  while (bcf_read(ifile, hdr, rec) == 0) {
    if (rec->n_allele != 2) continue;
    
    bool gqPresent = false;
    bool rcPresent = false;
    bool dvPresent = false;
    bool ficPresent = false;
    bool exchetPresent = false;
    bool hwePresent = false;
    bool ciPresent = false;
    bool svtPresent = false;
    bool svEndPresent = false;
    int32_t endsv = 0;
    
    ++siteCount;
    bcf_unpack(rec, BCF_UN_ALL);
    bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt);
    if (_isKeyPresent(hdr, "END")) {
      if (bcf_get_info_int32(hdr, rec, "END", &svend, &nsvend) > 0) {
	endsv = *svend;
	svEndPresent = true;
      }
    }
    if (!svEndPresent) {
      std::string refAllele = rec->d.allele[0];
      std::string altAllele = rec->d.allele[1];
      int32_t diff = refAllele.size() - altAllele.size();
      endsv = rec->pos + diff + 2;
    }
    if (_isKeyPresent(hdr, "INSLEN")) bcf_get_info_int32(hdr, rec, "INSLEN", &inslen, &ninslen);
    if (_isKeyPresent(hdr, "HOMLEN")) bcf_get_info_int32(hdr, rec, "HOMLEN", &homlen, &nhomlen);
    if (_isKeyPresent(hdr, "CIPOS")) {
      if (bcf_get_info_int32(hdr, rec, "CIPOS", &cipos, &ncipos) > 0) ciPresent = true;
    }
    if (_isKeyPresent(hdr, "FIC")) {
      if (bcf_get_info_float(hdr, rec, "FIC", &fic, &nfic) > 0) ficPresent = true;
    }
    if (_isKeyPresent(hdr, "CE")) bcf_get_info_float(hdr, rec, "CE", &ce, &nce);
    if (_isKeyPresent(hdr, "ExcHet")) {
      if (bcf_get_info_float(hdr, rec, "ExcHet", &exchet, &nexchet) > 0) exchetPresent = true;
    }
    if (_isKeyPresent(hdr, "HWE")) {
      if (bcf_get_info_float(hdr, rec, "HWE", &hwepval, &nhwepval) > 0) hwePresent = true;
    }
    if (_isKeyPresent(hdr, "SVTYPE")) {
      if (bcf_get_info_string(hdr, rec, "SVTYPE", &svt, &nsvt) > 0) svtPresent = true;
    }
    if (_isKeyPresent(hdr, "GQ")) {
      if (_getFormatType(hdr, "GQ") == BCF_HT_INT) {
	if (bcf_get_format_int32(hdr, rec, "GQ", &gqInt, &ngq) > 0) gqPresent = true;
      } else if (_getFormatType(hdr, "GQ") == BCF_HT_REAL) {
	if (bcf_get_format_float(hdr, rec, "GQ", &gqFloat, &ngq) > 0) gqPresent = true;
      }
    }
    bool passSite = (bcf_has_filter(hdr, rec, const_cast<char*>("PASS"))==1);
    bool precise = false;
    if (bcf_get_info_flag(hdr, rec, "PRECISE", 0, 0) > 0) precise = true;
    if (_isKeyPresent(hdr, "RC")) {
      if (bcf_get_format_int32(hdr, rec, "RC", &rc, &nrc) > 0) rcPresent = true;
    }
    if (_isKeyPresent(hdr, "RCL")) bcf_get_format_int32(hdr, rec, "RCL", &rcl, &nrcl);
    if (_isKeyPresent(hdr, "RCR")) bcf_get_format_int32(hdr, rec, "RCR", &rcr, &nrcr);
    if (_isKeyPresent(hdr, "DV")) {
      if (bcf_get_format_int32(hdr, rec, "DV", &dv, &ndv) > 0) dvPresent = true;
    }
    if (_isKeyPresent(hdr, "DR")) bcf_get_format_int32(hdr, rec, "DR", &dr, &ndr);
    if (_isKeyPresent(hdr, "RV")) bcf_get_format_int32(hdr, rec, "RV", &rv, &nrv);
    if (_isKeyPresent(hdr, "RR")) bcf_get_format_int32(hdr, rec, "RR", &rr, &nrr);
    std::string chr2Name(bcf_hdr_id2name(hdr, rec->rid));
    if (_isKeyPresent(hdr, "CHR2")) {
      if (bcf_get_info_string(hdr, rec, "CHR2", &chr2, &nchr2) > 0) chr2Name = std::string(chr2);
      if (_isKeyPresent(hdr, "POS2")) {
	if (bcf_get_info_int32(hdr, rec, "POS2", &pos2, &npos2) > 0) {
	  endsv = *pos2;
	}
      }
    }
    std::string ctval("NA");
    if (_isKeyPresent(hdr, "CT")) {
      if (bcf_get_info_string(hdr, rec, "CT", &ct, &nct) > 0) ctval = std::string(ct);
    }

    std::set<std::string> rareCarrierSet;
    uint32_t totalPE = 0;
    uint32_t totalSR = 0;
    uint64_t supportsum = 0;
    double gqsum = 0;
    double altgqsum = 0;
    typedef double TPrecision;
    typedef std::vector<TPrecision> TValueVector;
    TValueVector gqRef;   // GQ of non-carriers
    TValueVector gqAlt;   // GQ of het. carriers
    TValueVector ratioRef;  // SV support in non-carriers
    TValueVector ratioAlt;  // SV support in het. carriers
    TValueVector rcRefRatio;  // Read-depth ratio in non-carriers
    TValueVector rcAltRatio;  // Read-depth ratio in het. carriers
    TValueVector rcRef;  // Baseline read count in non-carriers
    int32_t ac[2];
    ac[0] = 0;
    ac[1] = 0;
    int32_t uncalled = 0;
    for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i) {
      if ((bcf_gt_allele(gt[i*2]) != -1) && (bcf_gt_allele(gt[i*2 + 1]) != -1)) {
	if (dvPresent) {
	  if (rr[i] + rv[i] + dr[i] + dv[i] == 0) {
	    // Imputed GT
	    ++uncalled;
	    continue;
	  }
	  if (precise) supportsum += rr[i] + rv[i];
	  else supportsum += dr[i] + dv[i];
	}
        int gt_type = bcf_gt_allele(gt[i*2]) + bcf_gt_allele(gt[i*2 + 1]);
	++ac[bcf_gt_allele(gt[i*2])];
	++ac[bcf_gt_allele(gt[i*2 + 1])];

	if (gt_type == 0) {
	  // Non-carrier
	  if (gqPresent) {
	    if (_getFormatType(hdr, "GQ") == BCF_HT_INT) {
	      if (_missing(gqInt[i])) gqRef.push_back(0);
	      else {
		gqRef.push_back( gqInt[i] );
		gqsum += gqInt[i];
	      }
	    } else if (_getFormatType(hdr, "GQ") == BCF_HT_REAL) {
	      if (_missing(gqFloat[i])) gqRef.push_back(0);
	      else {
		gqRef.push_back( gqFloat[i] );
		gqsum += gqFloat[i];
	      }
	    }
	  } else gqRef.push_back(0);
	  if (rcPresent) {
	    rcRef.push_back( rc[i] );
	    if (_isKeyPresent(hdr, "RCL")) {
	      if ((rcl[i] + rcr[i]) > 0) rcRefRatio.push_back( 2.0 * (double) rc[i] / (double) (rcl[i] + rcr[i]) );
	    } else rcRefRatio.push_back((double) rc[i]);
	  }
	  if (dvPresent) {
	    if (precise) ratioRef.push_back( (double) rv[i] / (double) (rr[i] + rv[i]) );
	    else ratioRef.push_back( (double) dv[i] / (double) (dr[i] + dv[i]) );
	  }
	} else {
	  if (dvPresent) {
	    totalSR += rv[i];
	    totalPE += dv[i];
	  }
	  
	  if (gt_type >= 1) {
	    if (rareCarrierSet.size() < 2) rareCarrierSet.insert(hdr->samples[i]);
	    if (gqPresent) {
	      if (_getFormatType(hdr, "GQ") == BCF_HT_INT) {
		if (_missing(gqInt[i])) gqAlt.push_back(0);
		else {
		  gqAlt.push_back( gqInt[i] );
		  gqsum += gqInt[i];
		  altgqsum += gqInt[i];
		}
	      } else if (_getFormatType(hdr, "GQ") == BCF_HT_REAL) {
		if (_missing(gqFloat[i])) gqAlt.push_back(0);
		else {
		  gqAlt.push_back( gqFloat[i] );
		  gqsum += gqFloat[i];
		  altgqsum += gqFloat[i];
		}
	      }
	    } else gqAlt.push_back(0);
	    if (rcPresent) {
	      if (_isKeyPresent(hdr, "RCL")) {
		if ((rcl[i] + rcr[i]) > 0) rcAltRatio.push_back( 2.0 * (double) rc[i] / (double) (rcl[i] + rcr[i]) );
	      } else rcAltRatio.push_back((double) rc[i]);
	    }
	    if (dvPresent) {
	      if (precise) ratioAlt.push_back( (double) rv[i] / (double) (rr[i] + rv[i]) );
	      else ratioAlt.push_back( (double) dv[i] / (double) (dr[i] + dv[i]) );
	    }
	  }
	}
      } else ++uncalled;
    }
    std::string rareCarrier = "NA";
    if (rareCarrierSet.size() == 1) rareCarrier = *(rareCarrierSet.begin());
    TPrecision af = 0;
    if ((ac[0] + ac[1]) > 0) af = (TPrecision) ac[1] / (TPrecision) (ac[0] + ac[1]);
    int32_t svlen = 1;
    if ((svt != NULL) && (std::string(svt) == "BND")) svlen = 0;
    else if (endsv != 0) svlen = endsv - rec->pos;
    if ((svt != NULL) && (std::string(svt) == "INS")) {
      if (inslen != NULL) svlen = *inslen;
    }
    int32_t ilen = 0;
    if ((precise) && (inslen != NULL)) ilen = *inslen;
    int32_t hlen = 0;
    if ((precise) && (homlen != NULL)) hlen = *homlen;
    TPrecision missingRate = (TPrecision) uncalled / (TPrecision) bcf_hdr_nsamples(hdr);
    
    TPrecision refratio = 0;
    _getMedian(ratioRef, refratio);
    TPrecision altratio = 0;
    _getMedian(ratioAlt, altratio);
    TPrecision maxaltratio = 0;
    _getMax(ratioAlt, maxaltratio);
    TPrecision refgq = 0;
    _getMedian(gqRef, refgq);
    TPrecision altgq = 0;
    _getMedian(gqAlt, altgq);
    TPrecision altRC = 0;
    _getMedian(rcAltRatio, altRC);
    TPrecision refRC = 0;
    _getMedian(rcRefRatio, refRC);
    TPrecision rdRatio = 0;
    if (refRC != 0) rdRatio = altRC / refRC;
    TPrecision rcMed = 0;
    _getMedian(rcRef, rcMed);

    // Write record
    for(TColumnHeader::const_iterator cHead = cHeader.begin(); cHead != cHeader.end(); ++cHead) {
      if (cHead != cHeader.begin()) std::cout << "\t";
      if (*cHead == "chr") std::cout << bcf_hdr_id2name(hdr, rec->rid);
      else if (*cHead == "start") std::cout << (rec->pos + 1);
      else if (*cHead == "chr2") std::cout << chr2Name;
      else if (*cHead == "end") std::cout << endsv;
      else if (*cHead == "id") std::cout << rec->d.id;
      else if (*cHead == "size") std::cout << svlen;
      else if (*cHead == "qual") std::cout << rec->qual;
      else if (*cHead == "vac") std::cout << ac[1];
      else if (*cHead == "vaf") std::cout << af;
      else if (*cHead == "pass") std::cout << passSite;
      else if (*cHead == "singleton") std::cout << rareCarrier;
      else if (*cHead == "missingrate") std::cout << missingRate;
      else if (*cHead == "svtype") {
	if (svtPresent) std::cout << svt;
	else std::cout << "NA";
      }
      else if (*cHead == "ct") std::cout << ctval;
      else if (*cHead == "precise") std::cout << precise;
      else if (*cHead == "ci") {
	if (ciPresent) std::cout << cipos[1];
	else std::cout << "NA";
      }
      else if (*cHead == "refratio") std::cout << refratio;
      else if (*cHead == "altratio") std::cout << altratio;
      else if (*cHead == "maxaltratio") std::cout << maxaltratio;
      else if (*cHead == "PEsupport") std::cout << totalPE;
      else if (*cHead == "SRsupport") std::cout << totalSR;
      else if (*cHead == "supportsum") std::cout << supportsum;
      else if (*cHead == "refgq") std::cout << refgq;
      else if (*cHead == "altgq") std::cout << altgq;
      else if (*cHead == "gqsum") std::cout << gqsum;
      else if (*cHead == "altgqsum") std::cout << altgqsum;
      else if (*cHead == "rdratio") std::cout << rdRatio;
      else if (*cHead == "medianrc") std::cout << rcMed;
      else if (*cHead == "inslen") std::cout << ilen;
      else if (*cHead == "homlen") std::cout << hlen;
      else if (*cHead == "fic") {
	if (ficPresent) std::cout << *fic;
	else std::cout << "NA";
      } else if (*cHead == "ce") {
	if ((precise) && (nce > 0)) std::cout << *ce;
	else std::cout << "0";
      }
      else if (*cHead == "exchet") {
	if (exchetPresent) std::cout << *exchet;
	else std::cout << "NA";
      } else if (*cHead == "hwepval") {
	if (hwePresent) std::cout << *hwepval;
	else std::cout << "NA";
      }
    }
    std::cout << std::endl;
  }

  // Clean-up
  if (svend != NULL) free(svend);
  if (inslen != NULL) free(inslen);
  if (pos2 != NULL) free(pos2);
  if (homlen != NULL) free(homlen);
  if (cipos != NULL) free(cipos);
  if (svt != NULL) free(svt);
  if (fic != NULL) free(fic);
  if (ce != NULL) free(ce);
  if (exchet != NULL) free(exchet);
  if (hwepval != NULL) free(hwepval);
  if (chr2 != NULL) free(chr2);
  if (ct != NULL) free(ct);
  if (gt != NULL) free(gt);
  if (rc != NULL) free(rc);
  if (rcl != NULL) free(rcl);
  if (rcr != NULL) free(rcr);
  if (dv != NULL) free(dv);
  if (dr != NULL) free(dr);
  if (rv != NULL) free(rv);
  if (rr != NULL) free(rr);
  if (gqInt != NULL) free(gqInt);
  if (gqFloat != NULL) free(gqFloat);
  bcf_hdr_destroy(hdr);
  bcf_close(ifile);
  bcf_destroy(rec);

  return 0;
}
