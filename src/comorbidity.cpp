// #include <unordered_set> // [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
//jw using namespace Rcpp;

typedef std::set<std::string> Icd9Set; // set of icd9 codes for comorbidities map
//typedef std::map<std::string, Icd9Set> CbdMapType; // comorbidity name to set of icd9 codes mapping
//typedef CbdMapType::iterator CbdMapTypeIt;
typedef std::pair<std::string, Icd9Set> CbdPairType;
typedef std::vector<CbdPairType> CbdPairVecType;
typedef CbdPairVecType::iterator CbdPairVecTypeIt;
typedef std::multimap<std::string, std::string> ViMapType; // visit:icd9 mapping - we have duplicate keys (visitids) // better as unordered multimap?
typedef ViMapType::iterator ViMapTypeIt;
typedef std::pair<ViMapTypeIt, ViMapTypeIt> ViPairItType;
typedef std::vector<std::string> StrVec;
typedef StrVec::iterator StrVecIt;

bool debug = false;

//' @name cvToVs
//' @title convert CharacterVector to C++ vector of strings
//' @description inverse of Rcpp::CharacterVector::create
//' @param cv vector of character strings
//' @return c++ vector of strings
//' @export
//' @import Rcpp
// [[Rcpp::export]]
std::vector<std::string> cvToVs(Rcpp::CharacterVector cv) {
  StrVec s(cv.size());
  for (int i=0; i<cv.size(); i++) {
    s[i] = std::string(cv[i]);
  }
  return(s);
}

// convert a mapping list to a map of 'sets'. Will enable fast look-up of an
// icd-9 code in each section of a mapping
//' @name icd9MappingToVectorSetsCpp
//' @title icd9MappingToVectorSetsCpp
//' @description purpose is to index each comorbidity so we can quickly tell
//' whether a given ICD-9 code is in that set.
//' @param icd9Mapping
//' @return C++ vector of sets of characters
//' @import Rcpp
std::vector<std::pair<std::string, std::set<std::string> > > icd9MappingToVectorSetsCpp(Rcpp::List icd9Mapping) {

  if (debug) { std::cout << "icd9MappingToVectorSetsCpp\n"; }

  CbdPairVecType cpv;
  Rcpp::StringVector svMapNames;
  Rcpp::CharacterVector oneMapOfIcd9s;

  int maplen = icd9Mapping.size();

  svMapNames = icd9Mapping.attr("names");

  StrVec cNames = cvToVs(svMapNames);

  for (StrVecIt cnm = cNames.begin(); cnm != cNames.end(); ++cnm) {
    if (debug) { std::cout << *cnm << ", "; }
    // get list for current comorbidity name
    oneMapOfIcd9s = icd9Mapping[*cnm];
    Icd9Set icd9s;
    StrVec svIcd9s = cvToVs(oneMapOfIcd9s); // convert to vec
    icd9s.insert(svIcd9s.begin(), svIcd9s.end()); // put all in set
    CbdPairType oneComorbidityPair = std::make_pair(*cnm, icd9s); // pair up name to icd9 codes for one comorbidity
    cpv.push_back(oneComorbidityPair); // put the pair in new c++ comorbidity map
  }
  if (debug) { std::cout << "\n"; }
  return(cpv);
}

//' @name icd9ComorbiditiesRaggedCpp
//' @title allocate co-morbidities based on list of vectors of ICD-9 codes
//' @param vipCodes list named by an Id (e.g. visit or patient), each with a character vector of short-form ICD-9 codes
//' @param icd9Mapping list named by  co-morbidity, with each item containing a character vector of ICD-9 codes for that co-morbidity
//' @return matrix of comorbidities, with row names being the visit or patients identifiers, column names being the comorbidities, and binary values.
//' @export
//' @import Rcpp
// [[Rcpp::export]]
Rcpp::LogicalMatrix icd9ComorbiditiesRaggedCpp(Rcpp::CharacterVector vipCodes, Rcpp::CharacterVector icd9Short, Rcpp::List icd9Mapping) {

  if (debug) { std::cout << "icd9ComorbiditiesRaggedCpp\n"; }

  Rcpp::CharacterVector comorbids = icd9Mapping.names();
  int nComorbids = comorbids.size();
  int nVipCodes = vipCodes.size();

  // todo setup hash map for each comorbidity then search may be quicker

  // initialize to zeroes
  Rcpp::LogicalMatrix out(nVipCodes, nComorbids);

  // loop through visits (i.e. row of output matrix)
  for (int v = 0; v < nVipCodes; v++) {
    SEXP li = vipCodes[v];
    Rcpp::CharacterVector icd9Codes(li);
    int nCodes = icd9Codes.size();
    // loop through each co-morbidity (i.e. column of output matrix)
    for (int c = 0; c < nComorbids; c++) {
      // now for each co-morbidity, e.g. HTN, check whether we have a match. No
      // need to continue checking once we have a single match. Could do this
      // with sugar 'match' but our lists are quite short, and I don't want all
      // the matches, just the first.
      Rcpp::CharacterVector comorbidSet = icd9Mapping[c];
      int nComorbidSet = comorbidSet.length();
      for (int d = 0; d < nComorbidSet; d++) {
        std::string comorbidToCheck = Rcpp::as<std::string>(comorbidSet[d]);
        // loop through icd-9 codes for current visit
        for (int i = 0; i < nCodes; i++) {
          std::string icd9CodeToCheck = Rcpp::as<std::string>(icd9Codes[i]);
          if (comorbidToCheck == icd9CodeToCheck) {
            out(v, c) = 1;
            //std::cout << "found a co-morbidity\n";
            goto doneThisComorbidity; // we have +ve so no need to look further
          }
        } // end loop through visit codes
      } // end loop through a comorbidity definition
      doneThisComorbidity: {}
    } // end loop through all comorbidities
  } // end loop through visits/patients

  // set the row and column names:
  Rcpp::List rcnames = Rcpp::List::create(vipCodes, comorbids);
  out.attr("dimnames") = rcnames;
  return(out);
}

//' @name icd9ComorbiditiesLongCpp
//' @title allocate co-morbidities based on vectors of IDs and ICD-9 codes
//' @param visitId character vector of visit or patient identifiers
//' @template icd9-short
//' @param icd9Mapping list named by  co-morbidity, with each item containing a character vector of ICD-9 codes for that co-morbidity
//' @return matrix of comorbidities, with row names being the visit or patients identifiers, column names being the comorbidities, and binary values.
//' @examples
//'  icd9ComorbiditiesLongCpp(c("pat1","pa2","three"), c("042","4011", "44179"), ahrqComorbid)
//' @export
//' @import Rcpp
// [[Rcpp::export]]
Rcpp::LogicalMatrix icd9ComorbiditiesLongCpp(Rcpp::CharacterVector visitId, Rcpp::CharacterVector icd9Short, Rcpp::List icd9Mapping) {
  if (debug) { std::cout << "icd9ComorbiditiesLongCpp\n"; }
  // convert my list of character vectors to a map of sets of strings, pure std/c++
  CbdPairVecType cMap;
  ViMapType viMap;
  Rcpp::CharacterVector uniqueVisitId;
  StrVec cbdNames;


  int nVisitId = visitId.size();
  // Rcpp::unique, then sort is horrible. what about a set? didn't do this before because I didn't want to deal with ordering a custom pair<string, string>
  // also note C++ unique is different to Rcpp::unique (the former only drops consecutive identical values)
  uniqueVisitId = Rcpp::unique(visitId);
  std::sort(uniqueVisitId.begin(), uniqueVisitId.end());

  int nUniqueVisitId = uniqueVisitId.size();
  int nMappingComorbids = icd9Mapping.size();

  // get C++ types for the input data: surely better way to do this.
  StrVec cppVisitId = cvToVs(visitId);
  StrVec cppIcd9Short = cvToVs(icd9Short);
  StrVec cppUniqueVisitId =  cvToVs(uniqueVisitId);

  cMap = icd9MappingToVectorSetsCpp(icd9Mapping);

  // turns out to be difficult to initialize maps from two vectors without boost or C++11:
  // TODO: consider using unordered_map C++11 or std::tr1 , but also available in std/tr1
  //std::unordered_map<Rcpp::CharacterVector> visitIdSet ( visitId );
  // this loop is probably expensive in time.
  if (debug) std::cout << "loop through visits to build vector of pairs for input data 'frame'\n";
  for (int ivi=0; ivi<nVisitId; ivi++) {
    if (debug) std::cout << ivi << ", ";
    viMap.insert(std::pair<std::string, std::string>(cppVisitId[ivi], cppIcd9Short[ivi]));
  }
  if (debug) { std::cout << "\n"; }
  // now we have a map to lookup icd9 codes for each unique visit Id. this is
  // likely to be more efficient with unordered list. Current assumption is that
  // the visitIds are to be returned sorted, which is reasonable, but should be
  // explicit (and the same as the R code).

  // straight to Rcpp matrix as our output. possible advantage in using native matrix first.
  Rcpp::LogicalMatrix out(nUniqueVisitId, nMappingComorbids); // rows,cols filled with zero
  //std::cout << "rows: " << out.nrow() << ", ";
  //std::cout << "cols: " << out.ncol() << "\n";
  // name the rows and cols up front, so we can reference the cols by name. This may be much slower...
  if (debug) { std::cout << "make list of the comorbidity names\n"; }
  for(CbdPairVecTypeIt it = cMap.begin(); it != cMap.end(); ++it) {
    if (debug) { std::cout << it->first << ", "; }
    cbdNames.push_back(it->first); // pair.second is the map.
  }
  if (debug) { std::cout << "\n"; }

  if (debug) { Rf_PrintValue(uniqueVisitId); }

// name the dimensions
  Rcpp::List rcnames = Rcpp::List::create(uniqueVisitId, cbdNames);
  out.attr("dimnames") = rcnames;
  std::string icd9CodeToCheck;

  // iterate through unique ids // should be able to avoid needing unique()
  if (debug) { std::cout << "main loop (through unique visit ids\n"; }

  for (StrVec::iterator vit=cppUniqueVisitId.begin(); vit!=cppUniqueVisitId.end(); ++vit) {
    if (debug) {std::cout << "visit id: " << *vit << "\n";}
    // now we can use the map structure to find the comorbidities, but not an iterator, since we need to make a matrix index... can't see a way to string index the matrix
    //for (CbdMapType::iterator it = cm.begin(); it != cm.end(); ++it) {

    // vit points to a unique visit id: now use 'cm' map to get list of icd-9 codes for that visit and loop through results:
    // *vit should be the (unique) visitId

    std::string icd9CodeToCheck;
    ViPairItType viPair;

    viPair = viMap.equal_range(*vit); // get pair of iterators that bracket the results (if any)

    // we then iterate through the range, with each iterator providing a multimap type
    for (ViMapTypeIt iit = viPair.first; iit != viPair.second; ++iit) {
      icd9CodeToCheck = iit->second; // equivalent to (*iit).second (first is the icd9 code in this case)
      if (debug) { std::cout << "icd code to check: " << icd9CodeToCheck << "\n"; }

      // have a single icd9 code for a single visit, now look it up in the comorbidities:
      //for (int cit = 0; cit < nMappingComorbids; cit++) { // use iterator and distance instead of numbers
      for (CbdPairVecTypeIt cit = cMap.begin(); cit!=cMap.end(); ++cit) {

        Icd9Set iSet = cit->second;
        if (iSet.find(icd9CodeToCheck) != iSet.end()) { // lookup the icd9 code in the icd9 'set' for current co-morbidity
        //std::cout << *vit << "\n";
        int vIndex = std::distance(cppUniqueVisitId.begin(), vit);
        int cIndex = std::distance(cMap.begin(), cit);
        if (debug) {
          std::cout << "v index = " << vIndex << ", ";
          std::cout << "c index = " << cIndex << "\n";
          //std::cout << icd9CodeToCheck << "\n";
        }
        out(vIndex, cIndex) = 1; // row x col
        } // end if matched
      } // end loop through mapping
    } // end loop through icd9 codes for a single visit
  } // end loop through visits
  return(out);
}

