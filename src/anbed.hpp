#ifndef BED_HPP
#define BED_HPP

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

class BED {
  public:
    std::string chr;
    uint32_t start;
    uint32_t end;
    BED(std::string, uint32_t, uint32_t);
    BED(std::string, uint32_t);
    BED(std::string&);
    int isOverlap(const BED&, const uint32_t&) const;
    void merge(const BED&);
    uint32_t size() const;
    int isEqual(const BED&) const;
    std::string toString() const;
    std::string toBEDstring() const;
};

inline BED::BED(
 std::string c, 
  uint32_t s, 
  uint32_t e): chr(c), start(s), end(e){};

inline BED::BED(std::string& line)
{
  if(line.empty()){
    fprintf(stderr,"Cannot parse empty GFF line\n");
    exit(1);
  } 
  else {
      int index = 0;
      std::string value;
      std::istringstream stream(line);
      while(std::getline(stream, value, ':')) {
        if(index == 0) chr = value;
        else if(index == 1){
          int index2 = 0;
          std::string value2;
          std::istringstream stream2(value);
          while(std::getline(stream2, value2, '-')) {
            if(index2 == 0) start = static_cast<uint32_t>(std::stoul(value2));
            else if(index2 == 1) end = static_cast<uint32_t>(std::stoul(value2));
            ++index2;
          }
        }
        ++index;
      }
  }
}

inline BED::BED(const std::string line, uint32_t flanksize){
  if(line.empty()){
    fprintf(stderr,"Cannot parse empty GFF line\n");
    exit(1);
  } else {
      int index = 0;
      std::string value;
      std::istringstream stream(line);
      while(std::getline(stream, value, '\t')){
        switch(index){
          case 0: {
            chr = value;
            break;
          }
          case 1: {
            start = static_cast<uint32_t>(std::stoul(value)) - flanksize;
            break;
          }
          case 2: {
            end = static_cast<uint32_t>(std::stoul(value)) + flanksize;
            break;
          }
          default: {
            break;
          }
        }
        ++index;
      }       
  }
}

inline int BED::isOverlap(const BED& that, const uint32_t& dist) const
{
  return (chr == that.chr) && (start <= that.end + dist)  &&  (end >= that.start - dist);
}

inline uint32_t BED::size() const
{
  return end - start;
}


inline void BED::merge(const BED& that)
{
  start = start < that.start ? start : that.start;
  end = end > that.end ? end : that.end;
}


inline int BED::isEqual( const BED& that ) const
{
    return chr == that.chr &&
    start == that.start && 
    end == that.end;
}

inline std::string BED::toString() const
{
  return chr + "\t" + std::to_string(start) + "\t" + std::to_string(end);
}

inline std::string BED::toBEDstring() const
{
  return chr + ":" + std::to_string(start) + "-" + std::to_string(end);
}

//type alias for map -> (chrm index, vector of Region objects)
typedef std::map<std::string, std::vector<BED>> BEDmap;


inline int sortbed(
  const BED &a,
  const BED &b)
{

  return a.start < b.start;
}

inline int sortbedsize(
  const BED &a,
  const BED &b)
{

  if(a.chr == b.chr){
    return a.size() > b.size();
  } else {
    return a.chr < b.chr;
  }
}

inline void mergebeds(
    std::vector<BED>& beds,
    std::vector<BED>& acc_beds,
    const uint32_t& dist
  )
{
  auto it = beds.begin();
  BED current_bed = *it;
  ++it;
  while(it != beds.end()){
    if(current_bed.isOverlap(*it, dist)) current_bed.merge(*it);
    else{
      acc_beds.push_back(current_bed);
      current_bed = *it;
    }
    ++it;
  }

  acc_beds.push_back(current_bed);
}

inline void uniquebed(
    std::vector<BED>& beds
  )
{
  auto it = beds.begin() + 1;
  while(it != beds.end()){
    if(!(it - 1)->isEqual(*it)) ++it;
    else it = beds.erase(it);
  }
}


inline void bedparser(
  std::string bedfile, 
  BEDmap& bedmap,
  uint32_t flanksize)
{
  std::ifstream inputbed(bedfile.c_str());
  std::string line;
  // Read the next line from File untill it reaches the end.
  while (std::getline(inputbed, line)){
    BED annotation = BED(line, flanksize);
    //check if chromosome already exists in index
      auto it = bedmap.find(annotation.chr);
      //exists, update key-value 
      if (it != bedmap.end()) {
        it->second.push_back(annotation);
      }
      //does not exist, create key,value entry
      else {
        //set temp vector of coords
        std::vector<BED> localbeds;
        //add current coord to temp vector
        localbeds.push_back(annotation);
        //add vector
        bedmap[annotation.chr] = localbeds;
      }
  }
  inputbed.close();
}

#endif