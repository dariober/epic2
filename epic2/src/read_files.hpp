#include <iostream>
#include <algorithm>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <vector>
#include <iterator>
// #include <cstdint>
// #include <unordered_map>
// #include <zlib.h>



typedef std::vector<uint32_t> intvec;
typedef std::pair <std::string, char> key;
typedef std::map<key, intvec> genome_map;

struct interval { uint32_t start; uint32_t end; };
typedef std::vector<interval> interval_vector;
// typedef std::map<key, intvec> genome_tags;
typedef std::map<key, interval_vector> genome_intervals;

uint32_t compare_by_start_end(const interval lhs, const interval rhs){
  if (lhs.start < rhs.start){
    return 1;
  } else if (rhs.start < lhs.start){
    return 0;
  } else if (lhs.end < rhs.end){
    return 1;
  } else {
    return 0;
      };

}


uint32_t start_end_equal(const interval lhs, const interval rhs){
  if ((lhs.start == rhs.start) && (lhs.end == rhs.end)){
      return 1;
    } else {
    return 0;
  }
}


genome_map intervals_to_five_end(genome_intervals genome, uint32_t drop_duplicates) {

  genome_intervals::iterator it;
  char strand;

  uint i = 0;

  key chrom_strand;
  interval_vector intervals;
  genome_map genome_tags;

  it = genome.begin();

  while(it != genome.end()){

    chrom_strand = it->first;
    strand = it->first.second;
    intervals = it->second;

    if (drop_duplicates){
      std::sort(intervals.begin(), intervals.end(), compare_by_start_end);
      intervals.erase(unique(intervals.begin(), intervals.end(), start_end_equal), intervals.end());
    }

    i = 0;
    if (strand == '+') {
      for (i = 0; i < intervals.size(); i++){
        genome_tags[chrom_strand].push_back(intervals[i].start);
      }
    } else {
      for (i = 0; i < intervals.size(); i++){
        genome_tags[chrom_strand].push_back(intervals[i].end);
      }
    }
    genome.erase(chrom_strand);
    it++;

  }

  return genome_tags;
}



genome_map read_bed(char const* fileName, uint32_t drop_duplicates)
{
  std::ifstream file(fileName);

  std::string   chromosome;
  std::string   junk;
  uint32_t start;
  uint32_t end;
  char strand;


  key chrom_strand;
  interval _interval;
  // intvec tags;
  interval_vector intervals;
  genome_intervals genome;


  while(file >> chromosome >> start >> end >> junk >> junk >> strand)
    {

      chrom_strand = std::make_pair(chromosome, strand);

      _interval = {start, end - 1};

      genome[chrom_strand].push_back(_interval);

    }

  return intervals_to_five_end(genome, drop_duplicates);
}





genome_map read_bedpe(char const* fileName, uint32_t drop_duplicates)
{
  std::ifstream file(fileName);



  std::string   chromosome;
  std::string   junk;
  uint32_t start;
  uint32_t end;
  char strand;


  key chrom_strand;
  interval _interval;
  // intvec tags;
  interval_vector intervals;
  genome_intervals genome;


  while(file >> chromosome >> start >> junk >> junk >> junk >> end >> junk >> junk >> strand >> junk)
    {

      chrom_strand = std::make_pair(chromosome, strand);

      _interval = {start, end - 1};

      genome[chrom_strand].push_back(_interval);

    }

  return intervals_to_five_end(genome, drop_duplicates);

}

#include "gzstream.h"
// using namespace gz;
genome_map read_bed_gz(char const* fileName, uint32_t drop_duplicates)
{
  igzstream in(fileName);


  std::string   chromosome;
  std::string   junk;
  uint32_t start;
  uint32_t end;
  char strand;


  key chrom_strand;
  interval _interval;
  // intvec tags;
  interval_vector intervals;
  genome_intervals genome;

  while(in >> chromosome >> start >> end >> junk >> junk >> strand)
    {

      chrom_strand = std::make_pair(chromosome, strand);

      _interval = {start, end - 1};

      genome[chrom_strand].push_back(_interval);

    }

  return intervals_to_five_end(genome, drop_duplicates);

}


genome_map read_bedpe_gz(char const* fileName, uint32_t drop_duplicates)
{

  igzstream in(fileName);


  std::string   chromosome;
  std::string   junk;
  uint32_t start;
  uint32_t end;
  char strand;


  key chrom_strand;
  interval _interval;
  // intvec tags;
  interval_vector intervals;
  genome_intervals genome;

  while(in >> chromosome >> start >> junk >> junk >> junk >> end >> junk >> junk >> strand >> junk)
    {

      chrom_strand = std::make_pair(chromosome, strand);

      _interval = {start, end - 1};

      genome[chrom_strand].push_back(_interval);

    }

  return intervals_to_five_end(genome, drop_duplicates);
}

// int main(){

//   read_bedpe_gz("examples/test.bedpe.gz", true);

//   return 0;
// }
