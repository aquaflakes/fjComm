#include <Rcpp.h>
#include <unordered_map>
using namespace Rcpp;



//' 不取反向互补
 //' 可支持gapped kmer的比对
 //' kmer个数不限
 // [[Rcpp::export]]
 SEXP score_genome_gkmers_cpp(std::vector<std::string> seqnames, std::vector<std::string> seqs, IntegerVector starts, std::vector<std::string> target_kmers, NumericVector kmer_scores, int gapmin=0, int gapmax=8)
 {
   Function Rle("Rle", Environment::namespace_env("IRanges"));

   std::unordered_map<std::string, double> kmer_score_map;
   for (int i = 0; i < target_kmers.size(); i++) {
     kmer_score_map[target_kmers[i]]= kmer_scores[i];
   }
   int k = target_kmers[0].length()-2;
   int k_half=k/2;

   List results;
   for (int j = 0; j < seqs.size(); j++)
   {
     int seq_length = seqs[j].length();
     int start_pos = starts[j];
     std::vector<double> currSeq_scores(seq_length,0);
     // std::cout << "unordered map" << seqs.size() << std::endl;

     for (int i = 0; i <= seq_length - k; i++)
     {
       for(int gap=gapmin; gap<=gapmax; gap++)
       {
         if(i+k+gap > seq_length){break;}
         // std::cout << seqs[j].substr(i, k_half)+std::to_string(gap)+"n"+seqs[j].substr(i+k_half+gap, k_half) << "\n";
         std::string curr_k=seqs[j].substr(i, k_half)+std::to_string(gap)+"n"+seqs[j].substr(i+k_half+gap, k_half);
         double curr_score=kmer_score_map[curr_k];
         if(curr_score){
           currSeq_scores[(i+floor((k+gap)/2))]+=curr_score;
         }
       }
     }
     results[seqnames[j]]=Rle(currSeq_scores);
   }

   if (seqs.size()>1){
     Function RleList("RleList", Environment::namespace_env("IRanges"));
     return RleList(results);
   }else{
     return results[0];
   }
 }

 //
 // /*** R
 // # Example usage
 // score_genome_with_iupac_kmers(seqnames = c("Chr1","Chr2","Chr3"),seqs = c("CAAAAAAAACCAC","CCCCTCCCCATT","TTTTTTTTTT"),starts = c(1,1,1),
 //                               target_kmers = c("AAAA0nAAAA","CCCC1nCCCC"),kmer_scores = c(1,4),gapmin = 0,gapmax = 8)
 // */
