// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <unordered_map>
#include <iostream>
#include <string.h>
#include <stdio.h>


using namespace Rcpp;
using namespace std;

// score across a stringset, each kmer with a custom score
// [[Rcpp::export]]
Rcpp::NumericVector countKmers(std::vector< std::string > strings, NumericVector score_val, std::vector< std::string > kmer_names, int k=10, int winsize=50)
{
  std::map<std::string, int> score;
  int n1 = score_val.size();
  for (int i = 0; i < n1; i++)
  {
    score[kmer_names[i]]=score_val[i];
  }

  int winNo= strings[0].length() / winsize;
  int kmerNo_win= winsize-k + 1;
  int n = strings.size();

  NumericMatrix out(n, winNo);

  for (int i = 0; i < n; i++) // each lig
  {
    cout << "ligand No." << i << "\n";
    for (int w=0; w< winNo; w++) // each win
    {
        for (int j=0; j<kmerNo_win; j++) // each pos of the win
        {
          std::string sub_s= strings[i].substr(winsize*w+j, k);
          if (  strchr(sub_s.c_str(), 'N') == NULL && strchr(sub_s.c_str(), 'n') == NULL )
          {    out(i,w) += score[sub_s];  }
        }
        out(i,w) /= kmerNo_win;
    }
  }

  return out;
}



/*** R
# countKmers finish
 */
