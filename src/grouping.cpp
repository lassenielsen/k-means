#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <list>
#include <string>
#include <fstream>
#include <streambuf>
#include <algorithm>
#include <cstdlib>


using namespace std;

template<class T> void k_means(const vector<T> &points, size_t k, double (*distance)(const T&, const T&), T (*mean)(const vector<T>&), vector<T> &centers, map<size_t,vector<T> > &groups) // {{{
{ if (k==0)
    throw string("k-means error: k must be greater than zero");
  // Randomly initialize k centers;
  centers.clear();
  for (size_t point=0; centers.size()<k; ++point)
  { if (point>points.size())
      throw string("k-means error: Not enough distinct points for creating k groups");
    bool distinct=true;
    for (size_t group=0; group<centers.size(); ++group)
    { if (distance(points[point],centers[group])==0)
        distinct=false;
    }
    if (distinct)
      centers.push_back(points[point]);
  }
  
  // Perform grouping
  double ref_err=(size_t)-1;
  for (;;) // Repeat until fixpoint
  { groups.clear();
    double sum_err=0;
    // Group with current centers
    for (size_t point=0; point<points.size(); ++point)
    { double min_dist=-1.0;
      size_t closest_group=0;
      for (size_t group=0; group<k; ++group)
      { double dist=distance(points[point],centers[group]);
        if (dist<min_dist || min_dist<0)
        { closest_group=group;
          min_dist=dist;
        }
      }
      groups[closest_group].push_back(points[point]);
      sum_err+=min_dist*min_dist;
    }
    cout << "sum_err: " << sum_err << endl;
    // Test if fixpoint
    if (ref_err<=sum_err && ref_err!=(size_t)-1)
      return;
    else
      ref_err=sum_err;

    // Recalculate centers
    for (size_t group=0; group<groups.size(); ++group)
    { centers[group]=mean(groups[group]);
    }
  }
} // }}}

double sqr(double x) // {{{
{ return x*x;
} // }}}
double pair_distance(const pair<double,double> &lhs, const pair<double,double> &rhs) // {{{
{ return sqrt(sqr(rhs.first-lhs.first)+sqr(rhs.second-lhs.second));
} // }}}
pair<double,double> pair_mean(const vector<pair<double,double> > &vals) // {{{
{ if (vals.empty())
    return pair<double,double>(0.0,0.0);
  pair<double,double> sum(0.0,0.0);
  for (auto it=vals.begin(); it!=vals.end(); ++it)
  { sum.first+=it->first;
    sum.second+=it->second;
  }
  return pair<double,double>(sum.first/vals.size(),sum.second/vals.size());
} // }}}
inline double min(const double &lhs, const double &rhs) // {{{
{ return lhs<rhs?lhs:rhs;
} // }}}
double lcs_distance(const string &lhs, const string &rhs) // {{{
{ double *dists=new double[rhs.size()+1];
  double *prev_dists=new double[rhs.size()+1];
  // Initialize prev_dists
  for (size_t i=0; i<rhs.size()+1; ++i)
    dists[i]=(double)i;

  for (size_t lhspos=0; lhspos<lhs.size(); ++lhspos)
  { double *tmp_dists=dists;
    dists=prev_dists;
    prev_dists=tmp_dists;
    // Calc dists
    dists[0]=prev_dists[0]+1.0;
    for (size_t rhspos=0; rhspos<rhs.size(); ++rhspos)
    { if (lhs[lhspos]==rhs[rhspos])
        dists[rhspos+1]=min(dists[rhspos]+1.0,min(prev_dists[rhspos+1]+1.0,prev_dists[rhspos]));
      else
        dists[rhspos+1]=min(dists[rhspos]+1.0,prev_dists[rhspos+1]+1.0);
    }
  }
  double result=dists[rhs.size()];
  delete [] dists;
  delete [] prev_dists;
  return result;
} // }}}
string lcs_mean(const vector<string> &vals) // {{{
{ // First compute shortest common supersequence with usage info
  list<char> scs_symbols;
  list<size_t> scs_usage;
  for (auto val=vals.begin(); val!=vals.end(); ++val)
  { // Do lcs table
    size_t scs_size=scs_symbols.size();
    double *dists=new double[(val->size()+1)*(scs_size+1)];
    size_t *prev=new size_t[(val->size()+1)*(scs_size+1)];
    // Initialize first row
    for (size_t i=0; i<scs_size+1; ++i)
    { dists[i]=(double)i;
      prev[i]=0; // meaning horizontal
    }
    for (size_t valpos=0; valpos<val->size(); ++valpos)
    { size_t scs_pos=0;
      dists[(scs_size+1)*(valpos+1)]=dists[(scs_size+1)*valpos]+1.0;
      prev[(scs_size+1)*(valpos+1)]=2; // meaning vertical
      for (auto scs_it=scs_symbols.begin(); scs_it!=scs_symbols.end(); ++scs_it, ++scs_pos)
      { if ((*val)[valpos]==*scs_it &&
            dists[scs_pos+(scs_size+1)*valpos]<=dists[scs_pos+(scs_size+1)*(valpos+1)]+1.0 &&
            dists[scs_pos+(scs_size+1)*valpos]<dists[scs_pos+1+(scs_size+1)*valpos]+1.0)
        { dists[scs_pos+1+(scs_size+1)*(valpos+1)]=dists[scs_pos+(scs_size+1)*valpos];
          prev[scs_pos+1+(scs_size+1)*(valpos+1)]=1; // meaning diagonal
        }
        else if (dists[scs_pos+(scs_size+1)*(valpos+1)]<=dists[scs_pos+1+(scs_size+1)*valpos])
        { dists[scs_pos+1+(scs_size+1)*(valpos+1)]=dists[scs_pos+(scs_size+1)*(valpos+1)]+1.0;
          prev[scs_pos+1+(scs_size+1)*(valpos+1)]=0; // meaning horizontal
        }
        else
        { dists[scs_pos+1+(scs_size+1)*(valpos+1)]=dists[scs_pos+1+(scs_size+1)*valpos]+1.0;
          prev[scs_pos+1+(scs_size+1)*(valpos+1)]=2; // meaning vertical
        }
      }
    }
    delete [] dists;
    // Update scs from prev
    auto symbols_it=scs_symbols.end();
    auto usage_it=scs_usage.end();
    size_t x=scs_size,y=val->size();
    //cout << "!!!" << endl;
    //cout << x << endl;
    //cout << y << endl;
    while (x!=0 || y!=0)
    { size_t dir=prev[x+(scs_size+1)*y];
      //cout << "dir=" << dir << endl;
      if (dir==0) // horizontal
      { --x;
        --symbols_it;
        --usage_it;
        // No changes to scs symbols or usage
      }
      else if (dir==1) //  diagonal
      { --x;
        --y;
        --symbols_it;
        --usage_it;
        // Update usage
        ++(*usage_it);
      }
      else // dir==2 vertical
      { --y;
        // Add char to scs symbols with one usage
        //cout << "Adding symbol: " << (*val)[y] << endl;
        symbols_it=scs_symbols.insert(symbols_it,(*val)[y]);
        usage_it=scs_usage.insert(usage_it,1);
      }
    }
    delete [] prev;
  }
  //cout << "SCS: ";
  //{ auto symbol_it=scs_symbols.begin();
  //  for (auto usage_it=scs_usage.begin(); usage_it!=scs_usage.end(); ++usage_it, ++symbol_it)
  //    cout << "(" << *symbol_it << "," << *usage_it << ")";
  //}
  //cout << endl;
  size_t meanlength=0;
  for (auto usage_it=scs_usage.begin(); usage_it!=scs_usage.end(); ++usage_it)
    if (*usage_it>vals.size()/2)
      ++meanlength;

  string result(meanlength,'.');
  size_t respos=0;
  auto symbol_it=scs_symbols.begin();
  for (auto usage_it=scs_usage.begin(); usage_it!=scs_usage.end(); ++usage_it, ++symbol_it)
    if (*usage_it>vals.size()/2)
      result[respos++]=*symbol_it;

  return result;
} // }}}

inline string readfile(string filename) // {{{
{ ifstream t(filename);
  string s=std::string((std::istreambuf_iterator<char>(t)),
                 std::istreambuf_iterator<char>()).substr(0,100);
  while (s.size()<100)
    s+=char(rand()%256);
  return s;
} // }}}
inline string findfile(const vector<string> &filenames, const vector<string> &strings, const string &s) // {{{
{ auto fit=filenames.begin();
  auto sit=strings.begin();
  while (fit!=filenames.end() && sit!=strings.end())
  { if (*sit==s)
      return *fit;
    ++fit;
    ++sit;
  }
  return "File not found";
} // }}}

int main(int argc, char** argv) // {{{
{ cout << "Start" << endl;
  cout << "Testing coordinate grouping" << endl;
  vector<pair<double,double> > points;
  points.push_back(pair<double,double>(-2.0,-2.0));
  points.push_back(pair<double,double>(-1.0,-1.0));
  points.push_back(pair<double,double>( 2.0,-2.0));
  points.push_back(pair<double,double>( 1.0,-1.0));
  points.push_back(pair<double,double>(-2.0, 2.0));
  points.push_back(pair<double,double>(-1.0, 1.0));
  points.push_back(pair<double,double>( 2.0, 2.0));
  points.push_back(pair<double,double>( 1.0, 1.0));
  vector<pair<double,double> > centers;
  map<size_t,vector<pair<double,double> > > groups;
  k_means(points,4,&pair_distance,&pair_mean,centers,groups);
  cout << "Results" << endl;
  for (size_t i=0; i<centers.size(); ++i)
  { cout << "Group center: " << centers[i].first << "," << centers[i].second << endl;
    for (auto point=groups[i].begin(); point!=groups[i].end(); ++point)
      cout << "  point: " << point->first << "," << point->second << endl;
  }
  cout << "Testing lcs" << endl;
  cout << "Distance between 'abcabc' and 'ababab' is " << lcs_distance("abcabc","ababab") << endl;
  cout << "Distance between 'abcabc' and 'cbabca' is " << lcs_distance("abcabc","cbacba") << endl;
  cout << "Distance between 'abcabc' and 'distant' is " << lcs_distance("abcabc","distant") << endl;
  vector<string> strs1;
  strs1.push_back("abcabc");
  strs1.push_back("cbacba");
  cout << "Subsequence for 'abcabc' and 'cbacba' is " << lcs_mean(strs1) << endl;
  vector<string> strs2;
  strs2.push_back("abcabc");
  strs2.push_back("ababab");
  cout << "Subsequence for 'abcabc' and 'ababab' is " << lcs_mean(strs2) << endl;
  vector<string> strs3;
  strs3.push_back("abcabc");
  strs3.push_back("ababab");
  strs3.push_back("cbacba");
  cout << "Subsequence is " << lcs_mean(strs3) << endl;

  vector<string> filenames;
  filenames.push_back("data/python1.txt");
  filenames.push_back("data/python2.txt");
  filenames.push_back("data/python3.txt");
  filenames.push_back("data/python4.txt");
  filenames.push_back("data/c1.txt");
  filenames.push_back("data/c2.txt");
  filenames.push_back("data/c3.txt");
  filenames.push_back("data/c4.txt");
  filenames.push_back("data/java1.txt");
  filenames.push_back("data/java2.txt");
  filenames.push_back("data/java3.txt");
  filenames.push_back("data/java4.txt");
  filenames.push_back("data/ocaml1.txt");
  filenames.push_back("data/ocaml2.txt");
  filenames.push_back("data/ocaml3.txt");
  filenames.push_back("data/ocaml4.txt");
  random_shuffle(filenames.begin(),filenames.end());
  vector<string> strings;
  for (auto filename=filenames.begin(); filename!=filenames.end(); ++filename)
    strings.push_back(readfile(*filename));
  vector<string> centers2;
  map<size_t,vector<string> > groups2;
  k_means(strings,4,&lcs_distance,&lcs_mean,centers2,groups2);
  cout << "Results" << endl;
  for (size_t i=0; i<centers2.size(); ++i)
  { cout << "Group: " << endl;
    for (auto point=groups2[i].begin(); point!=groups2[i].end(); ++point)
      cout << "  file: " << findfile(filenames,strings,*point) << endl;
  }
  return 0;
} // }}}
