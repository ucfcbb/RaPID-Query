//  * --------------------------------------------------------------------------------------------------------
//  * Name: PBWTHeadExtension.cpp
//  * Description: Implementation of Richard Durbin's Positional Burrows-Wheeler Transform (PBWT) algorithms.
//  * Author: Yuan Wei 
//  * Created on: Jun 22, 2021
//  * --------------------------------------------------------------------------------------------------------

#include "PBWTHead.cpp"
using namespace std;

class PBWTHeadExtension: public PBWTHead{
    protected:
    vector<int> upointer;
    vector<int> vpointer;

    void buildPrefixDivergenceUpointerVpointer(){
        int d0t = site_index + 1;
        int d1t = site_index + 1;
        vector<int> prefix_1;
        vector<int> divergence_1;
        int u = 0;
        int v = 0;
        for (int i = 0; i < number_of_samples; i++){
            if (divergence_prev[i] > d0t){
                d0t = divergence_prev[i];
            }
            if (divergence_prev[i] > d1t){
                d1t = divergence_prev[i];
            }
            upointer.push_back(u);
            vpointer.push_back(v);
            if (sites[i] == 0){
                prefix.push_back(prefix_prev[i]);
                divergence.push_back(d0t);
                d0t = 0;
                u++;
            }
            else {
                prefix_1.push_back(prefix_prev[i]);
                divergence_1.push_back(d1t);
                d1t = 0;
                v++;
            }
        }
        prefix.insert(prefix.end(), prefix_1.begin(), prefix_1.end());
        divergence.insert(divergence.end(), divergence_1.begin(), divergence_1.end());
        for (int i = 0; i < vpointer.size(); i++){
            vpointer[i] = vpointer[i] + u;
        }
    }

    public:
    PBWTHeadExtension(): PBWTHead(){
    }

    PBWTHeadExtension(int site_index_in, vector<bool> sites_in): PBWTHead(){
        site_index = site_index_in;
        number_of_samples = sites_in.size();

        //this is the first site of samples (no previous prefix or divergence array is provided): take sites as currently given order and initialize prefix and divergence arrays
        sites = sites_in;
        for (int i = 0; i < number_of_samples; i++){
            prefix_prev.push_back(i);
            divergence_prev.push_back(0);
        }
        buildPrefixDivergenceUpointerVpointer();
    }

    PBWTHeadExtension(int site_index_in, vector<bool> sites_in, vector<int> prefix_prev_in, vector<int> divergence_prev_in): PBWTHead(){
        site_index = site_index_in;
        number_of_samples = sites_in.size();

        //this is not the first site of samples (previous prefix and divergence arrays are provided): take prefix and divergence arrays and assign sites based on their prefix values
        prefix_prev = prefix_prev_in;
        divergence_prev = divergence_prev_in;
        sites.reserve(number_of_samples);
        for (int i = 0; i < sites_in.size(); i++){
            sites[i] = sites_in[prefix_prev[i]];
        }
        buildPrefixDivergenceUpointerVpointer();
    }

    vector<int> getUpointer(){
        return upointer;
    }

    vector<int> getVpointer(){
        return vpointer;
    }
};
