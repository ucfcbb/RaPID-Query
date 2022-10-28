//  * --------------------------------------------------------------------------------------------------------
//  * Name: PBWTHead.cpp
//  * Description: Implementation of Richard Durbin's Positional Burrows-Wheeler Transform (PBWT) algorithms.
//  * Author: Yuan Wei 
//  * Created on: Jun 22, 2021
//  * --------------------------------------------------------------------------------------------------------

#include <vector>
#include <tuple>
using namespace std;

class PBWTHead{
    protected:
    int number_of_samples;
    int site_index;
    vector<bool> sites;
    vector<int> prefix;
    vector<int> divergence;
    vector<int> prefix_prev;
    vector<int> divergence_prev;

    void buildPrefixDivergence(){
        int d0t = site_index + 1;
        int d1t = site_index + 1;
        vector<int> prefix_1;
        vector<int> divergence_1;
        for (int i = 0; i < number_of_samples; i++){
            if (divergence_prev[i] > d0t){
                d0t = divergence_prev[i];
            }
            if (divergence_prev[i] > d1t){
                d1t = divergence_prev[i];
            }
            if (sites[i] == 0){
                prefix.push_back(prefix_prev[i]);
                divergence.push_back(d0t);
                d0t = 0;
            }
            else {
                prefix_1.push_back(prefix_prev[i]);
                divergence_1.push_back(d1t);
                d1t = 0;
            }
        }
        prefix.insert(prefix.end(), prefix_1.begin(), prefix_1.end());
        divergence.insert(divergence.end(), divergence_1.begin(), divergence_1.end());
    }

    public:
    PBWTHead(){
    }

    PBWTHead(int site_index_in, vector<bool> sites_in){
        site_index = site_index_in;
        number_of_samples = sites_in.size();

        //this is the first site of samples (no previous prefix and divergence array is provided): take sites as currently given order and initialize prefix and divergence arrays
        sites = sites_in;
        for (int i = 0; i < number_of_samples; i++){
            prefix_prev.push_back(i);
            divergence_prev.push_back(0);
        }
        buildPrefixDivergence();
    }

    PBWTHead(int site_index_in, vector<bool> sites_in, vector<int> prefix_prev_in, vector<int> divergence_prev_in){
        site_index = site_index_in;
        number_of_samples = sites_in.size();

        //this is not the first site of samples (previous prefix and divergence arrays are provided): take prefix and divergence arrays and assign sites based on their prefix values
        prefix_prev = prefix_prev_in;
        divergence_prev = divergence_prev_in;
        sites.reserve(number_of_samples);
        for (int i = 0; i < sites_in.size(); i++){
            sites[i] = sites_in[prefix_prev[i]];
        }
        buildPrefixDivergence();
    }
    
    vector<int> getPrefix(){
        return prefix;
    }

    vector<int> getDivergence(){
        return divergence;
    }

    vector<tuple<int, int, int, int>> getLongMatches(int l, bool is_last_site){
        //matched pair is reported as a tuple: (sample_index_1, sample_index_2, starting_position_(inclusive), ending_position_(inclusive)), where sample_index_1 <= sample_index_2
        vector<tuple<int, int, int, int>> long_matches;
        if (l > 0){
            vector<int> list_0;
            vector<int> list_1;
            for (int i = 0; i < number_of_samples; i++){
                if (divergence_prev[i] > site_index - l){
                    if (list_0.size() > 0 && list_1.size() > 0){
                        int is = min(list_0[0], list_1[0]);
                        int ie = max(list_0[list_0.size() - 1], list_1[list_1.size() - 1]);
                        for (int i1 = is; i1 <= ie; i1++){
                            int dst = 0;
                            for (int i2 = i1 + 1; i2 <= ie; i2++){
                                if (divergence_prev[i2] > dst){
                                    dst = divergence_prev[i2];
                                }
                                if (sites[i1] != sites[i2]){
                                    if (prefix_prev[i1] <= prefix_prev[i2]){
                                        long_matches.push_back(tuple<int, int, int, int>(prefix_prev[i1], prefix_prev[i2], dst, site_index - 1));
                                    }
                                    else {
                                        long_matches.push_back(tuple<int, int, int, int>(prefix_prev[i2], prefix_prev[i1], dst, site_index - 1));
                                    }
                                }
                            }
                        }
                    }
                    list_0.clear();
                    list_1.clear();
                }
                if (sites[i] == 0){
                    list_0.push_back(i);
                }
                else {
                    list_1.push_back(i);
                }
            }

            //get matched pairs for the last sample case
            if (list_0.size() > 0 && list_1.size() > 0){
                int is = min(list_0[0], list_1[0]);
                int ie = max(list_0[list_0.size() - 1], list_1[list_1.size() - 1]);
                for (int i1 = is; i1 <= ie; i1++){
                    int dst = 0;
                    for (int i2 = i1 + 1; i2 <= ie; i2++){
                        if (divergence_prev[i2] > dst){
                            dst = divergence_prev[i2];
                        }
                        if (sites[i1] != sites[i2]){
                            if (prefix_prev[i1] <= prefix_prev[i2]){
                                long_matches.push_back(tuple<int, int, int, int>(prefix_prev[i1], prefix_prev[i2], dst, site_index - 1));
                            }
                            else {
                                long_matches.push_back(tuple<int, int, int, int>(prefix_prev[i2], prefix_prev[i1], dst, site_index - 1));
                            }
                        }
                    }
                }
            }
            list_0.clear();
            list_1.clear();

            //get matched pairs if this is the last site of the panel
            if (is_last_site){
                vector<int> list_t;
                for (int i = 0; i < number_of_samples; i++){
                    if (divergence[i] > site_index - l + 1){
                        if (list_t.size() > 1){
                            for (int i1 = 0; i1 < list_t.size(); i1++){
                                int dst = 0;
                                for (int i2 = i1 + 1; i2 < list_t.size(); i2++){
                                    if (divergence[list_t[i2]] > dst){
                                        dst = divergence[list_t[i2]];
                                    }
                                    if (prefix[list_t[i1]] <= prefix[list_t[i2]]){
                                        long_matches.push_back(tuple<int, int, int, int>(prefix[list_t[i1]], prefix[list_t[i2]], dst, site_index));
                                    }
                                    else {
                                        long_matches.push_back(tuple<int, int, int, int>(prefix[list_t[i2]], prefix[list_t[i1]], dst, site_index));
                                    }
                                }
                            }
                        }
                        list_t.clear();
                    }
                    list_t.push_back(i);
                }

                //report matched sample pair for the last sample case
                if (list_t.size() > 1){
                    for (int i1 = 0; i1 < list_t.size(); i1++){
                        int dst = 0;
                        for (int i2 = i1 + 1; i2 < list_t.size(); i2++){
                            if (divergence[list_t[i2]] > dst){
                                dst = divergence[list_t[i2]];
                            }
                            if (prefix[list_t[i1]] <= prefix[list_t[i2]]){
                                long_matches.push_back(tuple<int, int, int, int>(prefix[list_t[i1]], prefix[list_t[i2]], dst, site_index));
                            }
                            else {
                                long_matches.push_back(tuple<int, int, int, int>(prefix[list_t[i2]], prefix[list_t[i1]], dst, site_index));
                            }
                        }
                    }
                }
                list_t.clear();
            }
        }
        return long_matches;
    }

    vector<tuple<int, int, int, int>> getSetMaximalMatches(bool is_last_site){
        //matched pair is reported as a tuple: (sample_index_current, sample_index_other, starting_position_(inclusive), ending_position_(inclusive))
        vector<tuple<int, int, int, int>> set_maximal_matches;
        int i1 = -1;
        int i2 = number_of_samples;
        bool searchNext = false;
        bool searchedPrevious = false;
        bool searchedFollowing = false;
        for (int i = 0; i < number_of_samples; i++){
            i1 = i - 1;
            i2 = i + 1;
            searchNext = false;
            searchedPrevious = false;
            searchedFollowing = false;

            //search for all previous samples
            if (((i + 1 <= number_of_samples - 1 && divergence_prev[i] <= divergence_prev[i + 1]) || i == number_of_samples - 1) && i1 >= 0){
                while (i1 >= 0 && divergence_prev[i1 + 1] <= divergence_prev[i] && (divergence_prev[i] <= site_index - 1 && divergence_prev[i1 + 1] <= site_index - 1)){
                    searchedPrevious = true;
                    if (sites[i1] == sites[i]){
                        searchNext = true;
                        break;
                    }
                    else {
                        i1--;
                    }
                }
                if (searchNext){
                    continue;
                }
            }

            //search for all following samples
            if (i + 1 <= number_of_samples - 1 && divergence_prev[i + 1] <= divergence_prev[i] && i2 <= number_of_samples - 1){
                while (i2 <= number_of_samples - 1 && divergence_prev[i2] <= divergence_prev[i + 1] && (divergence_prev[i + 1] <= site_index - 1 && divergence_prev[i2] <= site_index - 1)){
                    searchedFollowing = true;
                    if (sites[i2] == sites[i]){
                        searchNext = true;
                        break;
                    }
                    else {
                        i2++;
                    }
                }
                if (searchNext){
                    continue;
                }
            }
            
            //report set maximal matches for current sample
            if (searchedPrevious && i1 + 1 >= 0){
                for (int j = i1 + 1; j <= i - 1; j++){
                    set_maximal_matches.push_back(tuple<int, int, int, int>(prefix_prev[i], prefix_prev[j], divergence_prev[i], site_index - 1));
                }
            }
            if (searchedFollowing && i2 - 1 <= number_of_samples - 1){
                for (int j = i + 1; j <= i2 - 1; j++){
                    set_maximal_matches.push_back(tuple<int, int, int, int>(prefix_prev[i], prefix_prev[j], divergence_prev[i + 1], site_index - 1));
                }
            }
        }

        //get set maximal matches if this is the last site of the panel
        if (is_last_site){
            for (int i = 0; i < number_of_samples; i++){
                i1 = i - 1;
                i2 = i + 1;
                searchedPrevious = false;
                searchedFollowing = false;

                //search for all previous samples
                if (((i + 1 <= number_of_samples - 1 && divergence[i] <= divergence[i + 1]) || i == number_of_samples - 1) && i1 >= 0){
                    while (i1 >= 0 && divergence[i1 + 1] <= divergence[i] && (divergence[i] <= site_index && divergence[i1 + 1] <= site_index)){
                        searchedPrevious = true;
                        i1--;
                    }
                }

                //search for all following samples
                if (i + 1 <= number_of_samples - 1 && divergence[i + 1] <= divergence[i] && i2 <= number_of_samples - 1){
                    while (i2 <= number_of_samples - 1 && divergence[i2] <= divergence[i + 1] && (divergence[i + 1] <= site_index && divergence[i2] <= site_index)){
                        searchedFollowing = true;
                        i2++;
                    }
                }

                //report set maximal matches for current sample
                if (searchedPrevious && i1 + 1 >= 0){
                    for (int j = i1 + 1; j <= i - 1; j++){
                        set_maximal_matches.push_back(tuple<int, int, int, int>(prefix[i], prefix[j], divergence[i], site_index));
                    }
                }
                if (searchedFollowing && i2 - 1 <= number_of_samples - 1){
                    for (int j = i + 1; j <= i2 - 1; j++){
                        set_maximal_matches.push_back(tuple<int, int, int, int>(prefix[i], prefix[j], divergence[i + 1], site_index));
                    }
                }
            }
        }
        return set_maximal_matches;
    }
};
