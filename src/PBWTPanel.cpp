//  * --------------------------------------------------------------------------------------------------------
//  * Name: PBWTPanel.cpp
//  * Description: Implementation of Richard Durbin's Positional Burrows-Wheeler Transform (PBWT) algorithms.
//  * Author: Yuan Wei 
//  * Created on: Jun 22, 2021
//  * --------------------------------------------------------------------------------------------------------

#include <stack>
#include <functional>
#include "PBWTHeadExtension.cpp"
#include "PBWTBucket.cpp"
using namespace std;

//Stack used for merging low resolution IBD segments on the fly with high resolution query
class PBWTStack{
    public:
    stack<tuple<int, int>> matches; //value is a tuple indicating a match [start_site_index, end_site_index]; the stack has an ascending position order matches, i.e. the top one is the most left match detected
    tuple<int, int> ur_match; //the last unreported match which overlaps with high resolution match

    PBWTStack(){
        ur_match = tuple<int, int>(-1, -1);
    }
};

class PBWTPanel{
    protected:
    int number_of_sites;
    int number_of_samples;
    vector<vector<bool>> panel_site;
    vector<vector<int>> panel_prefix;
    vector<vector<int>> panel_divergence;
    vector<vector<int>> panel_upointer;
    vector<vector<int>> panel_vpointer;

    void buildPBWTPanel(){
        for (int i = 0; i < panel_site.size(); i++){
            PBWTHeadExtension pbwt_head_extension;
            if (i == 0){
                pbwt_head_extension = PBWTHeadExtension(i, panel_site[i]);
            }
            else {
                pbwt_head_extension = PBWTHeadExtension(i, panel_site[i], panel_prefix[i - 1], panel_divergence[i - 1]);
            }
            panel_prefix.push_back(pbwt_head_extension.getPrefix());
            panel_divergence.push_back(pbwt_head_extension.getDivergence());
            panel_upointer.push_back(pbwt_head_extension.getUpointer());
            panel_vpointer.push_back(pbwt_head_extension.getVpointer());
        }
    }

    public:
    PBWTPanel(){
    }

    PBWTPanel(int number_of_sites_in, int number_of_samples_in, vector<vector<bool>> panel_site_in){
        number_of_sites = number_of_sites_in;
        number_of_samples = number_of_samples_in;
        panel_site = panel_site_in;
        buildPBWTPanel();
    }

    int getNumberOfSites(){
        return number_of_sites;
    }

    int getNumberOfSamples(){
        return number_of_samples;
    }

    vector<tuple<int, int, int>> getQuerySetMaximalMatches(const vector<bool>& query_site) const {
        //matched pair is reported as a tuple: (sample_index_matched_query, starting_position_(inclusive), ending_position_(inclusive))
        vector<tuple<int, int, int>> query_set_maximal_matches;

        //f or f_temp is the index of the first sample in the matching block between query and panel
        //g or g_temp is the index of the first sample not in the matching block between query and panel
        //e or e_temp is the start position of the sample in the matching block between query and panel (inclusive)
        int f = 0;
        int g = number_of_samples;
        int e = 0;
        int f_temp = f;
        int g_temp = g;
        int e_temp = e;
        for (int i = 0; i <= number_of_sites; i++){
            if (i < number_of_sites){
                if (f == number_of_samples){
                    if (query_site[i] == 0){
                        f_temp = panel_upointer[i][0];                    
                    }
                    else {
                        f_temp = panel_vpointer[i][0];
                    }
                }
                else {
                    if (query_site[i] == 0){
                        f_temp = panel_upointer[i][f];
                    }
                    else {
                        f_temp = panel_vpointer[i][f];
                    }
                }
                if (g == number_of_samples){
                    if (query_site[i] == 0){
                        g_temp = panel_vpointer[i][0];
                    }
                    else {
                        g_temp = number_of_samples;
                    }
                }
                else {
                    if (query_site[i] == 0){
                        g_temp = panel_upointer[i][g];
                    }
                    else {
                        g_temp = panel_vpointer[i][g];
                    }
                }
                if (f_temp < g_temp){
                    e_temp = e;
                }
                else {

                    //report set maximal matches in current matching block
                    if (!(i - 1 < 0 || i - 1 < e || e > number_of_sites - 1)){
                        for (int j = f; j < g; j++){
                            query_set_maximal_matches.push_back(tuple<int, int, int>(panel_prefix[i - 1][j], e, i - 1));
                        }
                    }

                    //find a new matching block
                    if (f_temp == 0 || f_temp == number_of_samples){
                        e_temp = i;
                    }
                    else {
                        e_temp = panel_divergence[i][f_temp] - 1;
                    }
                    if ((f_temp == number_of_samples && query_site[e_temp] == panel_site[e_temp][panel_prefix[i][number_of_samples - 1]]) || (f_temp > 0 && f_temp < number_of_samples && query_site[e_temp] == 0)){

                        //search for previous samples
                        g_temp = f_temp;
                        f_temp = f_temp - 1;
                        while (e_temp - 1 >= 0 && query_site[e_temp - 1] == panel_site[e_temp - 1][panel_prefix[i][f_temp]]){
                            e_temp--;
                        }
                        while (f_temp > 0 && panel_divergence[i][f_temp] <= e_temp){
                            f_temp--;
                        }
                    }
                    else if ((g_temp == 0 && query_site[e_temp] == panel_site[e_temp][panel_prefix[i][0]]) || (g_temp > 0 && g_temp < number_of_samples && query_site[e_temp] == 1)){

                        //search for following samples
                        f_temp = g_temp;
                        while (e_temp - 1 >= 0 && query_site[e_temp - 1] == panel_site[e_temp - 1][panel_prefix[i][g_temp]]){
                            e_temp--;
                        }
                        g_temp = g_temp + 1;
                        while (g_temp < number_of_samples && panel_divergence[i][g_temp] <= e_temp){
                            g_temp++;
                        }
                    }
                    else {
                        f_temp = 0;
                        g_temp = number_of_samples;
                        e_temp = i + 1;
                    }
                }
                f = f_temp;
                g = g_temp;
                e = e_temp;
            }
            else {
                
                //report set maximal matches in current matching block if all sites are processed
                if (!(i - 1 < 0 || i - 1 < e || e > number_of_sites - 1)){
                    for (int j = f; j < g; j++){
                        query_set_maximal_matches.push_back(tuple<int, int, int>(panel_prefix[i - 1][j], e, i - 1));
                    }
                }
            }
        }
        return query_set_maximal_matches;
    }

    vector<tuple<int, int, int>> getQueryLongMatches(const vector<bool>& query_site, int l) const {
        //matched pair is reported as a tuple: (sample_index_matched_query, starting_position_(inclusive), ending_position_(inclusive))
        vector<tuple<int, int, int>> query_long_matches;

        //track the closest cutoff length in sites; i_l is a site index
        int i_l = 0;

        //f_e or f_e_temp is the index of the first sample in the searching block between query and panel
        //g_e or g_e_temp is the index of the first sample not in the searching block between query and panel
        //e or e_temp is the starting position of the sample in the searching block between query and panel (inclusive)
        //f_l or f_l_temp is the index of the first sample in the matching block between query and panel
        //g_l or g_l_temp is the index of the first sample not in the matching block between query and panel
        //dz stores the starting positions of the samples in the matching block between query and panel (inclusive)
        int f_e = 0;
        int g_e = number_of_samples;
        int e = 0;
        int f_e_temp = f_e;
        int g_e_temp = g_e;
        int e_temp = e;
        int f_l = 0;
        int g_l = 0;
        int f_l_temp = f_l;
        int g_l_temp = g_l;
        unordered_map<int, int> dz;
        for (int i = 0; i <= number_of_sites; i++){
            if (i < number_of_sites){
                //track the closest cutoff length in sites
                if (i_l < 0){
                    i_l = 0;
                }
                while (i_l >= 0 && i_l < i && i + 1 - i_l > l){
                    i_l++;
                }

                //report long matches in the matching block if they break in current iteration
                if (!dz.empty()){
                    int f_r = 0;
                    int g_r = 0;
                    if (f_l == number_of_samples){
                        if (!(query_site[i] == 0)){
                            f_r = panel_upointer[i][0];                    
                        }
                        else {
                            f_r = panel_vpointer[i][0];
                        }
                    }
                    else {
                        if (!(query_site[i] == 0)){
                            f_r = panel_upointer[i][f_l];
                        }
                        else {
                            f_r = panel_vpointer[i][f_l];
                        }
                    }
                    if (g_l == number_of_samples){
                        if (!(query_site[i] == 0)){
                            g_r = panel_vpointer[i][0];
                        }
                        else {
                            g_r = number_of_samples;
                        }
                    }
                    else {
                        if (!(query_site[i] == 0)){
                            g_r = panel_upointer[i][g_l];
                        }
                        else {
                            g_r = panel_vpointer[i][g_l];
                        }
                    }
                    for (int j = f_r; j < g_r; j++){
                        if (dz.count(panel_prefix[i][j]) > 0){
                            int start_position = dz[panel_prefix[i][j]];
                            int end_position = i - 1;
                            if (!(end_position < 0 || end_position < start_position || start_position > number_of_sites - 1)){
                                query_long_matches.push_back(tuple<int, int, int>(panel_prefix[i][j], start_position, end_position));
                            }
                            dz.erase(panel_prefix[i][j]);
                        }
                    }
                }

                //update searching block
                if (f_e == number_of_samples){
                    if (query_site[i] == 0){
                        f_e_temp = panel_upointer[i][0];                    
                    }
                    else {
                        f_e_temp = panel_vpointer[i][0];
                    }
                }
                else {
                    if (query_site[i] == 0){
                        f_e_temp = panel_upointer[i][f_e];
                    }
                    else {
                        f_e_temp = panel_vpointer[i][f_e];
                    }
                }
                if (g_e == number_of_samples){
                    if (query_site[i] == 0){
                        g_e_temp = panel_vpointer[i][0];
                    }
                    else {
                        g_e_temp = number_of_samples;
                    }
                }
                else {
                    if (query_site[i] == 0){
                        g_e_temp = panel_upointer[i][g_e];
                    }
                    else {
                        g_e_temp = panel_vpointer[i][g_e];
                    }
                }
                if (f_e_temp < g_e_temp){
                    e_temp = e;
                }
                else {

                    //find a new searching block
                    if (f_e_temp == 0 || f_e_temp == number_of_samples){
                        e_temp = i;
                    }
                    else {
                        e_temp = panel_divergence[i][f_e_temp] - 1;
                    }
                    if ((f_e_temp == number_of_samples && query_site[e_temp] == panel_site[e_temp][panel_prefix[i][number_of_samples - 1]]) || (f_e_temp > 0 && f_e_temp < number_of_samples && query_site[e_temp] == 0)){

                        //search for previous samples
                        g_e_temp = f_e_temp;
                        f_e_temp = f_e_temp - 1;
                        while (e_temp - 1 >= 0 && query_site[e_temp - 1] == panel_site[e_temp - 1][panel_prefix[i][f_e_temp]]){
                            e_temp--;
                        }
                        while (f_e_temp > 0 && panel_divergence[i][f_e_temp] <= e_temp){
                            f_e_temp--;
                        }
                    }
                    else if ((g_e_temp == 0 && query_site[e_temp] == panel_site[e_temp][panel_prefix[i][0]]) || (g_e_temp > 0 && g_e_temp < number_of_samples && query_site[e_temp] == 1)){

                        //search for following samples
                        f_e_temp = g_e_temp;
                        while (e_temp - 1 >= 0 && query_site[e_temp - 1] == panel_site[e_temp - 1][panel_prefix[i][g_e_temp]]){
                            e_temp--;
                        }
                        g_e_temp = g_e_temp + 1;
                        while (g_e_temp < number_of_samples && panel_divergence[i][g_e_temp] <= e_temp){
                            g_e_temp++;
                        }
                    }
                    else {
                        f_e_temp = 0;
                        g_e_temp = number_of_samples;
                        e_temp = i + 1;
                    }
                }
                f_e = f_e_temp;
                g_e = g_e_temp;
                e = e_temp;

                //update matching block
                if (f_l < g_l){
                    if (f_l == number_of_samples){
                        if (query_site[i] == 0){
                            f_l_temp = panel_upointer[i][0];                    
                        }
                        else {
                            f_l_temp = panel_vpointer[i][0];
                        }
                    }
                    else {
                        if (query_site[i] == 0){
                            f_l_temp = panel_upointer[i][f_l];
                        }
                        else {
                            f_l_temp = panel_vpointer[i][f_l];
                        }
                    }
                    if (g_l == number_of_samples){
                        if (query_site[i] == 0){
                            g_l_temp = panel_vpointer[i][0];
                        }
                        else {
                            g_l_temp = number_of_samples;
                        }
                    }
                    else {
                        if (query_site[i] == 0){
                            g_l_temp = panel_upointer[i][g_l];
                        }
                        else {
                            g_l_temp = panel_vpointer[i][g_l];
                        }
                    }
                }
                else {
                    f_l_temp = f_l;
                    g_l_temp = g_l;
                }

                //find a new matching block
                if (f_l_temp == g_l_temp){
                    if (e <= i + 1 - l){
                        for (int j = f_e; j < g_e; j++){
                            dz[panel_prefix[i][j]] = e;
                        }
                        f_l_temp = f_e;
                        g_l_temp = g_e;
                    }
                }

                //expand existing matching block
                if (f_l_temp < g_l_temp){
                    while (f_l_temp > 0 && panel_divergence[i][f_l_temp] <= i + 1 - l){
                        f_l_temp--;
                        dz[panel_prefix[i][f_l_temp]] = i_l;
                    }
                    while (g_l_temp < number_of_samples && panel_divergence[i][g_l_temp] <= i + 1 - l){
                        dz[panel_prefix[i][g_l_temp]] = i_l;
                        g_l_temp++;
                    }
                }
                f_l = f_l_temp;
                g_l = g_l_temp;
            }
            else {

                //report long matches in the matching block if all sites are processed
                for (int j = f_l; j < g_l; j++){
                    if (dz.count(panel_prefix[i - 1][j]) > 0){
                        int start_position = dz[panel_prefix[i - 1][j]];
                        int end_position = i - 1;
                        if (!(end_position < 0 || end_position < start_position || start_position > number_of_sites - 1)){
                            query_long_matches.push_back(tuple<int, int, int>(panel_prefix[i - 1][j], start_position, end_position));
                        }
                        dz.erase(panel_prefix[i - 1][j]);
                    }
                }
            }
        }
        dz.clear();
        return query_long_matches;
    }

    vector<tuple<int, int, int>> getQueryLongMatchesWithGeneticDistance(const vector<bool>& query_site, double genetic_distance, int number_of_markers, const vector<double>& genetic_positions, int panel_window_size) const {
        //matched pair is reported as a tuple: (sample_index_matched_query, starting_position_(inclusive), ending_position_(inclusive))
        //genetic distance is the first requirement: minimal genetic distance a match should have; number of markers is the second requirement: minimal number of site index distance a match should have
        vector<tuple<int, int, int>> query_long_matches;
        int total_number_of_sites_in_panel = genetic_positions.size();

        //track the closest cutoff length in sites, equivalent to the cutoff genetic length in centiMorgans; i_l is a projected site index
        int i_l = 0;
        
        //f_e or f_e_temp is the index of the first sample in the searching block between query and panel
        //g_e or g_e_temp is the index of the first sample not in the searching block between query and panel
        //e or e_temp is the starting position of the sample in the searching block between query and panel (inclusive)
        //f_l or f_l_temp is the index of the first sample in the matching block between query and panel
        //g_l or g_l_temp is the index of the first sample not in the matching block between query and panel
        //dz stores the starting positions of the samples in the matching block between query and panel (inclusive)
        int f_e = 0;
        int g_e = number_of_samples;
        int e = 0;
        int f_e_temp = f_e;
        int g_e_temp = g_e;
        int e_temp = e;
        int f_l = 0;
        int g_l = 0;
        int f_l_temp = f_l;
        int g_l_temp = g_l;
        unordered_map<int, int> dz;
        for (int i = 0; i <= number_of_sites; i++){
            if (i < number_of_sites){
                //track the closest cutoff length in sites
                if (i_l < 0){
                    i_l = 0;
                }
                while (i_l >= 0 && i_l < i && isProjectedSiteIndicesAtLeastGeneticDistance(i_l, i, total_number_of_sites_in_panel, panel_window_size, genetic_distance, genetic_positions)){
                    i_l++;
                }
                i_l--;

                //report long matches in the matching block if they break in current iteration
                if (!dz.empty()){
                    int f_r = 0;
                    int g_r = 0;
                    if (f_l == number_of_samples){
                        if (!(query_site[i] == 0)){
                            f_r = panel_upointer[i][0];                    
                        }
                        else {
                            f_r = panel_vpointer[i][0];
                        }
                    }
                    else {
                        if (!(query_site[i] == 0)){
                            f_r = panel_upointer[i][f_l];
                        }
                        else {
                            f_r = panel_vpointer[i][f_l];
                        }
                    }
                    if (g_l == number_of_samples){
                        if (!(query_site[i] == 0)){
                            g_r = panel_vpointer[i][0];
                        }
                        else {
                            g_r = number_of_samples;
                        }
                    }
                    else {
                        if (!(query_site[i] == 0)){
                            g_r = panel_upointer[i][g_l];
                        }
                        else {
                            g_r = panel_vpointer[i][g_l];
                        }
                    }
                    for (int j = f_r; j < g_r; j++){
                        if (dz.count(panel_prefix[i][j]) > 0){
                            int start_position = dz[panel_prefix[i][j]];
                            int end_position = i - 1;
                            if (!(end_position < 0 || end_position < start_position || start_position > number_of_sites - 1)){
                                query_long_matches.push_back(tuple<int, int, int>(panel_prefix[i][j], start_position, end_position));
                            }
                            dz.erase(panel_prefix[i][j]);
                        }
                    }
                }

                //update searching block
                if (f_e == number_of_samples){
                    if (query_site[i] == 0){
                        f_e_temp = panel_upointer[i][0];                    
                    }
                    else {
                        f_e_temp = panel_vpointer[i][0];
                    }
                }
                else {
                    if (query_site[i] == 0){
                        f_e_temp = panel_upointer[i][f_e];
                    }
                    else {
                        f_e_temp = panel_vpointer[i][f_e];
                    }
                }
                if (g_e == number_of_samples){
                    if (query_site[i] == 0){
                        g_e_temp = panel_vpointer[i][0];
                    }
                    else {
                        g_e_temp = number_of_samples;
                    }
                }
                else {
                    if (query_site[i] == 0){
                        g_e_temp = panel_upointer[i][g_e];
                    }
                    else {
                        g_e_temp = panel_vpointer[i][g_e];
                    }
                }
                if (f_e_temp < g_e_temp){
                    e_temp = e;
                }
                else {

                    //find a new searching block
                    if (f_e_temp == 0 || f_e_temp == number_of_samples){
                        e_temp = i;
                    }
                    else {
                        e_temp = panel_divergence[i][f_e_temp] - 1;
                    }
                    if ((f_e_temp == number_of_samples && query_site[e_temp] == panel_site[e_temp][panel_prefix[i][number_of_samples - 1]]) || (f_e_temp > 0 && f_e_temp < number_of_samples && query_site[e_temp] == 0)){

                        //search for previous samples
                        g_e_temp = f_e_temp;
                        f_e_temp = f_e_temp - 1;
                        while (e_temp - 1 >= 0 && query_site[e_temp - 1] == panel_site[e_temp - 1][panel_prefix[i][f_e_temp]]){
                            e_temp--;
                        }
                        while (f_e_temp > 0 && panel_divergence[i][f_e_temp] <= e_temp){
                            f_e_temp--;
                        }
                    }
                    else if ((g_e_temp == 0 && query_site[e_temp] == panel_site[e_temp][panel_prefix[i][0]]) || (g_e_temp > 0 && g_e_temp < number_of_samples && query_site[e_temp] == 1)){

                        //search for following samples
                        f_e_temp = g_e_temp;
                        while (e_temp - 1 >= 0 && query_site[e_temp - 1] == panel_site[e_temp - 1][panel_prefix[i][g_e_temp]]){
                            e_temp--;
                        }
                        g_e_temp = g_e_temp + 1;
                        while (g_e_temp < number_of_samples && panel_divergence[i][g_e_temp] <= e_temp){
                            g_e_temp++;
                        }
                    }
                    else {
                        f_e_temp = 0;
                        g_e_temp = number_of_samples;
                        e_temp = i + 1;
                    }
                }
                f_e = f_e_temp;
                g_e = g_e_temp;
                e = e_temp;

                //update matching block
                if (f_l < g_l){
                    if (f_l == number_of_samples){
                        if (query_site[i] == 0){
                            f_l_temp = panel_upointer[i][0];                    
                        }
                        else {
                            f_l_temp = panel_vpointer[i][0];
                        }
                    }
                    else {
                        if (query_site[i] == 0){
                            f_l_temp = panel_upointer[i][f_l];
                        }
                        else {
                            f_l_temp = panel_vpointer[i][f_l];
                        }
                    }
                    if (g_l == number_of_samples){
                        if (query_site[i] == 0){
                            g_l_temp = panel_vpointer[i][0];
                        }
                        else {
                            g_l_temp = number_of_samples;
                        }
                    }
                    else {
                        if (query_site[i] == 0){
                            g_l_temp = panel_upointer[i][g_l];
                        }
                        else {
                            g_l_temp = panel_vpointer[i][g_l];
                        }
                    }
                }
                else {
                    f_l_temp = f_l;
                    g_l_temp = g_l;
                }

                //find a new matching block
                if (f_l_temp == g_l_temp){
                    if (isProjectedSiteIndicesAtLeastGeneticDistance(e, i, total_number_of_sites_in_panel, panel_window_size, genetic_distance, genetic_positions) && e <= i + 1 - number_of_markers){
                        //add sequences from searching block to form a new matching block
                        for (int j = f_e; j < g_e; j++){
                            dz[panel_prefix[i][j]] = e;
                        }
                        f_l_temp = f_e;
                        g_l_temp = g_e;
                    }
                }

                //expand existing matching block
                if (f_l_temp < g_l_temp){
                    while (f_l_temp > 0 && isProjectedSiteIndicesAtLeastGeneticDistance(panel_divergence[i][f_l_temp], i, total_number_of_sites_in_panel, panel_window_size, genetic_distance, genetic_positions) && panel_divergence[i][f_l_temp] <= i + 1 - number_of_markers){
                        f_l_temp--;
                        dz[panel_prefix[i][f_l_temp]] = i_l;
                    }
                    while (g_l_temp < number_of_samples && isProjectedSiteIndicesAtLeastGeneticDistance(panel_divergence[i][g_l_temp], i, total_number_of_sites_in_panel, panel_window_size, genetic_distance, genetic_positions) && panel_divergence[i][g_l_temp] <= i + 1 - number_of_markers){
                        dz[panel_prefix[i][g_l_temp]] = i_l;
                        g_l_temp++;
                    }
                }
                f_l = f_l_temp;
                g_l = g_l_temp;
            }
            else {

                //report long matches in the matching block if all sites are processed
                for (int j = f_l; j < g_l; j++){
                    if (dz.count(panel_prefix[i - 1][j]) > 0){
                        int start_position = dz[panel_prefix[i - 1][j]];
                        int end_position = i - 1;
                        if (!(end_position < 0 || end_position < start_position || start_position > number_of_sites - 1)){
                            query_long_matches.push_back(tuple<int, int, int>(panel_prefix[i - 1][j], start_position, end_position));
                        }
                        dz.erase(panel_prefix[i - 1][j]);
                    }
                }
            }
        }
        dz.clear();
        return query_long_matches;
    }

    void getQueryLongMatchesWithGeneticDistance(const vector<bool>& query_site, double genetic_distance, int number_of_markers, const vector<double>& genetic_positions, int panel_window_size, PBWTBucket& bucket) const {
        //matched pair is reported as a tuple: (sample_index_matched_query, starting_position_(inclusive), ending_position_(inclusive))
        //genetic distance is the first requirement: minimal genetic distance a match should have; number of markers is the second requirement: minimal number of site index distance a match should have
        int total_number_of_sites_in_panel = genetic_positions.size();

        //track the closest cutoff length in sites, equivalent to the cutoff genetic length in centiMorgans; i_l is a projected site index
        int i_l = 0;

        //f_e or f_e_temp is the index of the first sample in the searching block between query and panel
        //g_e or g_e_temp is the index of the first sample not in the searching block between query and panel
        //e or e_temp is the starting position of the sample in the searching block between query and panel (inclusive)
        //f_l or f_l_temp is the index of the first sample in the matching block between query and panel
        //g_l or g_l_temp is the index of the first sample not in the matching block between query and panel
        //dz stores the starting positions of the samples in the matching block between query and panel (inclusive)
        int f_e = 0;
        int g_e = number_of_samples;
        int e = 0;
        int f_e_temp = f_e;
        int g_e_temp = g_e;
        int e_temp = e;
        int f_l = 0;
        int g_l = 0;
        int f_l_temp = f_l;
        int g_l_temp = g_l;
        unordered_map<int, int> dz;
        for (int i = 0; i <= number_of_sites; i++){
            if (i < number_of_sites){
                //track the closest cutoff length in sites
                if (i_l < 0){
                    i_l = 0;
                }
                while (i_l >= 0 && i_l < i && isProjectedSiteIndicesAtLeastGeneticDistance(i_l, i, total_number_of_sites_in_panel, panel_window_size, genetic_distance, genetic_positions)){
                    i_l++;
                }
                i_l--;

                //report long matches in the matching block if they break in current iteration
                if (!dz.empty()){
                    int f_r = 0;
                    int g_r = 0;
                    if (f_l == number_of_samples){
                        if (!(query_site[i] == 0)){
                            f_r = panel_upointer[i][0];                    
                        }
                        else {
                            f_r = panel_vpointer[i][0];
                        }
                    }
                    else {
                        if (!(query_site[i] == 0)){
                            f_r = panel_upointer[i][f_l];
                        }
                        else {
                            f_r = panel_vpointer[i][f_l];
                        }
                    }
                    if (g_l == number_of_samples){
                        if (!(query_site[i] == 0)){
                            g_r = panel_vpointer[i][0];
                        }
                        else {
                            g_r = number_of_samples;
                        }
                    }
                    else {
                        if (!(query_site[i] == 0)){
                            g_r = panel_upointer[i][g_l];
                        }
                        else {
                            g_r = panel_vpointer[i][g_l];
                        }
                    }
                    for (int j = f_r; j < g_r; j++){
                        if (dz.count(panel_prefix[i][j]) > 0){
                            int start_position = dz[panel_prefix[i][j]];
                            int end_position = i - 1;
                            if (!(end_position < 0 || end_position < start_position || start_position > number_of_sites - 1)){
                                if (panel_window_size > 1){
                                    //rescale match
                                    if (start_position > end_position){
                                        throw invalid_argument("invalid range");
                                    }
                                    start_position = start_position * panel_window_size;
                                    end_position = (end_position + 1) * panel_window_size - 1 < total_number_of_sites_in_panel ? (end_position + 1) * panel_window_size - 1 : total_number_of_sites_in_panel - 1;
                                }
                                bucket.addMatch(panel_prefix[i][j], start_position, end_position, panel_window_size, false);
                            }
                            dz.erase(panel_prefix[i][j]);
                        }
                    }
                }

                //update searching block
                if (f_e == number_of_samples){
                    if (query_site[i] == 0){
                        f_e_temp = panel_upointer[i][0];                    
                    }
                    else {
                        f_e_temp = panel_vpointer[i][0];
                    }
                }
                else {
                    if (query_site[i] == 0){
                        f_e_temp = panel_upointer[i][f_e];
                    }
                    else {
                        f_e_temp = panel_vpointer[i][f_e];
                    }
                }
                if (g_e == number_of_samples){
                    if (query_site[i] == 0){
                        g_e_temp = panel_vpointer[i][0];
                    }
                    else {
                        g_e_temp = number_of_samples;
                    }
                }
                else {
                    if (query_site[i] == 0){
                        g_e_temp = panel_upointer[i][g_e];
                    }
                    else {
                        g_e_temp = panel_vpointer[i][g_e];
                    }
                }
                if (f_e_temp < g_e_temp){
                    e_temp = e;
                }
                else {

                    //find a new searching block
                    if (f_e_temp == 0 || f_e_temp == number_of_samples){
                        e_temp = i;
                    }
                    else {
                        e_temp = panel_divergence[i][f_e_temp] - 1;
                    }
                    if ((f_e_temp == number_of_samples && query_site[e_temp] == panel_site[e_temp][panel_prefix[i][number_of_samples - 1]]) || (f_e_temp > 0 && f_e_temp < number_of_samples && query_site[e_temp] == 0)){

                        //search for previous samples
                        g_e_temp = f_e_temp;
                        f_e_temp = f_e_temp - 1;
                        while (e_temp - 1 >= 0 && query_site[e_temp - 1] == panel_site[e_temp - 1][panel_prefix[i][f_e_temp]]){
                            e_temp--;
                        }
                        while (f_e_temp > 0 && panel_divergence[i][f_e_temp] <= e_temp){
                            f_e_temp--;
                        }
                    }
                    else if ((g_e_temp == 0 && query_site[e_temp] == panel_site[e_temp][panel_prefix[i][0]]) || (g_e_temp > 0 && g_e_temp < number_of_samples && query_site[e_temp] == 1)){

                        //search for following samples
                        f_e_temp = g_e_temp;
                        while (e_temp - 1 >= 0 && query_site[e_temp - 1] == panel_site[e_temp - 1][panel_prefix[i][g_e_temp]]){
                            e_temp--;
                        }
                        g_e_temp = g_e_temp + 1;
                        while (g_e_temp < number_of_samples && panel_divergence[i][g_e_temp] <= e_temp){
                            g_e_temp++;
                        }
                    }
                    else {
                        f_e_temp = 0;
                        g_e_temp = number_of_samples;
                        e_temp = i + 1;
                    }
                }
                f_e = f_e_temp;
                g_e = g_e_temp;
                e = e_temp;

                //update matching block
                if (f_l < g_l){
                    if (f_l == number_of_samples){
                        if (query_site[i] == 0){
                            f_l_temp = panel_upointer[i][0];                    
                        }
                        else {
                            f_l_temp = panel_vpointer[i][0];
                        }
                    }
                    else {
                        if (query_site[i] == 0){
                            f_l_temp = panel_upointer[i][f_l];
                        }
                        else {
                            f_l_temp = panel_vpointer[i][f_l];
                        }
                    }
                    if (g_l == number_of_samples){
                        if (query_site[i] == 0){
                            g_l_temp = panel_vpointer[i][0];
                        }
                        else {
                            g_l_temp = number_of_samples;
                        }
                    }
                    else {
                        if (query_site[i] == 0){
                            g_l_temp = panel_upointer[i][g_l];
                        }
                        else {
                            g_l_temp = panel_vpointer[i][g_l];
                        }
                    }
                }
                else {
                    f_l_temp = f_l;
                    g_l_temp = g_l;
                }

                //find a new matching block
                if (f_l_temp == g_l_temp){
                    if (isProjectedSiteIndicesAtLeastGeneticDistance(e, i, total_number_of_sites_in_panel, panel_window_size, genetic_distance, genetic_positions) && e <= i + 1 - number_of_markers){
                        //add sequences from searching block to form a new matching block
                        for (int j = f_e; j < g_e; j++){
                            dz[panel_prefix[i][j]] = e;
                        }
                        f_l_temp = f_e;
                        g_l_temp = g_e;
                    }
                }

                //expand existing matching block
                if (f_l_temp < g_l_temp){
                    while (f_l_temp > 0 && isProjectedSiteIndicesAtLeastGeneticDistance(panel_divergence[i][f_l_temp], i, total_number_of_sites_in_panel, panel_window_size, genetic_distance, genetic_positions) && panel_divergence[i][f_l_temp] <= i + 1 - number_of_markers){
                        f_l_temp--;
                        dz[panel_prefix[i][f_l_temp]] = i_l;
                    }
                    while (g_l_temp < number_of_samples && isProjectedSiteIndicesAtLeastGeneticDistance(panel_divergence[i][g_l_temp], i, total_number_of_sites_in_panel, panel_window_size, genetic_distance, genetic_positions) && panel_divergence[i][g_l_temp] <= i + 1 - number_of_markers){
                        dz[panel_prefix[i][g_l_temp]] = i_l;
                        g_l_temp++;
                    }
                }
                f_l = f_l_temp;
                g_l = g_l_temp;
            }
            else {

                //report long matches in the matching block if all sites are processed
                for (int j = f_l; j < g_l; j++){
                    if (dz.count(panel_prefix[i - 1][j]) > 0){
                        int start_position = dz[panel_prefix[i - 1][j]];
                        int end_position = i - 1;
                        if (!(end_position < 0 || end_position < start_position || start_position > number_of_sites - 1)){
                            if (panel_window_size > 1){
                                //rescale match
                                if (start_position > end_position){
                                    throw invalid_argument("invalid range");
                                }
                                start_position = start_position * panel_window_size;
                                end_position = (end_position + 1) * panel_window_size - 1 < total_number_of_sites_in_panel ? (end_position + 1) * panel_window_size - 1 : total_number_of_sites_in_panel - 1;
                            }
                            bucket.addMatch(panel_prefix[i - 1][j], start_position, end_position, panel_window_size, false);
                        }
                        dz.erase(panel_prefix[i - 1][j]);
                    }
                }
            }
        }
        dz.clear();
    }

    void getQueryLongMatchesWithGeneticDistanceWhileMerge(string sample_id_query, string haplotype_id_in_sample_id_query, const vector<bool>& query_site, double minimum_length_centimorgans, int number_of_markers, double maximum_gap_length_genetic, double minimum_length_centimorgans_low_resolution, const vector<string>& sample_ids, const vector<int>& physical_positions, const vector<double>& genetic_positions, int panel_window_size, string output_file_path_and_name, unordered_map<int, PBWTStack>& match_stacks, int& total_matches_of_the_query, double& total_match_length_genetic, function<void(int, int, int, unordered_map<int, PBWTStack>&, double, double, string, string, string, const vector<string>&, const vector<int>&, const vector<double>&, int&, double&)> merge_on_the_fly) const {
        //matched pair is reported as a tuple: (sample_index_matched_query, starting_position_(inclusive), ending_position_(inclusive))
        int total_number_of_sites_in_panel = genetic_positions.size();

        //track the closest cutoff length in sites, equivalent to the cutoff genetic length in centiMorgans; i_l is a projected site index
        int i_l = 0;

        //f_e or f_e_temp is the index of the first sample in the searching block between query and panel
        //g_e or g_e_temp is the index of the first sample not in the searching block between query and panel
        //e or e_temp is the starting position of the sample in the searching block between query and panel (inclusive)
        //f_l or f_l_temp is the index of the first sample in the matching block between query and panel
        //g_l or g_l_temp is the index of the first sample not in the matching block between query and panel
        //dz stores the starting positions of the samples in the matching block between query and panel (inclusive)
        int f_e = 0;
        int g_e = number_of_samples;
        int e = 0;
        int f_e_temp = f_e;
        int g_e_temp = g_e;
        int e_temp = e;
        int f_l = 0;
        int g_l = 0;
        int f_l_temp = f_l;
        int g_l_temp = g_l;
        unordered_map<int, int> dz;
        for (int i = 0; i <= number_of_sites; i++){
            if (i < number_of_sites){
                if (i_l < 0){
                    i_l = 0;
                }
                while (i_l >= 0 && i_l < i && isProjectedSiteIndicesAtLeastGeneticDistance(i_l, i, total_number_of_sites_in_panel, panel_window_size, minimum_length_centimorgans, genetic_positions)){
                    i_l++;
                }
                i_l--;

                //report long matches in the matching block if they break in current iteration
                if (!dz.empty()){
                    int f_r = 0;
                    int g_r = 0;
                    if (f_l == number_of_samples){
                        if (!(query_site[i] == 0)){
                            f_r = panel_upointer[i][0];                    
                        }
                        else {
                            f_r = panel_vpointer[i][0];
                        }
                    }
                    else {
                        if (!(query_site[i] == 0)){
                            f_r = panel_upointer[i][f_l];
                        }
                        else {
                            f_r = panel_vpointer[i][f_l];
                        }
                    }
                    if (g_l == number_of_samples){
                        if (!(query_site[i] == 0)){
                            g_r = panel_vpointer[i][0];
                        }
                        else {
                            g_r = number_of_samples;
                        }
                    }
                    else {
                        if (!(query_site[i] == 0)){
                            g_r = panel_upointer[i][g_l];
                        }
                        else {
                            g_r = panel_vpointer[i][g_l];
                        }
                    }
                    for (int j = f_r; j < g_r; j++){
                        if (dz.count(panel_prefix[i][j]) > 0){
                            int start_position = dz[panel_prefix[i][j]];
                            int end_position = i - 1;
                            if (!(end_position < 0 || end_position < start_position || start_position > number_of_sites - 1)){

                                if (panel_window_size > 1){
                                    //rescale match
                                    if (start_position > end_position){
                                        throw invalid_argument("invalid range");
                                    }
                                    start_position = start_position * panel_window_size;
                                    end_position = (end_position + 1) * panel_window_size - 1 < total_number_of_sites_in_panel ? (end_position + 1) * panel_window_size - 1 : total_number_of_sites_in_panel - 1;
                                }

                                //use current detected high resolution match (panel_prefix[i][j], start_position, end_position) to merge low resolution matches in PBWT stack
                                //merge_on_the_fly function parameters: int haplotype_id_index_other, int start_site_index, int end_site_index, unordered_map<int, PBWTStack>& match_stacks, double minimum_length_centimorgans_low_resolution, double maximum_gap_length_genetic, string output_file_path_and_name, string sample_id_query, string haplotype_id_in_sample_id_query, const vector<string>& sample_ids, const vector<int>& physical_positions, const vector<double>& genetic_positions, int& total_matches_of_the_query, double& total_match_length_genetic
                                merge_on_the_fly(panel_prefix[i][j], start_position, end_position, match_stacks, minimum_length_centimorgans_low_resolution, maximum_gap_length_genetic, output_file_path_and_name, sample_id_query, haplotype_id_in_sample_id_query, sample_ids, physical_positions, genetic_positions, total_matches_of_the_query, total_match_length_genetic);
                            }
                            dz.erase(panel_prefix[i][j]);
                        }
                    }
                }

                //update searching block
                if (f_e == number_of_samples){
                    if (query_site[i] == 0){
                        f_e_temp = panel_upointer[i][0];                    
                    }
                    else {
                        f_e_temp = panel_vpointer[i][0];
                    }
                }
                else {
                    if (query_site[i] == 0){
                        f_e_temp = panel_upointer[i][f_e];
                    }
                    else {
                        f_e_temp = panel_vpointer[i][f_e];
                    }
                }
                if (g_e == number_of_samples){
                    if (query_site[i] == 0){
                        g_e_temp = panel_vpointer[i][0];
                    }
                    else {
                        g_e_temp = number_of_samples;
                    }
                }
                else {
                    if (query_site[i] == 0){
                        g_e_temp = panel_upointer[i][g_e];
                    }
                    else {
                        g_e_temp = panel_vpointer[i][g_e];
                    }
                }
                if (f_e_temp < g_e_temp){
                    e_temp = e;
                }
                else {

                    //find a new searching block
                    if (f_e_temp == 0 || f_e_temp == number_of_samples){
                        e_temp = i;
                    }
                    else {
                        e_temp = panel_divergence[i][f_e_temp] - 1;
                    }
                    if ((f_e_temp == number_of_samples && query_site[e_temp] == panel_site[e_temp][panel_prefix[i][number_of_samples - 1]]) || (f_e_temp > 0 && f_e_temp < number_of_samples && query_site[e_temp] == 0)){

                        //search for previous samples
                        g_e_temp = f_e_temp;
                        f_e_temp = f_e_temp - 1;
                        while (e_temp - 1 >= 0 && query_site[e_temp - 1] == panel_site[e_temp - 1][panel_prefix[i][f_e_temp]]){
                            e_temp--;
                        }
                        while (f_e_temp > 0 && panel_divergence[i][f_e_temp] <= e_temp){
                            f_e_temp--;
                        }
                    }
                    else if ((g_e_temp == 0 && query_site[e_temp] == panel_site[e_temp][panel_prefix[i][0]]) || (g_e_temp > 0 && g_e_temp < number_of_samples && query_site[e_temp] == 1)){

                        //search for following samples
                        f_e_temp = g_e_temp;
                        while (e_temp - 1 >= 0 && query_site[e_temp - 1] == panel_site[e_temp - 1][panel_prefix[i][g_e_temp]]){
                            e_temp--;
                        }
                        g_e_temp = g_e_temp + 1;
                        while (g_e_temp < number_of_samples && panel_divergence[i][g_e_temp] <= e_temp){
                            g_e_temp++;
                        }
                    }
                    else {
                        f_e_temp = 0;
                        g_e_temp = number_of_samples;
                        e_temp = i + 1;
                    }
                }
                f_e = f_e_temp;
                g_e = g_e_temp;
                e = e_temp;

                //update matching block
                if (f_l < g_l){
                    if (f_l == number_of_samples){
                        if (query_site[i] == 0){
                            f_l_temp = panel_upointer[i][0];                    
                        }
                        else {
                            f_l_temp = panel_vpointer[i][0];
                        }
                    }
                    else {
                        if (query_site[i] == 0){
                            f_l_temp = panel_upointer[i][f_l];
                        }
                        else {
                            f_l_temp = panel_vpointer[i][f_l];
                        }
                    }
                    if (g_l == number_of_samples){
                        if (query_site[i] == 0){
                            g_l_temp = panel_vpointer[i][0];
                        }
                        else {
                            g_l_temp = number_of_samples;
                        }
                    }
                    else {
                        if (query_site[i] == 0){
                            g_l_temp = panel_upointer[i][g_l];
                        }
                        else {
                            g_l_temp = panel_vpointer[i][g_l];
                        }
                    }
                }
                else {
                    f_l_temp = f_l;
                    g_l_temp = g_l;
                }

                //find a new matching block
                if (f_l_temp == g_l_temp){
                    if (isProjectedSiteIndicesAtLeastGeneticDistance(e, i, total_number_of_sites_in_panel, panel_window_size, minimum_length_centimorgans, genetic_positions) && e <= i + 1 - number_of_markers){
                        //add sequences from searching block to form a new matching block
                        for (int j = f_e; j < g_e; j++){
                            dz[panel_prefix[i][j]] = e;
                        }
                        f_l_temp = f_e;
                        g_l_temp = g_e;
                    }
                }

                //expand existing matching block
                if (f_l_temp < g_l_temp){
                    while (f_l_temp > 0 && isProjectedSiteIndicesAtLeastGeneticDistance(panel_divergence[i][f_l_temp], i, total_number_of_sites_in_panel, panel_window_size, minimum_length_centimorgans, genetic_positions) && panel_divergence[i][f_l_temp] <= i + 1 - number_of_markers){
                        f_l_temp--;
                        dz[panel_prefix[i][f_l_temp]] = i_l;
                    }
                    while (g_l_temp < number_of_samples && isProjectedSiteIndicesAtLeastGeneticDistance(panel_divergence[i][g_l_temp], i, total_number_of_sites_in_panel, panel_window_size, minimum_length_centimorgans, genetic_positions) && panel_divergence[i][g_l_temp] <= i + 1 - number_of_markers){
                        dz[panel_prefix[i][g_l_temp]] = i_l;
                        g_l_temp++;
                    }
                }
                f_l = f_l_temp;
                g_l = g_l_temp;
            }
            else {

                //report long matches in the matching block if all sites are processed
                for (int j = f_l; j < g_l; j++){
                    if (dz.count(panel_prefix[i - 1][j]) > 0){
                        int start_position = dz[panel_prefix[i - 1][j]];
                        int end_position = i - 1;
                        if (!(end_position < 0 || end_position < start_position || start_position > number_of_sites - 1)){

                            if (panel_window_size > 1){
                                //rescale match
                                if (start_position > end_position){
                                    throw invalid_argument("invalid range");
                                }
                                start_position = start_position * panel_window_size;
                                end_position = (end_position + 1) * panel_window_size - 1 < total_number_of_sites_in_panel ? (end_position + 1) * panel_window_size - 1 : total_number_of_sites_in_panel - 1;
                            }

                            //use current detected high resolution match (panel_prefix[i - 1][j], start_position, end_position) to merge low resolution matches in PBWT stack
                            //merge_on_the_fly function parameters: int haplotype_id_index_other, int start_site_index, int end_site_index, unordered_map<int, PBWTStack>& match_stacks, double minimum_length_centimorgans_low_resolution, double maximum_gap_length_genetic, string output_file_path_and_name, string sample_id_query, string haplotype_id_in_sample_id_query, const vector<string>& sample_ids, const vector<int>& physical_positions, const vector<double>& genetic_positions, int& total_matches_of_the_query, double& total_match_length_genetic
                            merge_on_the_fly(panel_prefix[i - 1][j], start_position, end_position, match_stacks, minimum_length_centimorgans_low_resolution, maximum_gap_length_genetic, output_file_path_and_name, sample_id_query, haplotype_id_in_sample_id_query, sample_ids, physical_positions, genetic_positions, total_matches_of_the_query, total_match_length_genetic);
                        }
                        dz.erase(panel_prefix[i - 1][j]);
                    }
                }
            }
        }
        dz.clear();
    }

    //check the genetic distance between given projected site indices
    bool isProjectedSiteIndicesAtLeastGeneticDistance(int first_site_index_projected, int second_site_index_projected, int total_number_of_sites_in_panel, int panel_window_size, double genetic_distance, const vector<double>& genetic_positions) const {
        return genetic_positions[(second_site_index_projected + 1) * panel_window_size - 1 < total_number_of_sites_in_panel ? (second_site_index_projected + 1) * panel_window_size - 1 : total_number_of_sites_in_panel - 1] - genetic_positions[first_site_index_projected * panel_window_size] >= genetic_distance;
    }
};
