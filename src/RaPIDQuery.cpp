//  * --------------------------------------------------------------------------------------------------------
//  * Name: RaPIDQuery.cpp
//  * Description: Run RaPID query long match with refined original resolution by merging on the fly approach. 
//  * Author: Yuan Wei 
//  * Created on: Jun 22, 2021
//  * Modified on: Jul 04, 2024
//  * --------------------------------------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <map>
#include <random>
#include "PBWTPanel.cpp"
using namespace std;

unsigned long seed_selected = 199371;
default_random_engine my_random_engine(seed_selected);

//get uniform distributed random integer with given range [range_from, range_to] (both inclusive)
int get_random_integer(int range_from, int range_to){
    if (range_from > range_to){
        throw invalid_argument("invalid range");
    }
    uniform_int_distribution<int> my_uniform_distribution(range_from, range_to);
    return my_uniform_distribution(my_random_engine);
}

//get a list of site indices randomly (one site per window) with the weights on minor allele frequency of the site (per site in panel)
vector<int> get_site_indices_weighted_on_minor_alleles(int panel_window_size, int total_number_of_sites_in_panel, const vector<int>& numbers_of_minor_sites){
    int number_of_windows = total_number_of_sites_in_panel % panel_window_size == 0 ? total_number_of_sites_in_panel / panel_window_size : (total_number_of_sites_in_panel / panel_window_size) + 1;
    vector<int> random_indices;
    for (int i = 0; i < number_of_windows; i++){
        //select index within the window randomly by weight
        int random_index;
        int total_number_of_minor_sites_in_window = 0;
        vector<int> weight_of_sites;
        
        //calculate the total weight of the minor allele counts for each window
        if (i < number_of_windows - 1 || (i == number_of_windows - 1 && total_number_of_sites_in_panel % panel_window_size == 0)){
            weight_of_sites.resize(panel_window_size);
            for (int j = i * panel_window_size; j <= (i + 1) * panel_window_size - 1; j++){
                total_number_of_minor_sites_in_window += numbers_of_minor_sites[j];
                weight_of_sites[j - i * panel_window_size] = numbers_of_minor_sites[j];
            }
        }
        else {
            weight_of_sites.resize(total_number_of_sites_in_panel - i * panel_window_size);
            for (int j = i * panel_window_size; j <= total_number_of_sites_in_panel - 1; j++){
                total_number_of_minor_sites_in_window += numbers_of_minor_sites[j];
                weight_of_sites[j - i * panel_window_size] = numbers_of_minor_sites[j];
            }
        }

        //add a random weight to the weight of sites and track the site index having the maximum weight
        tuple<int, int> weight_of_sites_max(0, 0); //(site index, site weight)
        for (int j = 0; j < weight_of_sites.size(); j++){
            //add a random weight to the weight of site (random number from range [0, total_number_of_minor_sites_in_window])
            int random_weight = get_random_integer(0, total_number_of_minor_sites_in_window);
            weight_of_sites[j] += random_weight;

            //track the site index having the maximum weight
            if (get<1>(weight_of_sites_max) < weight_of_sites[j]){
                get<0>(weight_of_sites_max) = j;
                get<1>(weight_of_sites_max) = weight_of_sites[j];
            }
        }

        //select the site index in window with the largest weight and add it to the list
        random_index = get<0>(weight_of_sites_max) + i * panel_window_size;
        random_indices.push_back(random_index);
    }
    return random_indices;
}

//generate subpanel by selecting sites based on random site indices
vector<vector<bool>> get_subpanel(const vector<int>& site_indices, const vector<vector<bool>>& all_sites){
    vector<vector<bool>> all_sites_subpanel;
    for (int i = 0; i < site_indices.size(); i++){
        all_sites_subpanel.push_back(all_sites[site_indices[i]]);
    }
    return all_sites_subpanel;
}

//report matches
void report_match(string output_file_path_and_name, string sample_id_query, string haplotype_id_in_sample_id_query, int haplotype_id_index_other, int start_site_index, int end_site_index, const vector<string>& sample_ids, const vector<int>& physical_positions, const vector<double>& genetic_positions, int& total_matches_of_the_query, double& total_match_length_genetic){
    //convert match individual id and haplotype id format
    int sample_id_index_other = haplotype_id_index_other / 2;
    int haplotype_id_in_sample_id_other = haplotype_id_index_other % 2; //zero is the first haplotype id and one is the second haplotype id of the sample id
    string sample_id_other = sample_ids[sample_id_index_other];

    // //when output IBD matches, excluding the individual id which is the query itself
    // if (sample_id_query == sample_id_other){
    //     return;
    // }

    //collect match data
    int start_position_physical = physical_positions[start_site_index];
    int end_position_physical = physical_positions[end_site_index];
    double start_position_genetic = genetic_positions[start_site_index];
    double end_position_genetic = genetic_positions[end_site_index];
    double match_length_genetic = end_position_genetic - start_position_genetic;
    total_matches_of_the_query += 1;
    total_match_length_genetic += match_length_genetic;

    //write to match file
    ofstream output_onefile_matches;
    output_onefile_matches.open(output_file_path_and_name, ios::app);
    if (!output_onefile_matches){
        cout << "cannot create or open file " + output_file_path_and_name << endl;
        exit(1);
    }
    if (output_onefile_matches.is_open()){
        //format: Query individual one id, Query individual one id's haplotype id, Query individual two id, Query individual two id's haplotype id, Physical start position, Physical end position, Match genetic length, Site start index, Site end index
        string reported_matches_current_query = sample_id_query + "," + haplotype_id_in_sample_id_query + "," + sample_id_other + "," + to_string(haplotype_id_in_sample_id_other) + "," + to_string(start_position_physical) + "," + to_string(end_position_physical) + "," + to_string(match_length_genetic) + "," + to_string(start_site_index) + "," + to_string(end_site_index) + "\n";
        output_onefile_matches << reported_matches_current_query;
    }
    output_onefile_matches.close();
    return;
}

//merge current detected high resolution IBD with low resolution IBDs in PBWT stack on the fly
void merge_on_the_fly(int haplotype_id_index_other, int start_site_index, int end_site_index, unordered_map<int, PBWTStack>& match_stacks, double minimum_length_centimorgans_low_resolution, double maximum_gap_length_genetic, string output_file_path_and_name, string sample_id_query, string haplotype_id_in_sample_id_query, const vector<string>& sample_ids, const vector<int>& physical_positions, const vector<double>& genetic_positions, int& total_matches_of_the_query, double& total_match_length_genetic){
    if (match_stacks.count(haplotype_id_index_other) <= 0){
        //discard the detected high resolution IBD as no correlated low resolution IBDs exist
        return;
    }
    else if (match_stacks[haplotype_id_index_other].matches.empty()){
        //discard the detected high resolution IBD as no correlated low resolution IBDs exist; also remove the entry from the map
        match_stacks.erase(haplotype_id_index_other);
        return;
    }
    else {
        //merge
        tuple<int, int> match_low = match_stacks[haplotype_id_index_other].matches.top();
        int start_site_index_low = get<0>(match_low);
        int end_site_index_low = get<1>(match_low);
        tuple<int, int> last_match_unreported = match_stacks[haplotype_id_index_other].ur_match;
        int start_site_index_unreported = get<0>(last_match_unreported);
        int end_site_index_unreported = get<1>(last_match_unreported);
        if (end_site_index_low < start_site_index){
            //case 1
            if (start_site_index_unreported != -1 && end_site_index_unreported != -1){
                if (genetic_positions[end_site_index_unreported] - genetic_positions[start_site_index_unreported] >= minimum_length_centimorgans_low_resolution){
                    report_match(output_file_path_and_name, sample_id_query, haplotype_id_in_sample_id_query, haplotype_id_index_other, start_site_index_unreported, end_site_index_unreported, sample_ids, physical_positions, genetic_positions, total_matches_of_the_query, total_match_length_genetic);
                }
                match_stacks[haplotype_id_index_other].ur_match = tuple<int, int>(-1, -1);
                match_stacks[haplotype_id_index_other].matches.pop();
            }
            while (!match_stacks[haplotype_id_index_other].matches.empty() && get<1>(match_stacks[haplotype_id_index_other].matches.top()) < start_site_index){ //new_end_site_index_low < start_site_index
                match_stacks[haplotype_id_index_other].matches.pop();
            }
            if (match_stacks[haplotype_id_index_other].matches.empty()){
                //discard the detected high resolution IBD as no correlated low resolution IBDs exist; also remove the entry from the map
                match_stacks.erase(haplotype_id_index_other);
                return;
            }
            else {
                //refresh both low resolution and unreported matches
                match_low = match_stacks[haplotype_id_index_other].matches.top();
                start_site_index_low = get<0>(match_low);
                end_site_index_low = get<1>(match_low);
                tuple<int, int> last_match_unreported = match_stacks[haplotype_id_index_other].ur_match;
                start_site_index_unreported = get<0>(last_match_unreported);
                end_site_index_unreported = get<1>(last_match_unreported);
            }
        }
        if (start_site_index_unreported == -1 || end_site_index_unreported == -1){
            //no low resolution match segment overlapped with previous high resolution match was marked but not output
            if (end_site_index_low >= start_site_index && end_site_index_low <= end_site_index){
                if (start_site_index_low < start_site_index){
                    //case 2
                    match_stacks[haplotype_id_index_other].matches.pop();
                    match_stacks[haplotype_id_index_other].ur_match = tuple<int, int>(start_site_index, end_site_index_low);
                    if (genetic_positions[end_site_index_low] - genetic_positions[start_site_index] >= minimum_length_centimorgans_low_resolution){
                        report_match(output_file_path_and_name, sample_id_query, haplotype_id_in_sample_id_query, haplotype_id_index_other, start_site_index, end_site_index_low, sample_ids, physical_positions, genetic_positions, total_matches_of_the_query, total_match_length_genetic);
                    }
                    match_stacks[haplotype_id_index_other].ur_match = tuple<int, int>(-1, -1);
                }
                else {
                    //start_site_index_low >= start_site_index
                    //case 3
                    match_stacks[haplotype_id_index_other].matches.pop();
                    match_stacks[haplotype_id_index_other].ur_match = tuple<int, int>(start_site_index_low, end_site_index_low);
                    if (genetic_positions[end_site_index_low] - genetic_positions[start_site_index_low] >= minimum_length_centimorgans_low_resolution){
                        report_match(output_file_path_and_name, sample_id_query, haplotype_id_in_sample_id_query, haplotype_id_index_other, start_site_index_low, end_site_index_low, sample_ids, physical_positions, genetic_positions, total_matches_of_the_query, total_match_length_genetic);
                    }
                    match_stacks[haplotype_id_index_other].ur_match = tuple<int, int>(-1, -1);
                }
            }
            else {
                //end_site_index_low > end_site_index
                if (start_site_index_low < start_site_index){
                    //case 5
                    match_stacks[haplotype_id_index_other].ur_match = tuple<int, int>(start_site_index, end_site_index);
                }
                else if (start_site_index_low >= start_site_index && start_site_index_low <= end_site_index){
                    //case 4
                    match_stacks[haplotype_id_index_other].ur_match = tuple<int, int>(start_site_index_low, end_site_index);
                }
                else {
                    //start_site_index_low > end_site_index
                    //case 6, discard currently detected high resolution match
                }
            }
        }
        else {
            //a low resolution match segment overlapped with previous high resolution match was marked but not output
            if (end_site_index_low >= start_site_index && end_site_index_low <= end_site_index){
                //case 2
                if (start_site_index_low < start_site_index){
                    match_stacks[haplotype_id_index_other].matches.pop();
                    if (genetic_positions[start_site_index] - genetic_positions[end_site_index_unreported] <= maximum_gap_length_genetic){
                        //extend
                        match_stacks[haplotype_id_index_other].ur_match = tuple<int, int>(start_site_index_unreported, end_site_index_low);
                        if (genetic_positions[end_site_index_low] - genetic_positions[start_site_index_unreported] >= minimum_length_centimorgans_low_resolution){
                            report_match(output_file_path_and_name, sample_id_query, haplotype_id_in_sample_id_query, haplotype_id_index_other, start_site_index_unreported, end_site_index_low, sample_ids, physical_positions, genetic_positions, total_matches_of_the_query, total_match_length_genetic);
                        }
                        match_stacks[haplotype_id_index_other].ur_match = tuple<int, int>(-1, -1);
                    }
                    else {
                        //output unreported segment
                        if (genetic_positions[end_site_index_unreported] - genetic_positions[start_site_index_unreported] >= minimum_length_centimorgans_low_resolution){
                            report_match(output_file_path_and_name, sample_id_query, haplotype_id_in_sample_id_query, haplotype_id_index_other, start_site_index_unreported, end_site_index_unreported, sample_ids, physical_positions, genetic_positions, total_matches_of_the_query, total_match_length_genetic);
                        }
                        //reset unreported segment and output
                        match_stacks[haplotype_id_index_other].ur_match = tuple<int, int>(start_site_index, end_site_index_low);
                        if (genetic_positions[end_site_index_low] - genetic_positions[start_site_index] >= minimum_length_centimorgans_low_resolution){
                            report_match(output_file_path_and_name, sample_id_query, haplotype_id_in_sample_id_query, haplotype_id_index_other, start_site_index, end_site_index_low, sample_ids, physical_positions, genetic_positions, total_matches_of_the_query, total_match_length_genetic);
                        }
                        match_stacks[haplotype_id_index_other].ur_match = tuple<int, int>(-1, -1);
                    }
                }
                else {
                    //start_site_index_low >= start_site_index
                    //case 3, not possible
                }
            }
            else {
                //end_site_index_low > end_site_index
                if (start_site_index_low < start_site_index){
                    //case 5
                    if (genetic_positions[start_site_index] - genetic_positions[end_site_index_unreported] <= maximum_gap_length_genetic){
                        //extend
                        match_stacks[haplotype_id_index_other].ur_match = tuple<int, int>(start_site_index_unreported, end_site_index);
                    }
                    else {
                        //output unreported segment
                        if (genetic_positions[end_site_index_unreported] - genetic_positions[start_site_index_unreported] >= minimum_length_centimorgans_low_resolution){
                            report_match(output_file_path_and_name, sample_id_query, haplotype_id_in_sample_id_query, haplotype_id_index_other, start_site_index_unreported, end_site_index_unreported, sample_ids, physical_positions, genetic_positions, total_matches_of_the_query, total_match_length_genetic);
                        }
                        //reset unreported segment
                        match_stacks[haplotype_id_index_other].ur_match = tuple<int, int>(start_site_index, end_site_index);
                    }
                }
                else if (start_site_index_low >= start_site_index && start_site_index_low <= end_site_index){
                    //case 4, not possible
                }
                else {
                    //start_site_index_low > end_site_index
                    //case 6, not possible
                }
            }
        }
        return;
    }
}

//split a string into tokens with given delimiter
vector<string> split_string(string line, char delimiter){
    vector<string> tokens;
    string token;
    stringstream line_stream(line);
    while (getline(line_stream, token, delimiter)){
        if (!token.empty()){
            tokens.push_back(token);
        }
    }
    return tokens;
}

static void show_usage(string program_name){
    cout << "Usage: " << program_name << " <Option(s)>\n";
    cout << "Option(s):\n";
    cout << "\t-h,--help\t\t\t\t\t\t\t\t\tShow this help message\n";
    cout << "\t-p,--panel <INPUT PANEL FILE>\t\t\t\t\t\t\tInput panel path and file name\n";
    cout << "\t-q,--query <INPUT QUERY FILE>\t\t\t\t\t\t\tInput query path and file name\n";
    cout << "\t-g,--genetic <INPUT GENETIC MAPPING FILE>\t\t\t\t\tInput genetic mapping path and file name\n";
    cout << "\t-m,--match <OUTPUT MATCH FILE>\t\t\t\t\t\tOutput match path and file name\n";
    cout << "\t-lm,--length_marker <MINIMUM IBD MARKERS IN NUMBER OF SITES>\t\t\tMinimum IBD markers in number of sites\n";
    cout << "\t-lmh,--length_marker_high <MINIMUM IBD MARKERS IN NUMBER OF SITES (HIGH)>\tMinimum IBD markers in number of sites for high resolution\n";
    cout << "\t-d,--distance <MINIMUM IBD LENGTH IN CM (LOW)>\t\t\t\t\tMinimum IBD length in centiMorgans for low resolution\n";
    cout << "\t-dh,--distance_high <MINIMUM IBD LENGTH IN CM (HIGH)>\t\t\t\tMinimum IBD length in centiMorgans for high resolution\n";
    cout << "\t-dg,--distance_gap <MAXIMUM GAP BETWEEN IBDS IN CM>\t\t\t\tMaximum gap between IBDs in centiMorgans\n";
    cout << "\t-w,--window <WINDOW SIZE>\t\t\t\t\t\t\tPanel window size for sub sampling\n";
    cout << "\t-r,--run <NUMBER OF RUNS>\t\t\t\t\t\t\tNumber of runs (sub-panels)\n";
    cout << "\t-c,--count <MINIMUM NUMBER OF COUNT OF SUCCESSES>\t\t\t\tMinimum number of count of successes (hits) as a real match\n";
}

int main(int argc, char *argv[]){
    try {
        cout << "start program" << endl;

        //variables
        double minimum_length_centimorgans = 7.0; //minimum IBD length in centiMorgans (for low resolution), also as the target IBD length
        int minimum_number_of_markers = 700; //the minimum number of markers (sites) an IBD must have
        int panel_window_size = 13; //panel window size for sub sampling (select one site value out of the given number of consecutive sites in a window as down sampling)
        int number_of_subpanels = 5; //number of runs (each run correlates to a subpanel)
        int count_of_match_success = 1; //minimum number of count of successes (hits) as a real match
        double minimum_length_centimorgans_high_resolution = 1.0; //minimum IBD length in centiMorgans (for high resolution)
        int minimum_number_of_markers_high_resolution = 100; //the minimum number of markers (sites) an IBD must have
        double maximum_gap_length_genetic = 2.0; //maximum length as a gap between two IBDs (used in match trimming; not related to random projection)

        char vcf_delimiter = '\t'; //VCF format
        char map_delimiter = '\t'; //HapMap format
        map<int, double> genetic_maps; //genetic map read from input HapMap data; key: physical position (base pairs), value: genetic position (centiMorgans)
        vector<double> panel_genetic_positions;

        vector<vector<bool>> query_sites;
        int number_of_query_individuals = 0;
        int number_of_query_sites = 0;
        vector<string> query_individual_ids;
        vector<int> query_physical_positions;

        vector<vector<bool>> panel_sites;
        int number_of_panel_individuals = 0;
        int number_of_panel_sites = 0;
        vector<string> panel_individual_ids;
        vector<int> panel_physical_positions;
        vector<int> numbers_of_minor_sites_in_panel; //total number of haplotypes having minor allele values of the site (for each site in panel); the small the minor allele frequency (the more chance all haplotypes having the same allele value), the less chance of getting selected in downsampling

        string input_panel_path_and_file_name;
        string input_query_path_and_file_name;
        string input_map_path_and_file_name;
        string output_match_path_and_file_name;

        //get command line arguments
        string program_arguments = "#Program Arguments:";
        for (int i = 1; i < argc; i++){
            string argument = argv[i];
            if ((argument == "-h") || (argument == "--help")){
                show_usage(argv[0]);
                cout << "end program" << endl;
                return 0;
            }
            else if ((argument == "-p") || (argument == "--panel")){
                if (i + 1 < argc){
                    input_panel_path_and_file_name = argv[i + 1];
                    program_arguments += " p=" + input_panel_path_and_file_name;
                    i++;
                }
                else {
                    cout << "-p (or --panel) option requires one argument\n";
                    cout << "end program" << endl;
                    return 1;
                }
            }
            else if ((argument == "-q") || (argument == "--query")){
                if (i + 1 < argc){
                    input_query_path_and_file_name = argv[i + 1];
                    program_arguments += " q=" + input_query_path_and_file_name;
                    i++;
                }
                else {
                    cout << "-q (or --query) option requires one argument\n";
                    cout << "end program" << endl;
                    return 1;
                }
            }
            else if ((argument == "-g") || (argument == "--genetic")){
                if (i + 1 < argc){
                    input_map_path_and_file_name = argv[i + 1];
                    program_arguments += " g=" + input_map_path_and_file_name;
                    i++;
                }
                else {
                    cout << "-g (or --genetic) option requires one argument\n";
                    cout << "end program" << endl;
                    return 1;
                }
            }
            else if ((argument == "-m") || (argument == "--match")){
                if (i + 1 < argc){
                    output_match_path_and_file_name = argv[i + 1];
                    program_arguments += " m=" + output_match_path_and_file_name;
                    i++;
                }
                else {
                    cout << "-m (or --match) option requires one argument\n";
                    cout << "end program" << endl;
                    return 1;
                }
            }
            else if ((argument == "-lmh") || (argument == "--length_marker_high")){
                if (i + 1 < argc){
                    minimum_number_of_markers_high_resolution = stoi(argv[i + 1]);
                    program_arguments += " lmh=" + to_string(minimum_number_of_markers_high_resolution);
                    i++;
                }
                else {
                    cout << "-lmh (or --length_marker_high) option requires one argument\n";
                    cout << "end program" << endl;
                    return 1;
                }
            }
            else if ((argument == "-lm") || (argument == "--length_marker")){
                if (i + 1 < argc){
                    minimum_number_of_markers = stoi(argv[i + 1]);
                    program_arguments += " lm=" + to_string(minimum_number_of_markers);
                    i++;
                }
                else {
                    cout << "-l (or --length) option requires one argument\n";
                    cout << "end program" << endl;
                    return 1;
                }
            }
            else if ((argument == "-d") || (argument == "--distance")){
                if (i + 1 < argc){
                    minimum_length_centimorgans = stod(argv[i + 1]);
                    program_arguments += " d=" + to_string(minimum_length_centimorgans);
                    i++;
                }
                else {
                    cout << "-d (or --distance) option requires one argument\n";
                    cout << "end program" << endl;
                    return 1;
                }
            }
            else if ((argument == "-dh") || (argument == "--distance_high")){
                if (i + 1 < argc){
                    minimum_length_centimorgans_high_resolution = stod(argv[i + 1]);
                    program_arguments += " dh=" + to_string(minimum_length_centimorgans_high_resolution);
                    i++;
                }
                else {
                    cout << "-dh (or --distance_high) option requires one argument\n";
                    cout << "end program" << endl;
                    return 1;
                }
            }
            else if ((argument == "-dg") || (argument == "--distance_gap")){
                if (i + 1 < argc){
                    maximum_gap_length_genetic = stod(argv[i + 1]);
                    program_arguments += " dg=" + to_string(maximum_gap_length_genetic);
                    i++;
                }
                else {
                    cout << "-dg (or --distance_gap) option requires one argument\n";
                    cout << "end program" << endl;
                    return 1;
                }
            }
            else if ((argument == "-w") || (argument == "--window")){
                if (i + 1 < argc){
                    panel_window_size = stoi(argv[i + 1]);
                    program_arguments += " w=" + to_string(panel_window_size);
                    i++;
                }
                else {
                    cout << "-w (or --window) option requires one argument\n";
                    cout << "end program" << endl;
                    return 1;
                }
            }
            else if ((argument == "-r") || (argument == "--run")){
                if (i + 1 < argc){
                    number_of_subpanels = stoi(argv[i + 1]);
                    program_arguments += " r=" + to_string(number_of_subpanels);
                    i++;
                }
                else {
                    cout << "-r (or --run) option requires one argument\n";
                    cout << "end program" << endl;
                    return 1;
                }
            }
            else if ((argument == "-c") || (argument == "--count")){
                if (i + 1 < argc){
                    count_of_match_success = stoi(argv[i + 1]);
                    program_arguments += " c=" + to_string(count_of_match_success);
                    i++;
                }
                else {
                    cout << "-c (or --count) option requires one argument\n";
                    cout << "end program" << endl;
                    return 1;
                }
            }
            else {
                cout << "unrecognized arguments: " << argument << endl;
                show_usage(argv[0]);
                cout << "end program" << endl;
                return 1;
            }
        }

        cout << "program arguments: " << (program_arguments.size() > 20 ? program_arguments.substr(20) : "") << endl;

        //vcf file should only have: #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT (totally 9 fields) before the first individual id
        int chromosome_id_location_on_file = 0;
        int position_location_on_file = 1;
        int id_location_on_file = 2;
        int reference_base_location_on_file = 3;
        int alternate_base_location_on_file = 4;
        int quality_location_on_file = 5;
        int filter_location_on_file = 6;
        int information_location_on_file = 7;
        int format_location_on_file = 8;
        int individual_id_start_position = 9;

        //calculate the projected minimum long match length (if l < w, l_projected = 1, which below calculation also satisfies)
        int minimum_number_of_markers_projected = minimum_number_of_markers % panel_window_size == 0 ? minimum_number_of_markers / panel_window_size : ((minimum_number_of_markers / panel_window_size) + 1);
        double minimum_length_centimorgans_projected = minimum_length_centimorgans; //since random projection only selecting site indices, the difference between genetic positions of sites does not change

        //input and output
        ifstream input_file_data;
        ofstream output_file_data;

        cout << "start read panel data" << endl;

        //read panel data from file
        input_file_data.open(input_panel_path_and_file_name);
        if (!input_file_data){
            cout << "cannot open file " + input_panel_path_and_file_name << endl;
            exit(1);
        }
        if (input_file_data.is_open()){
            int line_number = 0;
            bool end_of_file = false;
            string line_from_file;
            string line_from_file_next;
            while (!end_of_file){
                line_from_file = line_from_file_next;
                if (!end_of_file){
                    getline(input_file_data, line_from_file_next);
                    if (input_file_data.eof()){
                        end_of_file = true;
                    }
                }
                if (line_number == 0){
                    line_number++;
                    continue;
                }
                vector<bool> site_of_panel_individuals;
                if (line_from_file.substr(0, 2).compare("##") == 0){
                    //do nothing; skip header lines
                }
                else if (line_from_file.substr(0, 1).compare("#") == 0){
                    //get panel individual ids
                    vector<string> tokens = split_string(line_from_file, vcf_delimiter);
                    for (int i = individual_id_start_position; i < tokens.size(); i++){
                        panel_individual_ids.push_back(tokens[i]);
                    }
                    number_of_panel_individuals = tokens.size() - individual_id_start_position;
                }
                else {
                    //get panel site values and minor allele frequency of the site
                    vector<string> tokens = split_string(line_from_file, vcf_delimiter);
                    int count_of_allele_value_zero = 0;
                    int count_of_allele_value_one = 0;
                    for (int i = 0; i < tokens.size(); i++){
                        //get physical position
                        if (i == position_location_on_file){
                            panel_physical_positions.push_back(stoi(tokens[i]));
                        }
                        //get panel individual first and second haplotype site values and minor allele frequency of the site
                        if (tokens[i].substr(1, 1).compare("|") == 0){
                            int site_of_first_haplotype = stoi(tokens[i].substr(0, 1));
                            int site_of_second_haplotype = stoi(tokens[i].substr(2, 1));
                            site_of_panel_individuals.push_back(site_of_first_haplotype);
                            site_of_panel_individuals.push_back(site_of_second_haplotype);
                            //calculate minor allele frequency of the site
                            if (site_of_first_haplotype == 0){
                                count_of_allele_value_zero++;
                            }
                            else {
                                count_of_allele_value_one++;
                            }
                            if (site_of_second_haplotype == 0){
                                count_of_allele_value_zero++;
                            }
                            else {
                                count_of_allele_value_one++;
                            }
                        }
                    }
                    if (site_of_panel_individuals.size() > 0){
                        //store current site values from panel individuals
                        panel_sites.push_back(site_of_panel_individuals);
                        site_of_panel_individuals.clear();
                        //store minor allele frequency (minor allele count) of the site
                        if (count_of_allele_value_zero > count_of_allele_value_one){
                            numbers_of_minor_sites_in_panel.push_back(count_of_allele_value_one);
                        }
                        else {
                            numbers_of_minor_sites_in_panel.push_back(count_of_allele_value_zero);
                        }
                        number_of_panel_sites++;
                    }
                    else {
                        cout << "data at line " + to_string(line_number) + " is not valid" << endl;
                    }
                }
                line_number++;
            }
            input_file_data.close();
        }

        cout << "end read panel data" << endl;

        cout << "start read query data" << endl;

        if (input_panel_path_and_file_name.compare(input_query_path_and_file_name) == 0){
            //no need to read as panel and query data is the same
            query_sites = panel_sites;
            number_of_query_individuals = number_of_panel_individuals;
            number_of_query_sites = number_of_panel_sites;
            query_individual_ids = panel_individual_ids;
            query_physical_positions = panel_physical_positions;
        }
        else {
            //read query data from file
            input_file_data.open(input_query_path_and_file_name);
            if (!input_file_data){
                cout << "cannot open file " + input_query_path_and_file_name << endl;
                exit(1);
            }
            if (input_file_data.is_open()){
                int line_number = 0;
                bool end_of_file = false;
                string line_from_file;
                string line_from_file_next;
                while (!end_of_file){
                    line_from_file = line_from_file_next;
                    if (!end_of_file){
                        getline(input_file_data, line_from_file_next);
                        if (input_file_data.eof()){
                            end_of_file = true;
                        }
                    }
                    if (line_number == 0){
                        line_number++;
                        continue;
                    }
                    vector<bool> site_of_query_individuals;
                    if (line_from_file.substr(0, 2).compare("##") == 0){
                        //do nothing; skip header lines
                    }
                    else if (line_from_file.substr(0, 1).compare("#") == 0){
                        //get query individual ids
                        vector<string> tokens = split_string(line_from_file, vcf_delimiter);
                        for (int i = individual_id_start_position; i < tokens.size(); i++){
                            query_individual_ids.push_back(tokens[i]);
                        }
                        number_of_query_individuals = tokens.size() - individual_id_start_position;
                    }
                    else {
                        //get query site values
                        vector<string> tokens = split_string(line_from_file, vcf_delimiter);
                        for (int i = 0; i < tokens.size(); i++){
                            //get physical position
                            if (i == position_location_on_file){
                                query_physical_positions.push_back(stoi(tokens[i]));
                            }
                            //get query individual first and second haplotype site values
                            if (tokens[i].substr(1, 1).compare("|") == 0){
                                site_of_query_individuals.push_back(stoi(tokens[i].substr(0, 1)));
                                site_of_query_individuals.push_back(stoi(tokens[i].substr(2, 1)));
                            }
                        }
                        if (site_of_query_individuals.size() > 0){
                            //store current site values from query individuals
                            query_sites.push_back(site_of_query_individuals);
                            site_of_query_individuals.clear();
                            number_of_query_sites++;
                        }
                        else {
                            cout << "data at line " + to_string(line_number) + " is not valid" << endl;
                        }
                    }
                    line_number++;
                }
                input_file_data.close();
            }
        }

        cout << "end read query data" << endl;

        cout << "start read map data" << endl;

        //read genetic map data from file
        input_file_data.open(input_map_path_and_file_name);
        if (!input_file_data){
            cout << "cannot open file " + input_map_path_and_file_name << endl;
            exit(1);
        }
        if (input_file_data.is_open()){
            int line_number = 0;
            string line_from_file;
            while (getline(input_file_data, line_from_file)){
                //get genetic map data
                vector<string> tokens = split_string(line_from_file, map_delimiter);
                if (tokens.size() >= 4){
                    //genetic map HapMap format: #Chromosome \t Position(bp) \t Rate(cM/Mb) \t Map(cM)
                    if (line_number > 0){
                        //skip the first header line
                        int physical_location = stoi(tokens[1]);
                        double genetic_location = stod(tokens[3]);
                        genetic_maps[physical_location] = genetic_location;
                    }
                }
                line_number++;
            }
            input_file_data.close();
        }

        //interpolate genetic map based on panel physical positions
        map<int, double>::iterator position_exact_iterator, position_lower_iterator, position_upper_iterator;
        for (int i = 0; i < panel_physical_positions.size(); i++){
            double panel_genetic_position = 0.0;
            //search current physical position from genetic map
            if (!genetic_maps.empty()){
                position_exact_iterator = genetic_maps.find(panel_physical_positions[i]);
                if (position_exact_iterator == genetic_maps.end()){
                    //exact physical position not found; need to estimate genetic position with neighbor positions in genetic map
                    position_upper_iterator = genetic_maps.upper_bound(panel_physical_positions[i]);
                    if (position_upper_iterator == genetic_maps.end()){
                        if (panel_physical_positions[i] > genetic_maps.rbegin()->first){
                            //physical position is larger than any physical position in map
                            panel_genetic_position = genetic_maps.rbegin()->second;
                        }
                        else { 
                            //physical position is less than any physical position in map
                            panel_genetic_position = genetic_maps.begin()->second;
                        }
                    }
                    else {
                        if (position_upper_iterator == genetic_maps.begin()){
                            //physical position is less than any physical position in map
                            panel_genetic_position = genetic_maps.begin()->second;
                        }
                        else {
                            //estimate genetic position by linear interpolation
                            int physical_position_next = position_upper_iterator->first;
                            double genetic_position_next = position_upper_iterator->second;
                            position_lower_iterator = --position_upper_iterator;
                            int physical_position_previous = position_lower_iterator->first;
                            double genetic_position_previous = position_lower_iterator->second;
                            panel_genetic_position = genetic_position_next + (double)(panel_physical_positions[i] - physical_position_next) * (genetic_position_next - genetic_position_previous) / (double)(physical_position_next - physical_position_previous);
                        }
                    }
                }
                else {
                    //exact physical position found; just add its corrlated genetic position
                    panel_genetic_position = position_exact_iterator->second;
                }
            }
            panel_genetic_positions.push_back(panel_genetic_position);
        }

        cout << "end read map data" << endl;

        cout << "start verify data" << endl;

        bool is_data_satisfied = false;

        //verify if query physical positions are the same as panel physical positions
        if (query_physical_positions == panel_physical_positions){
            is_data_satisfied = true;
        }
        else if (query_physical_positions.size() == 0){
            is_data_satisfied = false;
            cout << "query has no data" << endl;
        }
        else if (panel_physical_positions.size() == 0){
            is_data_satisfied = false;
            cout << "panel has no data" << endl;
        }
        else {
            is_data_satisfied = false;
            cout << "query and panel have mismatched POS values" << endl;
        }

        //verify parameters
        if (panel_window_size <= 0){
            is_data_satisfied = false;
            cout << "window size is invalid" << endl;
        }
        if (number_of_subpanels <= 0){
            is_data_satisfied = false;
            cout << "number of runs is invalid" << endl;
        }
        if (count_of_match_success <= 0){
            is_data_satisfied = false;
            cout << "count of success is invalid" << endl;
        }

        cout << "end verify data" << endl;

        if (!is_data_satisfied){
            return 0;
        }
        
        cout << "start build sub panels" << endl;

        //run RaPID query single-resolution with merge on the fly refinement
        int number_of_query_haplotypes = number_of_query_individuals * 2;
        int number_of_panel_haplotypes = number_of_panel_individuals * 2;
        cout << "number of query haplotypes: " << to_string(number_of_query_haplotypes) << endl;
        cout << "number of query sites: " << to_string(number_of_query_sites) << endl;
        cout << "number of panel haplotypes: " << to_string(number_of_panel_haplotypes) << endl;
        cout << "number of panel sites: " << to_string(number_of_panel_sites) << endl;

        //prepare panels
        PBWTPanel pbwt_full_panel = PBWTPanel(number_of_panel_sites, number_of_panel_haplotypes, panel_sites);
        vector<PBWTPanel> pbwt_sub_panels; //low resolution PBWT panels for random projection runs (one per run)
        vector<vector<int>> random_site_indices_list; //low resolution panel site indices list (to build sub panels and sub queries)
        if (panel_window_size <= 1){
            panel_window_size = 1;
        }
        if (panel_window_size >= number_of_panel_sites){
            panel_window_size = number_of_panel_sites;
        }
        
        //generate PBWT sub panels in low resolution
        if (panel_window_size > 1){
            //prepare all random indices and generate all sub panels and their pbwt sub panels
            for (int j = 0; j < number_of_subpanels; j++){
                //generate random indices (for getting current sub panel and later for sub query)
                vector<int> random_site_indices = get_site_indices_weighted_on_minor_alleles(panel_window_size, number_of_panel_sites, numbers_of_minor_sites_in_panel); //weighted by minor allele frequency

                //generate sub panel by selecting sites based on random site indices
                vector<vector<bool>> panel_sites_sub_panel = get_subpanel(random_site_indices, panel_sites);

                //use sub panel to generate PBWT sub panels
                PBWTPanel pbwt_sub_panel = PBWTPanel(panel_sites_sub_panel.size(), number_of_panel_haplotypes, panel_sites_sub_panel);

                //store random site indices, sub panel, pbwt sub panel, and sub queries
                random_site_indices_list.push_back(random_site_indices);
                pbwt_sub_panels.push_back(pbwt_sub_panel);
            }
        }

        cout << "end build sub panels" << endl;

        cout << "start run queries" << endl;

        //write match header line
        output_file_data.open(output_match_path_and_file_name, ios::app);
        if (!output_file_data){
            cout << "cannot create or open file " + output_match_path_and_file_name << endl;
            exit(1);
        }
        if (output_file_data.is_open()){
            string reported_matches = "#Query individual id, Query individual haplotype id, Panel individual id, Panel individual haplotype id, Physical start position, Physical end position, Genetic length, Site start index, Site end index\n";
            output_file_data << reported_matches;
        }
        output_file_data.close();

        //get matches for all queries
        if (panel_window_size == 1){
            //only run high resolution and no refinement is needed
            for (int i = 0; i < number_of_query_haplotypes; i++){
                //get matches from high resolution (full) panel
                PBWTBucket bucket_high = PBWTBucket(number_of_panel_sites);
                
                //get query
                vector<bool> query_site;
                for (int j = 0 ; j < number_of_query_sites; j++){
                    query_site.push_back(query_sites[j][i]);
                }

                //find long matches for the query; returned tuple is: haplotype_id_index_in_panel, starting_position_(inclusive), ending_position_(inclusive)
                vector<tuple<int, int, int>> query_long_matches_high = pbwt_full_panel.getQueryLongMatchesWithGeneticDistance(query_site, minimum_length_centimorgans_high_resolution, minimum_number_of_markers_high_resolution, panel_genetic_positions, 1);

                //convert query individual data
                int query_individual_haplotype_id_index = i;
                int query_individual_id_index = query_individual_haplotype_id_index / 2;
                int query_individual_haplotype_id = query_individual_haplotype_id_index % 2; //zero is the first haplotype id and one is the second haplotype id of the query individual
                string query_individual_id = query_individual_ids[query_individual_id_index];
                int total_matches_of_the_query = 0;
                double total_match_length_genetic = 0.0;

                //write match data
                output_file_data.open(output_match_path_and_file_name, ios::app);
                if (!output_file_data){
                    cout << "cannot create or open file " + output_match_path_and_file_name << endl;
                    exit(1);
                }
                if (output_file_data.is_open()){
                    //gather current query result
                    for (int j = 0; j < query_long_matches_high.size(); j++){
                        int panel_individual_haplotype_id_index = get<0>(query_long_matches_high[j]);
                        int panel_individual_id_index = panel_individual_haplotype_id_index / 2;
                        int panel_individual_haplotype_id = panel_individual_haplotype_id_index % 2; //zero is the first haplotype id and one is the second haplotype id of the panel individual
                        string panel_individual_id = panel_individual_ids[panel_individual_id_index];
                        int start_site_index = get<1>(query_long_matches_high[j]);
                        int end_site_index = get<2>(query_long_matches_high[j]);
                        int start_position_physical = panel_physical_positions[start_site_index];
                        int end_position_physical = panel_physical_positions[end_site_index];
                        double start_position_genetic = panel_genetic_positions[start_site_index];
                        double end_position_genetic = panel_genetic_positions[end_site_index];
                        double match_length_genetic = end_position_genetic - start_position_genetic;
                        total_matches_of_the_query++;
                        total_match_length_genetic += match_length_genetic;

                        //format: #Query individual id, Query individual haplotype id, Panel individual id, Panel individual haplotype id, Physical start position, Physical end position, Genetic length, Site start index, Site end index
                        string reported_matches_current_query = query_individual_id + "," + to_string(query_individual_haplotype_id) + "," + panel_individual_id + "," + to_string(panel_individual_haplotype_id) + "," + to_string(start_position_physical) + "," + to_string(end_position_physical) + "," + to_string(match_length_genetic) + "," + to_string(start_site_index) + "," + to_string(end_site_index) + "\n";
                        output_file_data << reported_matches_current_query;
                    }
                }
                output_file_data.close();
            }
        }
        else {
            //run low resolution and refine afterwards
            for (int i = 0; i < number_of_query_haplotypes; i++){
                //get query (for later merge on the fly with high resolution query)
                vector<bool> query_site;
                for (int j = 0 ; j < number_of_query_sites; j++){
                    query_site.push_back(query_sites[j][i]);
                }

                //get matches from low resolution sub panels
                PBWTBucket bucket_low = PBWTBucket(number_of_panel_sites);
                
                //generate sub query
                vector<vector<bool>> query_site_list;
                for (int j = 0; j < random_site_indices_list.size(); j++){
                    vector<bool> query_site_temp;
                    for (int j1 = 0; j1 < random_site_indices_list[j].size(); j1++){
                        query_site_temp.push_back(query_sites[random_site_indices_list[j][j1]][i]);
                    }
                    query_site_list.push_back(query_site_temp);
                }

                //perform random projection runs on sub queries sequentially
                for (int j = 0; j < number_of_subpanels; j++){
                    pbwt_sub_panels[j].getQueryLongMatchesWithGeneticDistance(query_site_list[j], minimum_length_centimorgans_projected, minimum_number_of_markers_projected, panel_genetic_positions, panel_window_size, bucket_low);
                }
                bucket_low.removeAllFalseMatchesAndMerge((panel_window_size <= 1 ? 1 : count_of_match_success), panel_window_size);

                //convert query individual data
                int query_individual_haplotype_id_index = i;
                int query_individual_id_index = query_individual_haplotype_id_index / 2;
                int query_individual_haplotype_id = query_individual_haplotype_id_index % 2; //zero is the first haplotype id and one is the second haplotype id of the query individual
                string query_individual_id = query_individual_ids[query_individual_id_index];
                int total_matches_of_the_query = 0;
                double total_match_length_genetic = 0.0;

                //convert low resolution PBWT bucket to PBWT stack
                unordered_map<int, PBWTStack> match_stacks;

                //iterate all low resolution matches and check each of them (from site n - 1 to site 0)
                unordered_map<int, PBWTMatches>::iterator matches_iterator;
                unordered_map<int, PBWTMatchField>::iterator match_field_iterator;
                for (int i = bucket_low.matches.size() - 1; i >= 0; i = bucket_low.getPreviousEndPositionLocation(i, panel_window_size)){
                    for (matches_iterator = bucket_low.matches[i].begin(); matches_iterator != bucket_low.matches[i].end(); matches_iterator++){
                        int panel_individual_haplotype_id_index = matches_iterator->first;
                        int panel_individual_id_index = panel_individual_haplotype_id_index / 2;
                        int panel_individual_haplotype_id = panel_individual_haplotype_id_index % 2; //zero is the first haplotype id and one is the second haplotype id of the panel individual
                        string panel_individual_id = panel_individual_ids[panel_individual_id_index];
                        for (match_field_iterator = (matches_iterator->second).fields_of_matches.begin(); match_field_iterator != (matches_iterator->second).fields_of_matches.end(); match_field_iterator++){
                            int start_site_index = match_field_iterator->first;
                            int end_site_index = i;

                            //store match
                            unordered_map<int, PBWTStack>::iterator match_stacks_iterator = match_stacks.find(panel_individual_haplotype_id_index);
                            if (match_stacks_iterator != match_stacks.end()){
                                //haplotype id index exists in current map
                                match_stacks[panel_individual_haplotype_id_index].matches.push(tuple<int, int>(start_site_index, end_site_index));
                            }
                            else {
                                //haplotype id index does not exist
                                match_stacks[panel_individual_haplotype_id_index] = PBWTStack();
                                match_stacks[panel_individual_haplotype_id_index].matches.push(tuple<int, int>(start_site_index, end_site_index));
                            }
                        }
                    }
                }

                //refine low resolution matches by merging them on the fly when querying full resolution panel
                pbwt_full_panel.getQueryLongMatchesWithGeneticDistanceWhileMerge(query_individual_id, to_string(query_individual_haplotype_id), query_site, minimum_length_centimorgans_high_resolution, minimum_number_of_markers_high_resolution, maximum_gap_length_genetic, minimum_length_centimorgans, panel_individual_ids, panel_physical_positions, panel_genetic_positions, 1, output_match_path_and_file_name, match_stacks, total_matches_of_the_query, total_match_length_genetic, merge_on_the_fly);

                //check the unreported match field in PBWT stack after completion of merging attempt with the last detected high resolution IBD
                unordered_map<int, PBWTStack>::iterator match_stacks_iterator;
                for (match_stacks_iterator = match_stacks.begin(); match_stacks_iterator != match_stacks.end();){
                    int panel_individual_haplotype_id_index = match_stacks_iterator->first;

                    //check unreported match field before discard
                    tuple<int, int> match_unreported = match_stacks[panel_individual_haplotype_id_index].ur_match;
                    int start_site_index_unreported = get<0>(match_unreported);
                    int end_site_index_unreported = get<1>(match_unreported);
                    if (start_site_index_unreported != -1 && end_site_index_unreported != -1){
                        //output match (panel_individual_haplotype_id_index, start_site_index_unreported, end_site_index_unreported), if its length is qualified
                        if (panel_genetic_positions[end_site_index_unreported] - panel_genetic_positions[start_site_index_unreported] >= minimum_length_centimorgans){
                            report_match(output_match_path_and_file_name, query_individual_id, to_string(query_individual_haplotype_id), panel_individual_haplotype_id_index, start_site_index_unreported, end_site_index_unreported, panel_individual_ids, panel_physical_positions, panel_genetic_positions, total_matches_of_the_query, total_match_length_genetic);
                        }
                    }
                    match_stacks[panel_individual_haplotype_id_index].ur_match = tuple<int, int>(-1, -1);

                    match_stacks_iterator++;
                    match_stacks.erase(panel_individual_haplotype_id_index);
                }
            }
        }

        cout << "end run queries" << endl;

        cout << "end program" << endl;
    }
    catch (exception& exception_output){
        cerr << exception_output.what() << endl;
    }
    return 0;
}
