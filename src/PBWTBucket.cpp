//  * --------------------------------------------------------------------------------------------------------
//  * Name: PBWTBucket.cpp
//  * Description: A data structure storing matches from PBWT algorithms.
//  * Author: Yuan Wei 
//  * Created on: Jun 22, 2021
//  * --------------------------------------------------------------------------------------------------------

#include <unordered_map>
using namespace std;

class PBWTMatchField{
    public:
    int count_of_matches; //count the total number of found matches having the same start position, end position, and sample id
    int previous_end_position; //track the last seen neighbor match's end position having the same sample id as this one's

    PBWTMatchField(){
        count_of_matches = 0;
        previous_end_position = -1;
    }
};

class PBWTMatches{
    public:
    unordered_map<int, PBWTMatchField> fields_of_matches; //dictionary: key is the start position, value is PBWTMatchField
    int start_position_min; //the minimum start position found in fields_of_matches

    PBWTMatches(){
        start_position_min = -1;
    }

    ~PBWTMatches(){
        fields_of_matches.clear();
    }
};

class PBWTBucket{
    public:
    vector<unordered_map<int, PBWTMatches>> matches; //dictionary: key is sample id (the index of one of the haplotype id of an individual, which has a region found as a match to the query), value is PBWTMatches with the sample id

    PBWTBucket(){
    }

    PBWTBucket(int total_number_of_sites_in){
        matches.resize(total_number_of_sites_in);
    }

    ~PBWTBucket(){
        for (int i = 0; i < matches.size(); i++){
            matches[i].clear();
        }
    }

    int getPreviousEndPositionLocation(int end_position_current, int window_size) const {
        int end_position_previous = 0;
        if (end_position_current == matches.size() - 1 && (matches.size() % window_size) > 0){
            end_position_previous = end_position_current - (matches.size() % window_size);
        }
        else { 
            //end_position_current < matches.size() - 1
            end_position_previous = end_position_current - window_size;
        }
        if (end_position_previous == end_position_current || end_position_previous < 0 || end_position_previous > matches.size() - 1){
            //end_position_current is the first position (i.e. window_size - 1) already
            return -1;
        }
        else {
            return end_position_previous;
        }
    }

    int getNextEndPositionLocation(int end_position_current, int window_size) const {
        //end_position_current in [window_size - 1, matches.size() - 1]
        int end_position_next = end_position_current + window_size > matches.size() - 1 ? matches.size() - 1 : end_position_current + window_size;
        if (end_position_next == end_position_current || end_position_next < 0 || end_position_next > matches.size() - 1){
            //end_position_current is the last position (i.e. matches.size() - 1) already
            return matches.size();
        }
        else {
            return end_position_next;
        }
    }

    void addMatch(int sample_id_in, int start_position_in, int end_position_in, int panel_window_size, bool site_barrier_used){
        if (start_position_in <= end_position_in && start_position_in >= 0 && end_position_in < matches.size()){
            unordered_map<int, PBWTMatches>::iterator matches_iterator = matches[end_position_in].find(sample_id_in);
            if (matches_iterator != matches[end_position_in].end()){
                //sample id exists in current end position list 

                //find if start_position_in exists
                unordered_map<int, PBWTMatchField>::iterator match_field_iterator = (matches_iterator->second).fields_of_matches.find(start_position_in);
                if (match_field_iterator != (matches_iterator->second).fields_of_matches.end()){
                    //start position exists; just update its counts
                    (match_field_iterator->second).count_of_matches++;
                }
                else {
                    //start position does not exist; add it
                    (matches_iterator->second).fields_of_matches[start_position_in] = PBWTMatchField();
                    (matches_iterator->second).fields_of_matches[start_position_in].count_of_matches = 1;
                }

                //get current minimum start position of current end position container
                int start_position_min_curr = (matches_iterator->second).start_position_min;

                //updates for matches ending at current end position
                //update count of matches for the new match
                if (start_position_in > start_position_min_curr){
                    matches[end_position_in][sample_id_in].fields_of_matches[start_position_in].count_of_matches++;
                }

                //update count of matches for other matches
                if (start_position_in <= start_position_min_curr){
                    //update other matches
                    for (match_field_iterator = (matches_iterator->second).fields_of_matches.begin(); match_field_iterator != (matches_iterator->second).fields_of_matches.end(); ++match_field_iterator){
                        if ((match_field_iterator->first) > start_position_in){
                            (match_field_iterator->second).count_of_matches++;
                        }
                    }
                    //update start_position_min
                    (matches_iterator->second).start_position_min = start_position_in;
                    start_position_min_curr = start_position_in;
                }
            }
            else {
                //sample id does not exist in current end position list
                //add new match with the new sample id for current end position list
                matches[end_position_in][sample_id_in].fields_of_matches[start_position_in] = PBWTMatchField();
                matches[end_position_in][sample_id_in].fields_of_matches[start_position_in].count_of_matches = 1;
                matches[end_position_in][sample_id_in].start_position_min = start_position_in;

                //no need to make updates for matches ending at current end position (neither update count of matches for the new match, nor update count of matches for other matches)
            }

            //updates for matches ending at less than current end position
            unordered_map<int, PBWTMatchField>::iterator match_field_iterator;
            
            if (!site_barrier_used){
                //update count of matches for the new match (this is not necessary if all subpanels are processed in parallel per site)
                for (int i = getNextEndPositionLocation(end_position_in, panel_window_size); i < matches.size(); i = getNextEndPositionLocation(i, panel_window_size)){ //for (int i = end_position_in + 1; i < matches.size(); i++)
                    matches_iterator = matches[i].find(sample_id_in);
                    if (matches_iterator != matches[i].end()){
                        for (match_field_iterator = (matches_iterator->second).fields_of_matches.begin(); match_field_iterator != (matches_iterator->second).fields_of_matches.end(); match_field_iterator++){
                            //update the new match if its start_position_in is greater or equal to the match's start position under the list
                            if (match_field_iterator->first <= start_position_in){
                                matches[end_position_in][sample_id_in].fields_of_matches[start_position_in].count_of_matches++;
                            }
                        }
                    }
                }
            }

            //update count of matches for other matches
            int i = getPreviousEndPositionLocation(end_position_in, panel_window_size);
            while (i > -1 && i >= start_position_in){
            //for (int i = end_position_in - 1; i >= start_position_in; i--)
                matches_iterator = matches[i].find(sample_id_in);
                if (matches_iterator != matches[i].end()){
                    for (match_field_iterator = (matches_iterator->second).fields_of_matches.begin(); match_field_iterator != (matches_iterator->second).fields_of_matches.end(); ++match_field_iterator){
                        //update the match under this list if its start position is greater or equal to the start_position_in
                        if (match_field_iterator->first >= start_position_in){
                            (match_field_iterator->second).count_of_matches++;
                        }
                    }
                }
                i = getPreviousEndPositionLocation(i, panel_window_size);
            }
        }
        else {
            throw invalid_argument("invalid range");
        }
    }

    void removeAllFalseMatchesAndMerge(int count_of_match_success_in, int panel_window_size){
        if (panel_window_size > 1){
            //for both removing all false matches and merging matches
            unordered_map<int, PBWTMatches>::iterator matches_iterator;
            unordered_map<int, PBWTMatchField>::iterator match_field_iterator;

            //for removing all false matches
            unordered_map<int, PBWTMatchField>::iterator match_field_update_iterator;

            //for merging matches
            unordered_map<int, PBWTMatches>::iterator matches_prior_iterator;
            unordered_map<int, PBWTMatchField>::iterator match_field_prior_iterator;
            unordered_map<int, PBWTMatchField>::iterator match_field_current_iterator;

            for (int i = panel_window_size - 1; i < matches.size(); i = getNextEndPositionLocation(i, panel_window_size)){ //for (int i = 0; i < matches.size(); i++)
                //remove all false matches
                for (matches_iterator = matches[i].begin(); matches_iterator != matches[i].end();){
                    for (match_field_iterator = (matches_iterator->second).fields_of_matches.begin(); match_field_iterator != (matches_iterator->second).fields_of_matches.end();){
                        if ((match_field_iterator->second).count_of_matches < count_of_match_success_in){
                            //remove the match
                            int start_position_curr = match_field_iterator->first;
                            match_field_iterator = (matches_iterator->second).fields_of_matches.erase(match_field_iterator);
                            //update start_position_min, if current removed match has the min start position and there are matches in directory
                            if (start_position_curr == (matches_iterator->second).start_position_min && (matches_iterator->second).fields_of_matches.size() > 0){
                                int start_position_min_curr = ((matches_iterator->second).fields_of_matches.begin())->first;
                                for (match_field_update_iterator = (matches_iterator->second).fields_of_matches.begin(); match_field_update_iterator != (matches_iterator->second).fields_of_matches.end(); match_field_update_iterator++){
                                    if (start_position_min_curr > match_field_update_iterator->first){
                                        start_position_min_curr = match_field_update_iterator->first;
                                    }
                                }
                                (matches_iterator->second).start_position_min = start_position_min_curr;
                            }
                        }
                        else {
                            match_field_iterator++;
                        }
                    }
                    //remove the sample id list if it does not contain any matches
                    if ((matches_iterator->second).fields_of_matches.size() == 0){
                        matches_iterator = matches[i].erase(matches_iterator);
                    }
                    else {
                        matches_iterator++;
                    }
                }

                //merge matches
                for (matches_iterator = matches[i].begin(); matches_iterator != matches[i].end(); matches_iterator++){
                    if (matches_iterator != matches[i].end()){
                        int sample_id_curr = matches_iterator->first;
                        int count_of_matches_curr = ((matches_iterator->second).fields_of_matches.begin()->second).count_of_matches;
                        int previous_end_position_curr = ((matches_iterator->second).fields_of_matches.begin()->second).previous_end_position;
                        int start_position_min_curr = (matches_iterator->second).start_position_min;

                        //merge matches ending at current site (only one match should remain) for sample_id_curr
                        for (match_field_iterator = (matches_iterator->second).fields_of_matches.begin(); match_field_iterator != (matches_iterator->second).fields_of_matches.end();){
                            if (match_field_iterator->first > start_position_min_curr){
                                //remove the match
                                match_field_iterator = (matches_iterator->second).fields_of_matches.erase(match_field_iterator);
                            }
                            else {
                                //keep the match
                                match_field_iterator++;
                            }
                        }

                        //merge matches ending at less than current site for sample_id_curr
                        int j = getPreviousEndPositionLocation(i, panel_window_size);
                        while (j > -1 && j >= start_position_min_curr - 1){
                            matches_prior_iterator = matches[j].find(sample_id_curr);
                            if (matches_prior_iterator != matches[j].end()){
                                int start_position_prev = (matches_prior_iterator->second).start_position_min;

                                //remove the prior match and the sample id list had the prior match as the match's start position has been stored (only one match should exist), and update the current match start position (if needed)
                                match_field_prior_iterator = (matches_prior_iterator->second).fields_of_matches.find(start_position_prev);
                                if (match_field_prior_iterator != (matches_prior_iterator->second).fields_of_matches.end()){
                                    //remove the prior match
                                    match_field_prior_iterator = (matches_prior_iterator->second).fields_of_matches.erase(match_field_prior_iterator);

                                    //remove the prior list which had prior match
                                    matches_prior_iterator = matches[j].erase(matches_prior_iterator);

                                    //update the current match with new start position if needed
                                    if (start_position_prev < start_position_min_curr){
                                        //update start position of the match (add a new match as the updated current match)
                                        (matches_iterator->second).fields_of_matches[start_position_prev] = PBWTMatchField();
                                        (matches_iterator->second).fields_of_matches[start_position_prev].count_of_matches = count_of_matches_curr;
                                        (matches_iterator->second).fields_of_matches[start_position_prev].previous_end_position = previous_end_position_curr;

                                        //update start position of the match (remove the outdated current match)
                                        match_field_current_iterator = (matches_iterator->second).fields_of_matches.find(start_position_min_curr);
                                        if (match_field_current_iterator != (matches_iterator->second).fields_of_matches.end()){
                                            match_field_current_iterator = (matches_iterator->second).fields_of_matches.erase(match_field_current_iterator);
                                            (matches_iterator->second).start_position_min = start_position_prev;
                                        }
                                        else {
                                            //should not be here as one and only one match with such start position exists
                                        }
                                        break;
                                    }
                                }
                                else {
                                    //should not be here as one and only one match with such start position exists
                                }
                            }
                            j = getPreviousEndPositionLocation(j, panel_window_size);
                        }
                    }
                }
            }
        }
    }
};
