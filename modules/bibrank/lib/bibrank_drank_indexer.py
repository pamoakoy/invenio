## -*- mode: python; coding: utf-8; -*-
##
## This file is part of Invenio.
## Copyright (C) 2007, 2008, 2010, 2011 CERN.
##
## Invenio is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 2 of the
## License, or (at your option) any later version.
##
## Invenio is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Invenio; if not, write to the Free Software Foundation, Inc.,
## 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.

"""This file defines functions/sub tasks used for processing Distributed Ranking(DRANK) in invenio. 
It includes normalization, generating of probability distributions, percentile ranks, 
interpolation(Newton's second order) and a Ranker function for aggregation.This processes are
applied to well exposed/less exposed records as well as scores from query dependent(relevance)
methods currently in invenio. Dictionaries of several look-up tables, and ranking scores using
counts/data from invenio logger(bibrank_drank_logger/queryanalyzer) are generated and stored 
in rnkDRANKLUT/rnkMETHODDATA tables. 
All configurations are currently read from drank-specific configuration file."""

__revision__ = "$Id$"

import datetime
import os, sys, cgi,time,ConfigParser,string,csv,glob
try:
    from scipy.stats.mstats import scoreatpercentile
    from scipy.stats.kde import gaussian_kde
    import_scipy=1
except ImportError:
    import_scipy=0
    
from invenio.dbquery import deserialize_via_marshal, serialize_via_marshal, run_sql
try:
    from numpy import mean, array, sort, append, divide, any, inf, column_stack, row_stack
    import_numpy=1
except ImportError:
    import_numpy=0
from invenio.config import CFG_ETCDIR,CFG_LOGDIR
from invenio.bibtask import task_init, write_message, task_set_option, \
        task_get_option, task_update_progress, task_update_status, \
        task_get_task_param, task_sleep_now_if_required
        
def get_current_time():
    """ Function generates the current time"""
    current_time = datetime.datetime.now()
    return current_time

def get_elaborate_submit_param(key, value, dummyopts, dummyargs):
    """
    Elaborate task submission parameter.  based on bibtask's
    task_submit_elaborate_specific_parameter_fnc.
    """
    if key in ('-d', '--wrd-score-folder'):
        try:
            task_set_option('wrd-score-folder', value)
        except ValueError:
            raise StandardError, "ERROR: Input '%s' is not valid." % \
                  value
    elif key in ('-c', '--insert-counts'):        
        try:
            task_set_option('insert-counts', bool(value))
        except ValueError:
            raise StandardError, "ERROR: Input '%s' is not valid." % \
                  value  

    elif key in ('-q', '--insert-quality'):        
        try:
            task_set_option('insert-quality', bool(value))
        except ValueError:
            raise StandardError, "ERROR: Input '%s' is not valid." % \
                  value  
    elif key in ('-w', '--insert-wrd'):
        try:
            task_set_option('insert-wrd', bool(value))
        except ValueError:
            raise StandardError, "ERROR: Input '%s' is not valid." % \
                  value                  
    elif key in ('-f', '--insert-fresh'):
        try:
            task_set_option('insert-fresh', bool(value))
        except ValueError:
            raise StandardError, "ERROR: Input '%s' is not valid." % \
                  value              
    else:
        return False
    return True
def get_config(dictname):
    """Read all drank specific configuration files and return as a tuple """
  
    drank_methods = {}
#    Trying reading the main config file.
    try:
        file = CFG_ETCDIR + "/bibrank/" + dictname + ".cfg"
#        file= "/home/pette/src/invenio/modules/bibrank/etc/" + dictname + ".cfg"
        config = ConfigParser.ConfigParser()
        config.readfp(open(file))
    except StandardError, e:
        pass
    #        cfg_parameters=config.get("parameters","")
    drank_bibrank_function = config.get("rank_method", "function")
    drank_methods["function"] = drank_bibrank_function
    if config.has_section(drank_bibrank_function): 
        drank_methods[drank_bibrank_function] = {}  
        drank_methods[drank_bibrank_function]["prefix"] = config.get(drank_bibrank_function, "relevance_number_output_prologue")
        drank_methods[drank_bibrank_function]["postfix"] = config.get(drank_bibrank_function, "relevance_number_output_epilogue")
        drank_methods[drank_bibrank_function]["ranked_by"] = config.get(drank_bibrank_function, "ranked_by")
        drank_methods[drank_bibrank_function]["description"] = config.get(drank_bibrank_function, "description")
        
##  Get DRANK weights and parameter         
    if config.has_section("drank_parameters"): 
        drank_methods["drank_parameters"] = {}           
        drank_methods["drank_parameters"]["drank_lut_table"] = config.get("drank_parameters", "drank_lut_table")
        drank_methods["drank_parameters"]["counts_table"]=config.get("drank_parameters", "counts_table")
        drank_methods["drank_parameters"]["max_nbr_of_displays"] = config.getfloat("drank_parameters", "max_nbr_of_displays")        
        drank_methods["drank_parameters"]["weight_of_relevance"] = config.getfloat("drank_parameters", "weight_of_relevance")
        drank_methods["drank_parameters"]["weight_of_exposed"] = config.getfloat("drank_parameters", "weight_of_exposed")
        drank_methods["drank_parameters"]["percentage_of_freshness"] = config.getfloat("drank_parameters", "percentage_of_freshness")
        drank_methods["drank_parameters"]["min_freshness_value"] = config.getfloat("drank_parameters", "min_freshness_value")      

#    Read DRANK specific Word Similarity parameters
    if config.has_section("relevance_parameters"):
        drank_methods["relevance_parameters"]={}
        drank_methods["relevance_parameters"]["max_word_occurence"] = float(config.get("relevance_parameters", "max_word_occurence"))
        drank_methods["relevance_parameters"]["min_word_occurence"] = float(config.get("relevance_parameters", "min_word_occurence"))
        drank_methods["relevance_parameters"]["min_word_length"] = int(config.get("relevance_parameters", "min_word_length"))
        drank_methods["relevance_parameters"]["min_nr_words_docs"] = int(config.get("relevance_parameters", "min_nr_words_docs"))
        drank_methods["relevance_parameters"]["max_nr_words_upper"] = int(config.get("relevance_parameters", "max_nr_words_upper"))
        drank_methods["relevance_parameters"]["max_nr_words_lower"] = int(config.get("relevance_parameters", "max_nr_words_lower"))
        drank_methods["relevance_parameters"]["default_min_relevance"] = int(config.get("relevance_parameters", "default_min_relevance"))


#  Extract all the individual weights/coefficient for the ranker function.
#  At the moment the ranker is just a linear combination of all normalized individually distributed ranking scores       
    if config.has_section("drank_coefficients"):    
        drank_methods["drank_coefficients"] = {}        
        drank_methods["drank_coefficients"]["coefficient_nbr_downloads"] = config.getfloat("drank_coefficients", "coefficient_nbr_downloads")
        drank_methods["drank_coefficients"]["coefficient_nbr_displays"] = config.getfloat("drank_coefficients", "coefficient_nbr_displays")
        drank_methods["drank_coefficients"]["coefficient_nbr_seens"] = config.getfloat("drank_coefficients", "coefficient_nbr_seens")
        drank_methods["drank_coefficients"]["coefficient_nbr_views"] = config.getfloat("drank_coefficients", "coefficient_nbr_views")
        drank_methods["drank_coefficients"]["coefficient_downloads_per_displays"] = config.getfloat("drank_coefficients", "coefficient_downloads_per_displays")
        drank_methods["drank_coefficients"]["coefficient_seens_per_displays"] = config.getfloat("drank_coefficients", "coefficient_seens_per_displays")
        drank_methods["drank_coefficients"]["coefficient_views_per_displays"] = config.getfloat("drank_coefficients", "coefficient_views_per_displays")
        drank_methods["drank_coefficients"]["coefficient_views_per_sum_downloads_and_displays"] = config.getfloat("drank_coefficients", "coefficient_views_per_sum_downloads_and_displays")
        drank_methods["drank_coefficients"]["coefficient_image_similarity"] = config.getfloat("drank_coefficients", "coefficient_image_similarity")
        drank_methods["drank_coefficients"]["coefficient_nbr_citation"] = config.getfloat("drank_coefficients", "coefficient_nbr_citation")
    
    return drank_methods


def insert_lut_intodb(drank_lut_table,drank_name,drank_description,drank_lut):
    """Save look-up dictionary scores in the LUT table. 
    Insert a new records if necessary"""
    identifier = run_sql("SELECT id from %s where drank_name=%%s" %drank_lut_table,(drank_name,))
    drank_lut=serialize_via_marshal(drank_lut)
    if identifier:          
        run_sql("UPDATE %s SET drank_name=%%s, drank_description=%%s, drank_lut=%%s WHERE id=%%s" %drank_lut_table,(drank_name,drank_description,drank_lut,identifier[0][0],))
    else:
        run_sql("INSERT INTO %s(drank_name,drank_description,drank_lut)VALUES(%%s,%%s,%%s)" %drank_lut_table,(drank_name,drank_description,drank_lut,))              
    return True

def normalize(scores, quality_scores=0,outlier=1):
    """Normalize scores in the ranking table by the method of kernel-density estimate using Gaussian kernels.
    The normalization also removes outliers."""
    
    scores= array(scores)
    scorelist = scores

    if outlier:
        outlier_above = scoreatpercentile(scorelist, 100-outlier)
        outlier_below = scoreatpercentile(scorelist, outlier)
        scorelist_new = []
        for item in scorelist:
            if float(item) >= float(outlier_below) and float(item) <= float(outlier_above):
                scorelist_new.append(float(item))
        scorelist = scorelist_new
#        print "scorelist-check:{0}".format(scorelist) 
    normalised_score  = {}
    
    #Do nothing(score all as 1) if scores are the same. A more appropriate solution is requires
    if min(scorelist)==max(scorelist):
        for item in scores:
            normalised_score[float(item)]=1
    else:            
    #   Compute the integral of a 1D pdf between two bounds. 
    #   This estimates a more robust normalized score but the range of scores must be greater than 0
        density_estimate = gaussian_kde(scorelist)
        
        if quality_scores:
            for item in scores:
                normalised_score[round(item,4)] = float(density_estimate.integrate_box_1d(-inf,item))
        else:
            for item in scores:
                normalised_score[float(item)] = float(density_estimate.integrate_box_1d(-inf,item))
    
    return normalised_score

def generate_views_lut(ranked_by,counts_table,drank_lut_table,indert_into_db=1):
    """ The function generates specific look-up table for number of page view counts
     and returns the normalized looked up scores"""    
    views_counts=run_sql("SELECT nbr_pageviews FROM %s" %counts_table)
#    views_normalised=get_percentiles_from_scores(views_counts,1)

    normalized_views=normalize(views_counts)
    
    if indert_into_db:
        success=insert_lut_intodb(drank_lut_table,ranked_by+"_views","look up table for page views",normalized_views)
        if not success:
            print "Error: could not insert views LUT into db"
            
    return normalized_views
            
def generate_downloads_lut(ranked_by,counts_table,drank_lut_table,indert_into_db=1):
    """The function generates specific look-up table for number of downloads counts 
     and returns the normalized looked up scores"""    
    downloads_counts=run_sql("SELECT nbr_downloads FROM %s" %counts_table )
#    views_normalised=get_percentiles_from_scores(views_counts,1)
    normalized_downloads=normalize(downloads_counts)
    
    if indert_into_db:
        success=insert_lut_intodb(drank_lut_table,ranked_by+"_downloads","look up table for downloads",normalized_downloads)
        if not success:
            print "Error: could not insert downloads LUT into db"
            
    return normalized_downloads          

def generate_seens_lut(ranked_by,counts_table,drank_lut_table,indert_into_db=1):
    """This function generates specific look-up table for counts of records seen by users
     and returns the normalized looked up scores """    
    seens_counts=run_sql("SELECT nbr_seens FROM %s" %counts_table )
#    views_normalised=get_percentiles_from_scores(views_counts,1)
    normalized_seens=normalize(seens_counts)
    
    if indert_into_db:
        success=insert_lut_intodb(drank_lut_table,ranked_by+"_seens","look up table for Seens",normalized_seens)
        if not success:
            print "Error: could not insert seens LUT into db"
            
    return normalized_seens       

def generate_displays_lut(ranked_by,counts_table,drank_lut_table,indert_into_db=1):
    """The function generates specific look-up table for total number of displays of all records 
    and returns the normalized looked up scores """    
    displays_counts=run_sql("SELECT nbr_displays FROM %s" %counts_table )
#    views_normalised=get_percentiles_from_scores(views_counts,1)
    normalized_displays=normalize(displays_counts)
    
    if indert_into_db:
        success=insert_lut_intodb(drank_lut_table,ranked_by+"_displays","look up table for displays",normalized_displays)
        if not success:
            print "Error: could not insert displays LUT into db"
            
    return normalized_displays 

def get_quality_records(max_nbr_of_displays,counts_table,all=1):
    """ The Function retrieves well exposed data from the pre-calculated counts 
    from rnkUSAGEDATA based on pre-configured maximum number of displays."""

    if max_nbr_of_displays:
        if all:
            exposed_records=array(run_sql("SELECT * FROM %s WHERE nbr_displays>=%%s" %counts_table, (max_nbr_of_displays,)))
        else:
            exposed_records=array(run_sql("SELECT id_bibrec,nbr_pageviews,nbr_downloads,nbr_seens,nbr_displays FROM %s WHERE nbr_displays>=%%s" %counts_table, (max_nbr_of_displays,)))
 
    return exposed_records

def generate_quality_lut(ranked_by,description,quality_scores,drank_lut_table,indert_into_db=1):
    """The function generates specific look-up table based on scores from the ranker function, 
    defined by DRANK specific configuration and inserts into the LUT table. 
    It also returns the normalized looked up scores """ 
       
#   Normalize scores first
    normalized_quality_scores=normalize(quality_scores,1,1)
    
    if indert_into_db:
        success=insert_lut_intodb(drank_lut_table,ranked_by,"look up table for " + description,normalized_quality_scores)
        if not success:
            print "Error: could not insert well exposed records LUT into db"
            
    return normalized_quality_scores


def get_freshness_score(recid):
    """  get freshness scores from bibrec table. Currently the value is estimated by 
    the number of days(the life/age of the record) since the record was added to the collection"""

    current_time = datetime.datetime.now()
    date_created=run_sql("select creation_date from bibrec where id=%s",(recid,))
    freshness_score=((current_time-date_created[0][0]).days)
    return freshness_score

def get_freshness_records(max_nbr_of_displays,counts_table):
    """ The Function retrieves less exposed data from the pre-calculated counts 
    from rnkUSAGEDATA table. based on the pre-configured max_nbr_of_displays parameter"""
    fresh_records={}
    if max_nbr_of_displays:        
        scores_record_list=run_sql("SELECT id_bibrec FROM %s WHERE nbr_displays<%%s" %counts_table,(max_nbr_of_displays,))
    
    if scores_record_list:
        for recid in scores_record_list:
            fresh_records[recid[0]]=get_freshness_score(recid[0])

    return fresh_records  

def generate_freshness_lut(ranked_by,max_nbr_of_displays,counts_table,drank_lut_table,indert_into_db=1):
    """This function generates lookup table for less exposed records, measured by the record's freshness value
    The function should be called once per each indexing cycle since its a complete disjoint set from well exposed records """
    freshness_scores=get_freshness_records(max_nbr_of_displays,counts_table)
    normalized_freshness_scores=normalize(freshness_scores.values())
    
    if indert_into_db:
        try:            
            success=insert_lut_intodb(drank_lut_table,ranked_by+"_freshness","look up table for freshness",normalized_freshness_scores)
        except:
         print "Error:Could not insert freshness into lookup table"
         pass
            
    return normalized_freshness_scores, freshness_scores


def insert_fresh_records(max_nbr_of_displays,counts_table):
    """Fresh records are separated from well exposed calculations. 
    Therefore is best to index it separately from the well-exposed records. 
    This function will be called once per each indexing cycle"""
    
    fresh_scores=get_freshness_records(max_nbr_of_displays,counts_table)
    freshness_percentile_scores, lut_fresh =get_percentiles_from_scores(fresh_scores,"both")
    
    insert_dictionary_intodb(freshness_percentile_scores, "freshness")

def get_wrd_sim_scores_from_file(wrd_sim_lut_folder): 
    """This function retrieves ALL cvs files in the wrd_sim_lut folder, 
    and reads all word similarity scores generated from user query. It also
    Attempts to remove unwanted characters([,],None) which gets read after csv-read
    
    The function requires only the wrd_sim dump folder and assumes that every csv file
    contains scores from word similarity scores  """
    
    scores_in_file=[]
    try:        
        for file in glob.glob(os.path.join(wrd_sim_lut_folder,'*.csv')):              
            scores_file=csv.reader(open(file, 'rb'))
            for line in scores_file:
                for score in line:
                    if not score =='None':
                        scores_in_file.append(float(str(score).lstrip('[').rstrip(']')))
        return scores_in_file
    except:
        print "Could not read csv files in wrd_sim dumb folder"
        pass
    
     
def generate_wrd_sim_lut(ranked_by,wrd_sim_scores,drank_lut_table,indert_into_db=1):
    """The function generates look-up table word similarity scores extracted from query-analyzer dumps
     """    
    normalized_wrdsim_scores=normalize(wrd_sim_scores,1,1)
    
    if indert_into_db:
        success=insert_lut_intodb(drank_lut_table,ranked_by+"_wrd","look up table for word Similarity",normalized_wrdsim_scores)
        if not success:
            print "Error:insert_lut_intodb"
            
    return normalized_wrdsim_scores

def get_interpolated_value(x, x1, x2, y1, y2):
    """Function that utilizes a second order newtons interpulant.    
    Get the interpolated value base on the two point (x1, y1) and (x2, y2) that x1 < x < x2
    Input:
        - x: new value
        - (x1, y1): the first point of interpolation
        - (x2, y2): the second point of interpolation
    Output:
        - interpolated value of y
    """
    y = float(y1) + ((float(y2) - float(y1))*(float(x) - float(x1)))/(float(x2)-float(x1))

    return y

def lookup_in_lut(lut,reclist):
    """This is a mapping function which looks into existing look-up table/dictionary, 
    the corresponding scores for a record list.  It also interpulates using 
    newton's second order to estimate scores not in the look up table/dictionary """
    freq_values=sorted(lut.keys()) 
    looked_up_scores={}         
    for key,score in reclist.items():
        if lut.has_key(reclist[key]):
            looked_up_scores[key] = lut[reclist[key]]

        else:
            for k in range(0, len(freq_values)-1):

                if (freq_values[k] < score) and (freq_values[k+1] > score):
                    looked_up_scores[int(key)]= get_interpolated_value(score, freq_values[k], freq_values[k+1], lut[freq_values[k]], lut[freq_values[k+1]])
                    break
#    print "record_pscore:{0}".format(record_pscore)
    return looked_up_scores

def drank_quality_function(drank_config,view,downloads,seens,displays):
    """ This is the main DRANk ranker function. At the momment, it aggregates linearly, all the quality based ranking methods.
    Note that image_similarity and nbr_citation are ignored due to lack of data for these methods
    """
    
    views_wt=drank_config["drank_coefficients"]["coefficient_nbr_views"]
    downloads_wt=drank_config["drank_coefficients"]["coefficient_nbr_downloads"]
    seens_wt=drank_config["drank_coefficients"]["coefficient_nbr_seens"]
    displays_wt=drank_config["drank_coefficients"]["coefficient_nbr_displays"]  
    views_per_displays_wt=drank_config["drank_coefficients"]["coefficient_views_per_displays"]   
    downloads_per_displays_wt=drank_config["drank_coefficients"]["coefficient_downloads_per_displays"]
    seens_per_displays_wt=drank_config["drank_coefficients"]["coefficient_seens_per_displays"]
    views_per_sum_downloads_and_displays_wt=drank_config["drank_coefficients"]["coefficient_views_per_sum_downloads_and_displays"]
#    image_similarity_wt=drank_config["drank_coefficients"]["coefficient_image_similarity"]
#    nbr_citation_wt=drank_config["drank_coefficients"]["coefficient_nbr_citation"]
        
    quality_ranker=views_wt*view+views_per_displays_wt*(view/displays)+\
                    downloads_wt*downloads+downloads_per_displays_wt*(downloads/displays)+\
                    seens_wt*seens+seens_per_displays_wt*(seens/displays)+\
                    displays_wt*displays+\
                    views_per_sum_downloads_and_displays_wt*(view/(downloads+displays))
                                                             
    return quality_ranker

def get_percentiles(values, score, method='weak'):
    """The percentile rank of a score relative to a list of scores(values). 
    the exact definition depends on the optional keyword, `method`. The function takes values, an array of scores to which a value(score)
    is compared.
    The method keyword has four variants. the first is by 'rank' which gives Average percentage ranking of score and in a case of
    multiple matches, average the percentage rankings of all matching scores. The second method and default is "weak" which
    corresponds to the definition of values by a cumulative distribution function. The third is "strict": Similar to "weak", 
    except that only values that are strictly less than the given score are counted. The last method is "mean": which gives the
    average of the "weak" and "strict" scores for the purpose of testing. """
    values = array(values)
    n = len(values)
    
    if method == 'rank':
        if not(any(values == score)):
            values = append(values, score)
            values_len = array(range(len(values)))
        else:
            values_len = array(range(len(values))) + 1.0
    
        values = sort(values)
        idx = [values == score]
        pct = (mean(values_len[idx]) / n) * 100.0
        return pct
    
    elif method == 'strict':
        return int(sum(values < score) / float(n) * 100)
    elif method == 'weak':
        return sum(values <= score) / float(n) * 100
    elif method == 'mean':
        return (sum(values < score) + sum(values <= score)) * 50 / float(n)
    else:
        raise ValueError, "method can only be 'rank', 'strict', 'weak' or 'mean'"  
    
    
def get_percentiles_from_scores(scores,kind='both'):
    """This function examines the scores being passed and calls and assigns the appropriate parameters. 
    The function may return both the percentile scores and value/percentile scores for use as lookup table, or ranking.
    Note that, I round up the percentile score allowing a normalization between 0 and 100. 
    I also choose the strict option which finds the percentile by using cumulative frequency distribution """

    percentile_scores={} 
    lut={}
    if kind=='only_lut':
        for v in scores:
            lut[int(v[0])]=get_percentiles(scores, v[0],'strict')
        return lut
    elif kind=='only_scores':
        for k,v in scores.iteritems():
            percentile_scores[k]=round(get_percentiles(scores.values(), scores[k],'strict'),2)  
        return percentile_scores
    
    elif kind=='both':    
        for k,v in scores.iteritems():
            percentile_scores[k]=round(get_percentiles(scores.values(), scores[k],'strict'),2)
            lut[v]=round(get_percentiles(scores.values(), scores[k],'strict'))
        return percentile_scores,lut
    else:
        raise ValueError, "You can only use: only_lut, only_scores or both"

def add_new_drank_method(dictname):
    """Adds a new drank method to rnkMETHOD Table
    dictname - the "name" for the rank method, to be used by drank daemon """

    try:
        run_sql("INSERT INTO rnkMETHOD (name) VALUES (%s)", (dictname,))
        identifier = run_sql("SELECT id FROM rnkMETHOD WHERE name=%s", (dictname,))
        if identifier:
            return identifier
        else:
            raise StandardError
    except StandardError, e:
        return (0, e)

def delete_from_db(dictname):
    """Delete the ranking dictionary fron ranking table."""
    
    identifier = run_sql("SELECT id from rnkMETHOD where name=\"%s\"" % dictname)
    if identifier:
        run_sql("DELETE FROM rnkMETHODDATA WHERE id_rnkMETHOD=\"%s\"" % identifier[0][0])
    return   

def insert_dictionary_intodb(dictionary, dictname):
    """Save ranking scores into the Ranking methods table."""
    identifier = run_sql("SELECT id from rnkMETHOD where name=\"%s\"" % dictname)
    if identifier:        
        delete_from_db(dictname)
    else:
        identifier=add_new_drank_method(dictname)
        
    serdata = serialize_via_marshal(dictionary)
    midstr  = str(identifier[0][0])
    run_sql("INSERT INTO rnkMETHODDATA(id_rnkMETHOD, relevance_data) VALUES (%s,%s)", (midstr, serdata,))

    return

def check_drank(max_displays,counts_table):
    """ """
    if int(run_sql("SELECT count(nbr_displays) from %s where nbr_pageviews>0"%counts_table)[0][0])<max_displays:
        good_to_run_drank=0
    elif int(run_sql("SELECT count(nbr_seens) from %s where nbr_seens>0" %counts_table)[0][0])<max_displays:
        good_to_run_drank=0
    elif int(run_sql("SELECT count(nbr_pageviews) from %s where nbr_pageviews>0" %counts_table)[0][0])<max_displays:
        good_to_run_drank=0
    elif int(run_sql("SELECT count(nbr_downloads) from %s where nbr_downloads>0" %counts_table)[0][0])<max_displays:
        good_to_run_drank=0 
    else:    
        good_to_run_drank=1    
    return good_to_run_drank

def drank_indexer():  
    """ The function accepts as input the name given to the configuration file during setup, 
    mostly, drank, and calls all offline based calculations """  
    if not import_numpy:
        write_message('The numpy package could not be imported. \
                        This package is compulsory for running the Distributed ranking daemon.')
        return
    elif not import_scipy:
        write_message('The Scipy package could not be imported. \
                        This package is compulsory for running the Distributed ranking Daemon.')
        return
    elif not good_to_run_drank:
        write_message('Limited amounts of counts in counts table. \
                        Distributed ranking daemon requires queryanalizer to run with enough counts')
    else: 

        startime = time.time()
        wrd_sim_folder = task_get_option('wrd-sim-folder','/drank_dumps/wrd_sim_lut/' )
        insert_counts = task_get_option('insert-counts', bool(True))
        insert_quality = task_get_option('insert-quality', bool(True))
        insert_wrd = task_get_option('insert-wrd', bool(True))
        insert_fresh = task_get_option('insert-fresh', bool(True))
        verbose = task_get_option('verbose', 0)
        
        task_update_progress("Initializing drank processes ...")
        write_message("Initializing drank processes ...")
        try:
            wrd_sim_lut_folder=os.path.join(CFG_LOGDIR + wrd_sim_folder)  
        except StandardError, e:
            pass     
        task_sleep_now_if_required()
        rank_methods = run_sql("SELECT name from rnkMETHOD where name not like %s"\
                               " and name not like %s and name not like %s" \
                               " and name not like %s and name not like %s", \
                               ('%_displays','%_seens','%_views','%_downloads','%_freshness'))
        for (rank_method,) in rank_methods:
    
            try:
                file = CFG_ETCDIR + "/bibrank/" + rank_method + ".cfg"
    #            file = "/home/pette/src/invenio/modules/bibrank/etc/" + rank_method + ".cfg"
                config = ConfigParser.ConfigParser()
                config.readfp(open(file))
            except StandardError, e:
                pass
            if os.path.exists(file) and config.get("rank_method", "function")== "distributed_ranking":        
                if config.has_section("distributed_ranking"): 
                    #get configuration
                    ranked_by = config.get("distributed_ranking", "ranked_by")
    
                    drank_config=get_config(ranked_by)
                    drank_lut_table=drank_config["drank_parameters"]["drank_lut_table"]
                    
                    counts_table=drank_config["drank_parameters"]["counts_table"]   
                    max_nbr_of_displays=drank_config["drank_parameters"]["max_nbr_of_displays"]        
                    description=drank_config[drank_config["function"]]["description"]    
        
                    ## Normalize all Counts(currently only four(4) counts are considered
                    write_message("Indexing Drank Method and %s " % config.get("distributed_ranking", "description"))
                    write_message("Starting Task #1/14: Generating View LUT from counts")
                    task_update_progress("Starting Task #1/14: Generating View LUT from counts")                  
                    #generate normalized lookup table for counts of pageviews per record 
                    normalized_pageviews=generate_views_lut(ranked_by,counts_table,drank_lut_table)
                    
                    write_message("Done Task #1/14: View counts generated and its LUT inserted successfully")
                    task_update_progress("Done Task #1/14: View counts generated and its LUT inserted successfully")  
                                   
                    #generate normalized lookup table for number of downloads count
                    write_message("Starting Task #2/14: Generating Download LUT from counts")
                    task_update_progress("Starting Task #2/14: Generating Download LUT from counts")
                    normalized_downloads=generate_downloads_lut(ranked_by,counts_table,drank_lut_table)
                    write_message("Done Task #2/14: Download counts generated and its LUT inserted successfully")
                    task_update_progress("Done Task #2/14: Download counts generated and its LUT inserted successfully") 
                    
                    #generate normalized lookup table for counts of displays per record
                    write_message("Starting Task #3/14: Generating Displays LUT from counts")
                    task_update_progress("Starting Task #3/14: Generating Displays LUT from counts")                
                    normalized_displays=generate_displays_lut(ranked_by,counts_table,drank_lut_table)
                    write_message("Done Task #3/14: Displays counts generated and its LUT inserted successfully")
                    task_update_progress("Done Task #3/14: Displays counts generated and its LUT inserted successfully") 
                                    
                    #generate normalized lookup table for counts of record seens
                    write_message("Starting Task #4/14: Generating Seens LUT from counts")
                    task_update_progress("Starting Task #4/14: Generating Seens LUT from counts")
                    normalized_seens=generate_seens_lut(ranked_by,counts_table,drank_lut_table)                 
                    write_message("Done Task #4/14: Seens counts generated and its LUT inserted successfully")
                    task_update_progress("Done Task #4/14: Seens counts generated and its LUT inserted successfully") 
                    
                    #Retrieve well exposed records based on pre-configured maximum number of displays to qualify a record to be well exposed
                    write_message("Starting Task #5/14: Retrieving Well exposed Records from counts table")
                    task_update_progress("Starting Task #5/14: Retrieving Well exposed Records from counts table")                
                    well_exposed_records=get_quality_records(max_nbr_of_displays,counts_table,0)             
                    write_message("Done Task #5/14: Well exposed Records retrieved successfully")
                    task_update_progress("Done Task #5/14: Well exposed Records retrieved successfully")  
                                   
                    ##Regenerate scores for records based based on generated lookup table for all counts
                    write_message("Starting Task #6/14: Looking up scores for quality counts...")
                    task_update_progress("Starting Task #6/14:  Looking up scores for quality counts...")                 
                    lookedup_pageviews=lookup_in_lut(normalized_pageviews,dict(well_exposed_records[:,[0,1]]))
                    lookedup_downloads=lookup_in_lut(normalized_downloads,dict(well_exposed_records[:,[0,2]]))
                    lookedup_seens=lookup_in_lut(normalized_seens,dict(well_exposed_records[:,[0,3]]))
                    lookedup_displays=lookup_in_lut(normalized_displays,dict(well_exposed_records[:,[0,4]]))  
                    write_message("Done Task #6/14: All scores for all counts generated")
                    task_update_progress("Done Task #6/14: All scores for all counts generated")  
                                 
                    #Calculate quality scores from counts based on current(linear) DRANK function 
                    write_message("Starting Task #7/14: Calculating Quality scores using current DRANK ranker function...")
                    task_update_progress("Starting Task #7/14:  Calculating Quality scores using current DRANK ranker function...")               
                    quality_scores= drank_quality_function(drank_config,array(lookedup_pageviews.values()),array(lookedup_downloads.values()),array(lookedup_seens.values()),array(lookedup_displays.values()))
                    write_message("Done Task #7/14: Quality scores generated")
                    task_update_progress("Done Task #7/14: Quality scores generated") 
                                 
                    #Generate normalized lookup table for quality scores
                    write_message("Starting Task #8/14: Generating normalized lookup table for quality scores")
                    task_update_progress("Starting Task #8/14:  Generating normalized lookup table for quality scores")                 
                    normalized_quality_scores=generate_quality_lut(ranked_by,description,quality_scores,drank_lut_table,1)
                    write_message("Done Task #8/14: Normalized lookup table for quality scores generated and its LUT inserted successfully")
                    task_update_progress("Done Task #8/14: normalized lookup table for quality scores generated and its LUT inserted successfully ") 
                    
                    #Attach/stack generated scores to the respective records
                    write_message("Starting Task #9/14: Stacking quality scores to respective records")
                    task_update_progress("Starting Task #9/14: Stacking quality scores to respective records")                
                    quality_scores_generated=dict(column_stack([well_exposed_records[:,[0]],quality_scores]))
                    write_message("Done Task #9/14: Quality dictionary scores generated")
                    task_update_progress("Done Task #9/14: Quality dictionary scores generated")
                    
                    #Generate look up table for generated quality scores
                    write_message("Starting Task #10/14: Calculating normalized percentile scores for generated quality scores")
                    task_update_progress("Starting Task #10/14:  Calculating normalized percentile scores for generated quality scores")
                    lookedup_quality_scores=lookup_in_lut(normalized_quality_scores,quality_scores_generated)                
                    write_message("Done Task #10/14: Percentile scores successfully looked up")
                    task_update_progress("Done Task #10/14: Percentile scores successfully looked up")              
                    #Calculate the percentile scores for each quality score
                    if insert_counts:
                        write_message("Starting Task #11/14: Inserting percentile scores for all counts")
                        task_update_progress("Starting Task #11/14:  Inserting percentile scores for all counts") 
                        
                        #insert views scores dictionary                 
                        views_scores_dictionary=get_percentiles_from_scores(lookedup_pageviews,"only_scores")
                        insert_dictionary_intodb(views_scores_dictionary,ranked_by+"_views")
                                           
                        #insert views scores dictionary  
                        downloads_scores_dictionary=get_percentiles_from_scores(lookedup_downloads,"only_scores")
                        insert_dictionary_intodb(downloads_scores_dictionary,ranked_by+"_downloads")
                        
                        #insert views scores dictionary  
                        seens_scores_dictionary=get_percentiles_from_scores(lookedup_seens,"only_scores")
                        insert_dictionary_intodb(seens_scores_dictionary,ranked_by+"_seens")
                        
                        #insert views scores dictionary  
                        displays_scores_dictionary=get_percentiles_from_scores(lookedup_displays,"only_scores")
                        insert_dictionary_intodb(displays_scores_dictionary,ranked_by+"_displays")
                        
                        write_message("Done Task #11/14: Percentile scores for all counts successfully inserted into database")
                        task_update_progress("Done Task #11/14: Percentile scores for all counts successfully inserted into database")
                    
                    if insert_quality:
                        #generate percentile scores for aggregated quality records
                        write_message("Starting Task #12/14: Calculating percentile scores for aggregated quality records")
                        task_update_progress("Starting Task #12/14:  Calculating percentile scores for aggregated quality records")
                        quality_scores_dictionary=get_percentiles_from_scores(lookedup_quality_scores,"only_scores")
        
                        insert_dictionary_intodb(quality_scores_dictionary, ranked_by)
                        write_message("Done Task #12/14: Percentile scores for aggregated quality records successfully inserted into database")
                        task_update_progress("Done Task #12/14: Percentile scores for aggregated quality records successfully inserted into database")
                    
                    #Read scores from file generated by DRANK Logger(queryanalyzer) and generate LUT from scores
                    if insert_wrd:
                        write_message("Starting Task #13/14: Generating normalized lookup table for word similarity scores")
                        task_update_progress("Starting Task #13/14:  Generating normalized lookup table for word similarity scores")
                        #get scores from file
                        wrd_sim_scores=get_wrd_sim_scores_from_file(wrd_sim_lut_folder)
                        
                        #normalize and insert LUT for word similarity
                        normalized_wrdsim_scores=generate_wrd_sim_lut(ranked_by,wrd_sim_scores,drank_lut_table,1)
                        write_message("Done Task #13/14: Lookup table for word similarity successfully inserted into database")
                        task_update_progress("Done Task #13/14: Lookup table for word similarity successfully inserted into database")
                        
                        #generate freshness scores(age of record in days) for less exposed records. 
                        #Less exposed records cannot be ranked by the ranker function which is based on displays
                    if insert_fresh:
                        write_message("Starting Task #14/14: Generating normalized lookup table for less exposed records")
                        task_update_progress("Starting Task #14/14:  Generating normalized lookup table for less exposed records")
                        normalized_freshness_scores, freshness_scores=generate_freshness_lut(ranked_by,max_nbr_of_displays,counts_table,drank_lut_table,indert_into_db=1)                
                        lookedup_freshness_scores=lookup_in_lut(normalized_freshness_scores, freshness_scores)
                        freshness_scores_dictionary=get_percentiles_from_scores(lookedup_freshness_scores,"only_scores")                
                        insert_dictionary_intodb(freshness_scores_dictionary,ranked_by+"_freshness")
                        write_message("Done Task #14/14: Less exposed records successfully scored with freshness values")
                        task_update_progress("Done Task #14/14:  Less exposed records successfully scored with freshness values")
                write_message("All task completed for DRANK METHOD:%s" %description)
        
        time_taken=round(time.time()-startime,1)
        task_update_progress("All Done in %s seconds" %time_taken)
        task_sleep_now_if_required(can_stop_too=True)
        write_message("All Done in %s seconds" %time_taken)
        return True 

def main():
    """main function for constructing bibtask for DRANK."""

    task_init(authorization_action='run_drank_indexer',
              authorization_msg="Running distributed ranking indexer",
              description="""Description:
    Distributed Ranking Indexer(dranker): This task/daemon pre-generates normalized scores 
    and look up tables for all counts recorded by DRANK Logger. It also pre-calculates the quality scores
    based on the current ranking function and populates a ranking dictionary for these scores.\n""",
    help_specific_usage="""
    -d, --wrd-score-folder=FOLDER location of scores dump folder
    -c, --insert-counts=COUNTS    insert percentile scores for all counts,
                                  insert(True), don't insert(False)[default=True]
    -q, --insert-quality=QUALITY  insert percentile scores for well exposed records
                                  using ranker function
                                  insert(True), don't insert(False)[default=True]
    -w, --insert-wrd=WRD          insert lookup table for word similarity using files from logger
                                  insert(True), don't insert(False)[default=True]
    -f, --insert-fresh=FRESHNESS  insert percentile scores(freshness) for less exposed records
                                  insert(True), don't insert(False)[default=True]
    
    Examples:
    Run periodically (eg. Daily) and insert dictionaries into Ranking Table:
    $ dranker -f1  -u admin -s 1d \n""",
              specific_params=("d:c:q:w:f:",[
                               "wrd-score-folder=",
                               "insert-counts=",
                               "insert-quality="
                               "insert-wrd=",
                               "insert-fresh="]),                            
            task_submit_elaborate_specific_parameter_fnc = get_elaborate_submit_param,
            task_run_fnc = drank_indexer)

if __name__ == '__main__':
    main()