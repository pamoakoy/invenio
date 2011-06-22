# -*- coding: utf-8 -*-
##
## This file is part of CDS Invenio.
## Copyright (C) 2002, 2003, 2004, 2005, 2006, 2007, 2008 CERN.
##
## CDS Invenio is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 2 of the
## License, or (at your option) any later version.
##
## CDS Invenio is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with CDS Invenio; if not, write to the Free Software Foundation, Inc.,
## 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA

__revision__ = "$Id$"

""" Calculation of rnkUSAGEDATA and
    rnkDOWNLOADS/rnkEXTLINKS/rnkPAGEVIEWS ranks for rnkPAGEVIEWS and rnkDOWNLOADS """
#import cPickle
#import zlib
from invenio.config import \
     CFG_SITE_LANG, CFG_SITE_NAME,\
     CFG_LOGDIR, CFG_PATH_MYSQL, CFG_DATABASE_HOST, \
     CFG_DATABASE_USER, CFG_DATABASE_PASS, CFG_DATABASE_NAME, \
     CFG_SITE_URL, CFG_SITE_SECURE_URL, CFG_CACHEDIR, \
     CFG_BIBRANK_DRANK_POSTPROCESSING_TIME_WINDOW, \
     CFG_ETCDIR, \
     CFG_PREFIX
import datetime
import os, sys, cgi,time,ConfigParser,string
from invenio.errorlib import register_exception
from invenio.dbquery import run_sql, DatabaseError, serialize_via_marshal, deserialize_via_marshal
#from invenio.bibrank_record_sorter import rank_records
from invenio.intbitset import intbitset
from invenio.shellutils import run_shell_command, escape_shell_arg              
from invenio.bibtask import task_init, write_message, task_set_option, \
        task_get_option, task_update_progress, task_update_status, \
        task_get_task_param, task_sleep_now_if_required
import numpy as np

#from decimal import *
#from invenio.search_engine import get_fieldvalues as get_format
#from invenio.search_engine import perform_request_search as search
# Update of the seens will happen from the "time_limit" time up to the current time
#last_run_time_file = os.path.join(CFG_CACHEDIR, 'queryanalyzer_runtime')
#from invenio.bibrank_record_sorter import rank_records as ws

def get_current_time():
    current_time = datetime.datetime.now()
    return current_time

def itemgetter(*items):
    if len(items) == 1:
        item = items[0]
        def g(obj):
            return obj[item]
    else:
        def g(obj):
            return tuple(obj[item] for item in items)
    return g  

def get_rnkMethod_id(table,name):
    """Returns the name of the rank method, either in default language or given language.
    rank_method_function = short name or function name defined in config file
    ln - the language to get the name in
    type - which name "type" to get."""
    #TODO: complete it by adding table 
    try:
        rmid = run_sql("SELECT id FROM %s where name='%s'" % rank_method_code)
        if rmid_id:
            id = id[0][0]
            return id
        else:
            raise Exception
    except Exception, e:
        write_message("Cannot run function specified, either the given Ranking function is wrongly spelt or has not been added to DRank.")
        raise Exception
    
def get_config(rank_method_code,wrd="wrd"):
    """ Get all configuration from drank config file"""

#    global methods
#    bibrank_meths = run_sql("SELECT name from rnkMETHOD")
    drank_methods = {}
    try:
        file = CFG_ETCDIR + "/bibrank/" + rank_method_code + ".cfg"
        config = ConfigParser.ConfigParser()
        config.readfp(open(file))
    except StandardError, e:
        pass
    #        cfg_parameters=config.get("parameters","")
    drank_bibrank_function = config.get("rank_method", "function")
    drank_methods["function"] = drank_bibrank_function
    if config.has_section(drank_bibrank_function): 
        drank_methods[drank_bibrank_function] = {}  
        if config.has_option(drank_bibrank_function, "table"):
            drank_methods[drank_bibrank_function]["table"]= config.get(drank_bibrank_function,"table")  
        drank_methods[drank_bibrank_function]["prefix"] = config.get(drank_bibrank_function, "relevance_number_output_prologue")
        drank_methods[drank_bibrank_function]["postfix"] = config.get(drank_bibrank_function, "relevance_number_output_epilogue")
        
    if config.has_section("drank_parameters"): 
        drank_methods["drank_parameters"] = {}           
        drank_methods["drank_parameters"]["drank_dictionary_table"] = config.get("drank_parameters", "drank_dictionary_table")
        drank_methods["drank_parameters"]["counts_table"]=config.get("drank_parameters", "counts_table")         
        drank_methods["drank_parameters"]["default_method"] = config.get("drank_parameters", "default_method") 
        drank_methods["drank_parameters"]["max_nbr_of_displays"] = config.getfloat("drank_parameters", "max_nbr_of_displays")        
        drank_methods["drank_parameters"]["weight_of_relevance_exposed"] = config.getfloat("drank_parameters", "weight_of_relevance_exposed")
        drank_methods["drank_parameters"]["weight_of_relevance_fresh"] = config.getfloat("drank_parameters", "weight_of_relevance_fresh")
        drank_methods["drank_parameters"]["weight_of_exposed"] = config.getfloat("drank_parameters", "weight_of_exposed")
        drank_methods["drank_parameters"]["weight_of_fresh"] = config.getfloat("drank_parameters", "weight_of_fresh")
        drank_methods["drank_parameters"]["percentage_of_freshness"] = config.getfloat("drank_parameters", "percentage_of_freshness")
        drank_methods["drank_parameters"]["min_freshness_value"] = config.getfloat("drank_parameters", "min_freshness_value")      
     
    if config.has_section("relevance_parameters"):
        drank_methods["relevance_parameters"]={}
        drank_methods["relevance_parameters"]["max_word_occurence"] = float(config.get("relevance_parameters", "max_word_occurence"))
        drank_methods["relevance_parameters"]["min_word_occurence"] = float(config.get("relevance_parameters", "min_word_occurence"))
        drank_methods["relevance_parameters"]["min_word_length"] = int(config.get("relevance_parameters", "min_word_length"))
        drank_methods["relevance_parameters"]["min_nr_words_docs"] = int(config.get("relevance_parameters", "min_nr_words_docs"))
        drank_methods["relevance_parameters"]["max_nr_words_upper"] = int(config.get("relevance_parameters", "max_nr_words_upper"))
        drank_methods["relevance_parameters"]["max_nr_words_lower"] = int(config.get("relevance_parameters", "max_nr_words_lower"))
        drank_methods["relevance_parameters"]["default_min_relevance"] = int(config.get("relevance_parameters", "default_min_relevance"))
        
    if config.has_section("drank_methods"):
        i = 1
        drank_methods["drank_methods"] ={}
        while config.has_option("drank_methods", "drank%s" % i):
            mi=config.get("drank_methods", "drank%s" % i)
            drank_methods["drank_methods"][i]=mi
            i += 1       
    
        for k,drank_method in drank_methods["drank_methods"].iteritems():    
            if config.has_section(drank_method):
                drank_method_section=drank_method
                write_message("Using custom configuration")
            else:
                write_message("Using default configuration")
                if config.has_section("default_config"):
                    drank_method_section="default_config" 
                else:
                    raise Exception("Error in configuration file: %s" % (CFG_ETCDIR + "/bibrank/" + rank_method_code + ".cfg"))  
            
            drank_methods[drank_method] = {}       
            drank_methods[drank_method]["description"] = config.get(drank_method_section, "description")   
            drank_methods[drank_method]["coefficient_nbr_downloads"] = config.getfloat(drank_method_section, "coefficient_nbr_downloads")
            drank_methods[drank_method]["coefficient_nbr_displays"] = config.getfloat(drank_method_section, "coefficient_nbr_displays")
            drank_methods[drank_method]["coefficient_nbr_seens"] = config.getfloat(drank_method_section, "coefficient_nbr_seens")
            drank_methods[drank_method]["coefficient_nbr_views"] = config.getfloat(drank_method_section, "coefficient_nbr_views")
            drank_methods[drank_method]["coefficient_downloads_per_displays"] = config.getfloat(drank_method_section, "coefficient_downloads_per_displays")
            drank_methods[drank_method]["coefficient_seens_per_displays"] = config.getfloat(drank_method_section, "coefficient_seens_per_displays")
            drank_methods[drank_method]["coefficient_views_per_displays"] = config.getfloat(drank_method_section, "coefficient_views_per_displays")
            drank_methods[drank_method]["coefficient_views_per_sum_downloads_and_displays"] = config.getfloat(drank_method_section, "coefficient_views_per_sum_downloads_and_displays")
            drank_methods[drank_method]["coefficient_image_similarity"] = config.getfloat(drank_method_section, "coefficient_image_similarity")
            drank_methods[drank_method]["coefficient_nbr_citation"] = config.getfloat(drank_method_section, "coefficient_nbr_citation")
    

    try:
        file = CFG_ETCDIR + "/bibrank/" + wrd + ".cfg"
        wrd_config = ConfigParser.ConfigParser()
        wrd_config.readfp(open(file))
    except Error, e:
        wrd_config=False
        
    if wrd_config:
        word_similarity = wrd_config.get("rank_method", "function")
        if wrd_config.has_section(word_similarity):
            drank_methods["wrd"] = {}
            drank_methods["wrd"]["function"] = word_similarity
            drank_methods["wrd"]["prefix"] = wrd_config.get(word_similarity, "relevance_number_output_prologue")
            drank_methods["wrd"]["postfix"] = wrd_config.get(word_similarity, "relevance_number_output_epilogue")
            drank_methods["wrd"]["chars_alphanumericseparators"] = r"[1234567890\!\"\#\$\%\&\'\(\)\*\+\,\-\.\/\:\;\<\=\>\?\@\[\\\]\^\_\`\{\|\}\~]"
        else:
            raise Exception("Error in configuration file: %s" % (CFG_ETCDIR + "/bibrank/" + wrd + ".cfg"))
    
        i8n_names = run_sql("""SELECT ln,value from rnkMETHODNAME,rnkMETHOD where id_rnkMETHOD=rnkMETHOD.id and rnkMETHOD.name=%s""", (wrd,))
        for (ln, value) in i8n_names:
            drank_methods["wrd"][ln] = value
    
        if wrd_config.has_option(word_similarity, "table"):
            drank_methods["wrd"]["rnkWORD_table"] = wrd_config.get(word_similarity, "table")
            drank_methods["wrd"]["col_size"] = run_sql("SELECT count(*) FROM %sR" % drank_methods[wrd]["rnkWORD_table"][:-1])[0][0]
    
        if wrd_config.has_option(word_similarity, "stemming") and wrd_config.get(word_similarity, "stemming"):
            try:
                drank_methods["wrd"]["stemmer"] = wrd_config.get(word_similarity, "stemming")
            except Exception,e:
                pass
    
        if wrd_config.has_option(word_similarity, "stopword"):
            drank_methods["wrd"]["stopwords"] = wrd_config.get(word_similarity, "stopword")
    return drank_methods


def insert_unique_records(bibrec_id,current_time):
    """ Function to populate rnkDRANKADATA table with unique data per record from rnkUSAGEDATA
    This is created temporally to address issues with rnkUSAGEDATA dumps, which currently provides redundant data
    NB: The way this calculation is done is done so as to get more data for testing """
    #TODO: get the correct counts and edit the meaning of the above function
    #
    #Calculate from rnkUSAGEDATA, the respective maximum or recent values
    res=run_sql("SELECT id_bibrec,nbr_displays,nbr_downloads,nbr_seens,nbr_pageviews,time_stamp FROM rnkUSAGEDATA where id_bibrec=%s",(bibrec_id,))
    max_nbr_displays=int(sum([row[1] for row in res]))
    max_nbr_downloads=int(max([row[2] for row in res]))
    max_nbr_seens=int(max([row[3] for row in res]))
    max_nbr_pageviews=int(max([row[4] for row in res]))
    max_time_stamp=max([row[5] for row in res])
    date_created=[row[0] for row in run_sql("select creation_date from bibrec where id=%s",(bibrec_id,))]
    freshness=((current_time-date_created[0]).days)
    max_nbr_pageviews_decayed=int(0)
    max_nbr_downloads_decayed=int(0)
    max_nbr_seens_decayed=int(0)
    max_nbr_displays_decayed=int(0)
    #
    # insert records uniquely into the table "rnkDRANKDATA"
    insert_str="""INSERT INTO rnkDRANKDATA(id_bibrec,
    nbr_pageviews,
    nbr_downloads,
    nbr_seens,
    nbr_displays,
    freshness,
    nbr_pageviews_decayed,
    nbr_downloads_decayed,
    nbr_displays_decayed,
    nbr_seens_decayed,
    time_stamp) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s) ON DUPLICATE KEY UPDATE
     nbr_pageviews=VALUES(nbr_pageviews),nbr_downloads=VALUES(nbr_downloads),
    nbr_seens=VALUES(nbr_seens),nbr_displays=VALUES(nbr_displays),freshness=VALUES(freshness),time_stamp=VALUES(time_stamp)"""
    arg_values=bibrec_id,max_nbr_pageviews,max_nbr_downloads,\
           max_nbr_seens,max_nbr_displays,freshness,max_nbr_pageviews_decayed,\
           max_nbr_downloads_decayed,max_nbr_seens_decayed,max_nbr_displays_decayed,\
           current_time,

#    run_sql("DELETE FROM rnkDRANKDATA WHERE id_bibrec=%s", (bibrec_id, ))
    run_sql(insert_str,(arg_values))
    
    nbr_rec_updates=[row[0] for row in run_sql("select count(id_bibrec) from rnkDRANKDATA")]
    #
    return nbr_rec_updates

def repopulate_data_from_queryanalyzer():
    """ The function prepares the records and then inserts by calling the correct function """      
    recid_list=[row[0] for row in run_sql("select distinct id_bibrec from rnkUSAGEDATA")]
   
    for recid in recid_list:
        current_time=get_current_time()
        nbr_rec_updates=insert_unique_records(recid,current_time)
        
    print "Done:  updated " +  str(nbr_rec_updates) + "Records"
    
    
def get_quality_records(max_nbr_of_displays):
    """ The Function retrieves fresh and exposed data from the pre-calculated counts 
    from rnkDRANKDATA and not rnkUSAGEDATA.
    TODO: Merge both tables to have a single data table"""
    
    if max_nbr_of_displays:
        exposed_records=np.array(run_sql("SELECT * FROM rnkDRANKDATA WHERE nbr_displays>=%s",(max_nbr_of_displays,)))
        scores_fresh_records=dict(run_sql("SELECT id_bibrec, freshness FROM rnkDRANKDATA WHERE nbr_displays<%s",(max_nbr_of_displays,)))
    print "quality records retrieved"  
#     print "exposed_records:{0}".format(exposed_records)   
    return exposed_records, scores_fresh_records    
    
def get_score_from_exposed_record(drank_config,drank_method,nbr_of_seens, nbr_of_displays,\
                                   nbr_of_downloads,nbr_of_pageviews,nbr_of_citations):
    """The function calculates the exposed scores based on predefined formula and the various
    coefficients defined in the configuration file . Note, you will need to change drank formula physically here.""" 
                                   
    coefficient_nbr_seens                   =drank_config[drank_method]["coefficient_nbr_seens"]
    coefficient_nbr_views                   =drank_config[drank_method]["coefficient_nbr_views"]
    coefficient_nbr_citation                =drank_config[drank_method]["coefficient_nbr_citation"]
    coefficient_nbr_displays                =drank_config[drank_method]["coefficient_nbr_displays"]
    coefficient_nbr_downloads               =drank_config[drank_method]["coefficient_nbr_downloads"]
    coefficient_image_similarity            =drank_config[drank_method]["coefficient_image_similarity"]
    coefficient_seens_per_displays          =drank_config[drank_method]["coefficient_seens_per_displays"]
    coefficient_views_per_displays          =drank_config[drank_method]["coefficient_views_per_displays"]
    coefficient_downloads_per_displays      =drank_config[drank_method]["coefficient_downloads_per_displays"]
    coefficient_views_per_sum_downloads_and_displays=drank_config[drank_method]["coefficient_views_per_sum_downloads_and_displays"]                                                               
#    print "this is:{0}".format(coefficient_views_per_sum_downloads_and_displays)
    #Method 1    
    raw_score=coefficient_nbr_seens*(nbr_of_seens)+\
    coefficient_nbr_views*(nbr_of_pageviews)+\
    coefficient_nbr_citation*(nbr_of_citations)+\
    coefficient_nbr_displays*(nbr_of_displays)+\
    coefficient_nbr_downloads*(nbr_of_downloads)+\
    coefficient_seens_per_displays*(np.divide(nbr_of_seens,nbr_of_displays))+\
    coefficient_views_per_displays*(np.divide(nbr_of_pageviews,nbr_of_displays))+\
    coefficient_downloads_per_displays*(np.divide(nbr_of_downloads,nbr_of_displays))+\
    coefficient_views_per_sum_downloads_and_displays*(np.divide(nbr_of_pageviews,(nbr_of_downloads+nbr_of_displays)))    
    
#        coefficient_image_similarity*()+\    
    return raw_score
def look_up(res,reclist):
    """ """
    if res:
        lut= deserialize_via_marshal(res[0][0])  
        freq_values=sorted(lut.keys()) 
#    print "reclist:{0}".format(reclist)
    record_pscore={}         
    for key,score in reclist.iteritems():
        print score
        s=score
        if lut.has_key(reclist[key]):
            record_pscore[s] = lut[reclist[key]]
            print "reclist:{0}".format(record_pscore)
        else:
            for k in range(0, len(freq_values)-1):
                print "using this"
                if (freq_values[k] < score) and (freq_values[k+1] > score):
                    record_pscore[k]= get_interpolated_value(score, freq_values[k], freq_values[k+1], lut[freq_values[k]], lut[freq_values[k+1]])
                    break
    print "record_pscore:{0}".format(record_pscore)
    return record_pscore
    
def catch_lookup_table(drank_dictionary_table,exposed_records):
    """ """
    
    views_res=run_sql("""SELECT dictionary FROM %s WHERE name=%%s and dictionary_type=%%s""" %drank_dictionary_table,("views","views",))
    views=look_up(views_res,dict(exposed_records[:,[0,1]]))
    downloads_res=run_sql("""SELECT dictionary FROM %s WHERE name=%%s and dictionary_type=%%s""" %drank_dictionary_table,("downloads","downloads",))
    downloads=look_up(downloads_res,dict(exposed_records[:,[0,2]]))    
    seens_res=run_sql("""SELECT dictionary FROM %s WHERE name=%%s and dictionary_type=%%s""" %drank_dictionary_table,("seens","seens",))
    seens=look_up(seens_res,dict(exposed_records[:,[0,3]]))
    displays_res=run_sql("""SELECT dictionary FROM %s WHERE name=%%s and dictionary_type=%%s""" %drank_dictionary_table,("displays","displays",))
    displays=look_up(displays_res,dict(exposed_records[:,[0,4]]))
    
#    res=run_sql("SELECT dictionary FROM %s where dictionary_type='%s' and name='%s'"% (count_table,dict_type,name,))
#    views=
#    res=run_sql("SELECT dictionary FROM rnkDRANKLUT where dictionary_type='%s' and name='%%s'" % 'wrd_lut',('word similarity',))
    return views,downloads,seens,displays
        
def generate_new_scores(exposed_records,drank_method,drank_dictionary_table,drank_config):
    """ The function prepares the parameters from the config file and calls the correct function"""  

    #print orec
#    res=exposed_records[0]
#    views=res[:,[1,2]]
    
    views,downloads,seens,displays=catch_lookup_table(drank_dictionary_table,exposed_records)
#    pscore=lookup_percentile_scores(relevance_wrdsim,"wrd_lut","wrd") 
#    pscore=lookup_percentile_scores(relevance_wrdsim,"wrd_lut","wrd") 

    nbr_of_citations=0.0
    
    exposed_scores={}    
    for row in exposed_records:  
        print "rows:{0},{1},{2},{3},{4}".format(row[0],row[1],row[2],row[3],row[4])
        print row[1]
#        pscore=lookup_percentile_scores(relevance_wrdsim,"wrd_lut","wrd") 
        print views
        recid=row[0]        
        nbr_of_pageviews=views[row[1]]
        nbr_of_downloads=downloads[row[2]]
        nbr_of_seens=seens[row[3]]
        nbr_of_displays=displays[row[4]]
        print "norms:{0},{1},{2},{3}".format(nbr_of_pageviews,nbr_of_downloads,nbr_of_seens,nbr_of_displays)
        raw_score=\
        get_score_from_exposed_record(drank_config,drank_method,nbr_of_seens, nbr_of_displays,\
                                   nbr_of_downloads,nbr_of_pageviews,nbr_of_citations)
        
        exposed_scores[recid]=float(raw_score)  
    return exposed_scores


    
 
def get_percentiles(values, score, method='weak'):
    """The percentile rank of a score relative to a list of scores(values). 
    the exact definition depends on the optional keyword, `method`. The function takes values, an array of scores to which a value(score)
    is compared.
    The method keyword has four variants. the first is by 'rank' which gives Average percentage ranking of score and in a case of
    multiple matches, average the percentage rankings of all matching scores. The second method and default is "weak" which
    corresponds to the definition of values by a cumulative distribution function. The third is "strict": Similar to "weak", 
    except that only values that are strictly less than the given score are counted. The last method is "mean": which gives the
    average of the "weak" and "strict" scores for the purpose of testing. """
    values = np.array(values)
    n = len(values)
    
    if method == 'rank':
        if not(np.any(values == score)):
            values = np.append(values, score)
            values_len = np.array(range(len(values)))
        else:
            values_len = np.array(range(len(values))) + 1.0
    
        values = np.sort(values)
        idx = [values == score]
        pct = (np.mean(values_len[idx]) / n) * 100.0
        return pct
    
    elif method == 'strict':
#        print "here:{0}".format(int(sum(values < score) / float(n) * 100))
        return int(sum(values < score) / float(n) * 100)
    elif method == 'weak':
        return sum(values <= score) / float(n) * 100
    elif method == 'mean':
        return (sum(values < score) + sum(values <= score)) * 50 / float(n)
    else:
        raise ValueError, "method can only be 'rank', 'strict', 'weak' or 'mean'"  
    
    
def get_percentiles_from_scores(scores,only_lut=0):
    """Prepare data and call the correct method. 
    The function will return both the percentile scores and value/percentile scores for use as lookup table.
    Note that, I round up the percentile score allowing a normalization between 0 and 100. 
    I also choose the strict option which finds the percentile by using cfd """
    print 'scores:{0}'.format(scores)
    percentile_scores={} 
    lut={}
    if only_lut:
        for v in scores:
#            print v[0]
#        percentile_scores[k]=round(get_percentiles(pageviews, scores[k],'strict'),0)
            lut[v[0]]=get_percentiles(scores, v[0],'strict')
        print "luts:{0}".format(lut)
        return lut
    else:    
        for k,v in scores.iteritems():
            percentile_scores[k]=round(get_percentiles(scores.values(), scores[k],'strict'),0)
            lut[v]=round(get_percentiles(scores.values(), scores[k],'strict'))
        return percentile_scores,lut

def get_doc_format(recid):
    """This function is not currently being used, but its there to be called 
    to get the format of a records when we consider format based ranking """
    doc_format = get_format(k, '980__a')
    return doc_format
                
def drank_aggregator(run):
    """This function calls set of functions Run the indexing task.  The row argument is the BibSched task
    queue row, containing if, arguments, etc.
    Return 1 in case of success and 0 in case of failure.
    """
    #TODO read from config file rank_method_code,alpha,max_nbr_of_displays
    new_dict=False
    alpha=0.5
    max_nbr_of_displays=2
    rank_method_code='drank'
    id=get_rnkMethod_id('drank')
    populate_data()
    scores_exposed_records,score_fresh_records=generate_new_scores(alpha,max_nbr_of_displays)
    exposed_percentile_scores=get_percentiles_from_scores(scores_exposed_records)
    freshness_percentile_scores=get_percentiles_from_scores(score_fresh_records)
    
    merge_2dictionaries(freshness_percentile_scores,exposed_percentile_scores)
    insert_dict_intodb(dictionary,rank_method_code,new_dict)
    ## import optional modules:


def merge_2dictionaries(dictionary_1, dictionary_2):
    """
    Input:
        - Dictionary 1
        - Dictionary 2
    Output:
        Merging of two dictionaries -> dictionary_2
    """
    for k in dictionary_1.keys():
        if dictionary_2.has_key(k): continue
#            dictionary_2[k]=int(round((dictionary_2[k]*(dictionary_2[k]+dictionary_1[k]))/100,0))
        else: dictionary_2[k] = dictionary_1[k]
#    print dictionary_2
    return dictionary_2

def merge_relevant_quality(relevance,quality):
    """ """
    quality_per_relevance={}   
    for k in relevance.keys():
        if quality.has_key(k):
            quality_per_relevance[k]=round(quality[k]*relevance[k]/100,0)
#            dictionary_2[k]=int(round((dictionary_2[k]*(dictionary_2[k]+dictionary_1[k]))/100,0))
        else: quality_per_relevance[k] = relevance[k]
    
    return quality_per_relevance   

def merge_exposed_fresh_dictionary(freshness_percentile_scores,exposed_percentile_scores):
    """ """
    if len(freshness_percentile_scores)>len(exposed_percentile_scores):
        dictionary=merge_2dictionaries(freshness_percentile_scores,exposed_percentile_scores)
#    print 'using this'
    else: 
        dictionary=merge_2dictionaries(exposed_percentile_scores,freshness_percentile_scores)
    return dictionary 

def dynamic_ranking(ranktype,run):
    """Call correct method"""
    if(ranktype==True):
        return drank_aggregator(run)
    else:
        return drank_sorter(run)  


       
def _interpolate(a, b, fraction):
    """Returns the point at the given fraction between a and b, where
    'fraction' must be between 0 and 1.
    """
    return a + (b - a)*fraction;

def get_interpolated_value(x, x1, x2, y1, y2):
    """
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

def lookup_percentile_scores(reclist,dict_type,name):
    """Look up information and if it does not exist interpolate using 
    second order newtons interpolation polynomial """
    
    res=run_sql("SELECT dictionary FROM rnkDRANKLUT where dictionary_type='%s' and name='%s'"% (dict_type,name,))
    
#    res=run_sql("SELECT dictionary FROM rnkDRANKLUT where dictionary_type='%s' and name='%%s'" % 'wrd_lut',('word similarity',))
    if res:
        lut= deserialize_via_marshal(res[0][0])  
        freq_values=sorted(lut.keys()) 

    record_pscore={}         
    for key,score in reclist.iteritems():

        if lut.has_key(reclist[key]):
            record_pscore[key] = lut[reclist[key]]

        else:
            for k in range(0, len(freq_values)-1):

                if (freq_values[k] < score) and (freq_values[k+1] > score):
                    record_pscore[k]= get_interpolated_value(score, freq_values[k], freq_values[k+1], lut[freq_values[k]], lut[freq_values[k+1]])
                    break
    return record_pscore

def look_up_scores(a, per, limit=()):
    """Calculate the score at the given 'per' percentile of the
    sequence a. If the desired quantile lies between two data points, we
    interpolate between them. If the parameter 'limit' is provided, it should be a tuple (lower,
    upper) of two values.  Values of 'a' outside this (closed)
    interval will be ignored.
    """    
    values = np.sort(a,axis=0)
    if limit:
        values = values[(limit[0] <= values) & (values <= limit[1])]

    idx = per /100. * (values.shape[0] - 1)
    if (idx % 1 == 0):
        return values[idx]
    else:
        return _interpolate(values[int(idx)], values[int(idx) + 1], idx % 1)

def get_percentile_from_lut(scores, score):
    """ """
    percentile_score=look_up_scores(scores, score)
    return percentile_score

def get_nb_zero_values(result):
    """
    Get number of zero-values based on rnkMETHODDATA and bibrec
    Input:
        - rank_method_code
    Output:
        - number of zero-values
    """

    nb_zero_values = 0

#    id_rnkMETHOD = get_id_rnkMETHOD(rank_method_code)

#    result = run_sql("SELECT relevance_data FROM rnkMETHODDATA WHERE id_rnkMETHOD = %s", (id_rnkMETHOD,))

    data_dict = {}
    if result and result[0] and result[0][0]:
        data_dict = deserialize_via_marshal(result[0][0])

        for (key, value) in data_dict.iteritems():
            # do not accept the zero-value in bibrec_dict
            if value == 0:
                nb_zero_values = nb_zero_values + 1
    print "#zero:{0}".format(nb_zero_values)
    result = run_sql("SELECT id FROM bibrec")
    for i in range(0, len(result)):
        id = result[i][0]
        if not data_dict.has_key(id):
            nb_zero_values = nb_zero_values + 1
    print "#zero after bibrec:{0}".format(nb_zero_values)
    return nb_zero_values

def histogram(a, numbins=10, defaultlimits=None, printextras=True):
    # fixme: use numpy.histogram() to implement
    """
Returns (i) an array of histogram bin counts, (ii) the smallest value
of the histogram binning, and (iii) the bin width (the last 2 are not
necessarily integers).  Default number of bins is 10.  Defaultlimits
can be None (the routine picks bins spanning all the numbers in the
a) or a 2-sequence (lowerlimit, upperlimit).  Returns all of the
following: array of bin values, lowerreallimit, binsize, extrapoints.

Returns: (array of bin counts, bin-minimum, min-width, #-points-outside-range)
"""
    a = np.ravel(a)               # flatten any >1D arrays
    if (defaultlimits is not None):
        lowerreallimit = defaultlimits[0]
        upperreallimit = defaultlimits[1]
        binsize = (upperreallimit-lowerreallimit) / float(numbins)
    else:
        Min = a.min()
        Max = a.max()
        estbinwidth = float(Max - Min)/float(numbins - 1)
        binsize = (Max-Min+estbinwidth)/float(numbins)
        lowerreallimit = Min - binsize/2.0  #lower real limit,1st bin
    bins = np.zeros(numbins)
    extrapoints = 0
    for num in a:
        try:
            if (num-lowerreallimit) < 0:
                extrapoints += 1
            else:
                bintoincrement = int((num-lowerreallimit) / float(binsize))
                bins[bintoincrement] = bins[bintoincrement] + 1
        except:                           # point outside lower/upper limits
            extrapoints += 1
    if extrapoints > 0 and printextras:
        # fixme: warnings.warn()
        print '\nPoints outside given histogram range =',extrapoints
    return (bins, lowerreallimit, binsize, extrapoints)


def cumfreq(a, numbins=10, defaultreallimits=None):
    """
Returns a cumulative frequency histogram, using the histogram function.
Defaultreallimits can be None (use all data), or a 2-sequence containing
lower and upper limits on values to include.

Returns: array of cumfreq bin values, lowerreallimit, binsize, extrapoints
"""
    h,l,b,e = histogram(a,numbins,defaultreallimits)
    cumhist = np.cumsum(h*1, axis=0)
    return cumhist,l,b,e

def abut(source, *args):
    # comment: except for the repetition, this is equivalent to hstack.
    """\nLike the |Stat abut command.  It concatenates two arrays column-wise
    and returns the result.  CAUTION:  If one array is shorter, it will be
    repeated until it is as long as the other.

    Format:  abut (source, args)    where args=any # of arrays
    Returns: an array as long as the LONGEST array past, source appearing on the
    'left', arrays in <args> attached on the 'right'.\n"""

    source = np.asarray(source)
    if len(source.shape)==1:
        width = 1
        source = np.resize(source,[source.shape[0],width])
    else:
        width = source.shape[1]
    for addon in args:
        if len(addon.shape)==1:
            width = 1
            addon = np.resize(addon,[source.shape[0],width])
        else:
            width = source.shape[1]
        if len(addon) < len(source):
            addon = np.resize(addon,[source.shape[0],addon.shape[1]])
        elif len(source) < len(addon):
            source = np.resize(source,[addon.shape[0],source.shape[1]])
        source = np.concatenate((source,addon),1)
    return source
def unique(inarray):
    """Returns unique items in the FIRST dimension of the passed array. Only
    works on arrays NOT including string items (e.g., type 'O' or 'c').
    """
    inarray = np.asarray(inarray)
    uniques = np.array([inarray[0]])
    if len(uniques.shape) == 1:            # IF IT'S A 1D ARRAY
        for item in inarray[1:]:
            if np.add.reduce(np.equal(uniques,item).flat) == 0:
                try:
                    uniques = np.concatenate([uniques,np.array[np.newaxis,:]])
                except TypeError:
                    uniques = np.concatenate([uniques,np.array([item])])
    else:                                  # IT MUST BE A 2+D ARRAY
        if inarray.dtype.char != 'O':  # not an Object array
            for item in inarray[1:]:
                if not np.sum(np.alltrue(np.equal(uniques,item),1),axis=0):
                    try:
                        uniques = np.concatenate( [uniques,item[np.newaxis,:]] )
                    except TypeError:    # the item to add isn't a list
                        uniques = np.concatenate([uniques,np.array([item])])
                else:
                    pass  # this item is already in the uniques array
        else:   # must be an Object array, alltrue/equal functions don't work
            for item in inarray[1:]:
                newflag = 1
                for unq in uniques:  # NOTE: cmp --> 0=same, -1=<, 1=>
                    test = np.sum(abs(np.array(map(cmp,item,unq))),axis=0)
                    if test == 0:   # if item identical to any 1 row in uniques
                        newflag = 0 # then not a novel item to add
                        break
                if newflag == 1:
                    try:
                        uniques = np.concatenate( [uniques,item[np.newaxis,:]] )
                    except TypeError:    # the item to add isn't a list
                        uniques = np.concatenate([uniques,np.array([item])])
    return uniques  

          
def itemfreq(a):
    # fixme: I'm not sure I understand what this does. The docstring is
    # internally inconsistent.
    # comment: fortunately, this function doesn't appear to be used elsewhere
    """Returns a 2D array of item frequencies.

    Column 1 contains item values, column 2 contains their respective counts.
    Assumes a 1D array is passed.

    Parameters
    ----------
    a : array

    Returns
    -------
    A 2D frequency table (col [0:n-1]=scores, col n=frequencies)
    """
    scores = unique(a)
    scores = np.sort(scores)
    freq = np.zeros(len(scores))
    for i in range(len(scores)):
        freq[i] = np.add.reduce(np.equal(a,scores[i]))
        
    return np.array(abut(scores, freq))



def dranker(config):
    """Connect the given tag with the data from the kb file given"""
    write_message("Loading knowledgebase file", verbose=9)
    kb_data = {}
    records = []

    write_message("Reading knowledgebase file: %s" % \
                   config.get(config.get("rank_method", "function"), "kb_src"))
    input = open(config.get(config.get("rank_method", "function"), "kb_src"), 'r')
    data = input.readlines()
    for line in data:
        if not line[0:1] == "#":
            kb_data[string.strip((string.split(string.strip(line), "---"))[0])] = (string.split(string.strip(line), "---"))[1]
    write_message("Number of lines read from knowledgebase file: %s" % len(kb_data))

    tag = config.get(config.get("rank_method", "function"), "tag")
    tags = config.get(config.get("rank_method", "function"), "check_mandatory_tags").split(", ")
    if tags == ['']:
        tags = ""

    records = []
    for (recids, recide) in options["recid_range"]:
        task_sleep_now_if_required(can_stop_too=True)
        write_message("......Processing records #%s-%s" % (recids, recide))
        recs = run_sql("SELECT id_bibrec, value FROM bib%sx, bibrec_bib%sx WHERE tag=%%s AND id_bibxxx=id and id_bibrec >=%%s and id_bibrec<=%%s" % (tag[0:2], tag[0:2]), (tag, recids, recide))
        valid = HitSet(trailing_bits=1)
        valid.discard(0)
        for key in tags:
            newset = HitSet()
            newset += [recid[0] for recid in (run_sql("SELECT id_bibrec FROM bib%sx, bibrec_bib%sx WHERE id_bibxxx=id AND tag=%%s AND id_bibxxx=id and id_bibrec >=%%s and id_bibrec<=%%s" % (tag[0:2], tag[0:2]), (key, recids, recide)))]
            valid.intersection_update(newset)
        if tags:
            recs = filter(lambda x: x[0] in valid, recs)
        records = records + list(recs)
        write_message("Number of records found with the necessary tags: %s" % len(records))

    records = filter(lambda x: x[0] in options["validset"], records)
    rnkset = {}
    for key, value in records:
        if kb_data.has_key(value):
            if not rnkset.has_key(key):
                rnkset[key] = float(kb_data[value])
            else:
                if kb_data.has_key(rnkset[key]) and float(kb_data[value]) > float((rnkset[key])[1]):
                    rnkset[key] = float(kb_data[value])
        else:
            rnkset[key] = 0

    write_message("Number of records available in rank method: %s" % len(rnkset))
    return rnkset  

def del_rank_method_codeDATA(rank_method_code):
    """Delete the data for a rank method"""
    
    id = run_sql("SELECT id from rnkMETHOD where name=%s", (rank_method_code, ))
    run_sql("DELETE FROM rnkMETHODDATA WHERE id_rnkMETHOD=%s", (id[0][0], ))
    
def fromDB(rank_method_code):
    """Get the data for a rank method"""
    id = run_sql("SELECT id from rnkMETHOD where name=%s", (rank_method_code, ))
    res = run_sql("SELECT relevance_data FROM rnkMETHODDATA WHERE id_rnkMETHOD=%s", (id[0][0], ))
    if res:
        return deserialize_via_marshal(res[0][0])
    else:
        return {}

def intoDB(dict, last_updated, rank_method_code):
    """Insert the rank method data into the database"""
    mid = run_sql("SELECT id from rnkMETHOD where name=%s", (rank_method_code, ))
    del_rank_method_codeDATA(rank_method_code)
    serdata = serialize_via_marshal(dict);
    midstr = str(mid[0][0]);
    run_sql("INSERT INTO rnkMETHODDATA(id_rnkMETHOD, relevance_data) VALUES (%s,%s)", (midstr, serdata,))
    run_sql("UPDATE rnkMETHOD SET last_updated=%s WHERE name=%s", (last_updated, rank_method_code))  

def intoDRANKDB(drank_dictionary_table,id,name,description,dict_type,dictionary,last_updated):
    """Insert the drank method data into the database"""

    run_sql("DELETE FROM %s WHERE name=%%s and dictionary_type=%%s" % drank_dictionary_table,(name,dict_type,))
    serialized_dictionary = serialize_via_marshal(dictionary)
    run_sql("""INSERT INTO %s(id,name,description,dictionary_type,dictionary,last_updated) VALUES (%%s,%%s,%%s,%%s,%%s,%%s)"""%\
             drank_dictionary_table, (id,name,description,dict_type,serialized_dictionary,last_updated,))

def insert_dict_intodb(drank_dictionary_table,id,name,description,dict_type,dictionary,dmethod):
    """ """
    last_updated = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    
    if dmethod=="drank":
        intoDRANKDB(drank_dictionary_table,id,name,description,dict_type,dictionary,last_updated)
    elif dmethod=="bibrank":
        intoDB(dictionary,last_updated,"drank") 
   
def drank_indexer(rank_method_code):  
    """ The function accepts as input the name given to the config file during setup, 
    mostly, drank, and calls all offline based calculations """  
    
    #repopulte data from queryanalyzer
    repopulate_data_from_queryanalyzer()
    
    #get configuration
    drank_config=get_config(rank_method_code)
    print "config:{0}".format(drank_config)
    
    #TODO: put the following into a function
    drank_dictionary_table=drank_config["drank_parameters"]["drank_dictionary_table"]
    counts_table=drank_config["drank_parameters"]["counts_table"]   
    default_method=drank_config["drank_parameters"]["default_method"]
    max_nbr_of_displays=drank_config["drank_parameters"]["max_nbr_of_displays"]        
    weight_of_relevance_exposed=drank_config["drank_parameters"]["weight_of_relevance_exposed"]
    weight_of_relevance_fresh=drank_config["drank_parameters"]["weight_of_relevance_fresh"]
    weight_of_exposed=drank_config["drank_parameters"]["weight_of_exposed"]
    weight_of_fresh=drank_config["drank_parameters"]["weight_of_fresh"]
    percentage_of_freshness=drank_config["drank_parameters"]["percentage_of_freshness"]
    min_freshness_value=drank_config["drank_parameters"]["min_freshness_value"]

    run_sql("DELETE FROM %s" % drank_dictionary_table)
    n=0
#    res=np.array(run_sql("SELECT nbr_pageviews,nbr_downloads,nbr_seens,nbr_displays,freshness,nbr_pageviews_decayed,nbr_downloads_decayed,nbr_seens_decayed,nbr_displays_decayed FROM %s" %counts_table ))
#    res=np.array(run_sql("""SELECT nbr_pageviews,nbr_downloads,nbr_seens,nbr_displays,freshness,nbr_pageviews_decayed,nbr_downloads_decayed,nbr_seens_decayed,nbr_displays_decayed FROM %s""" %drank_dictionary_table))
    
    views=run_sql("SELECT nbr_pageviews FROM %s" %counts_table )
    views=get_percentiles_from_scores(views,1)
    print "views:{0}".format(views)
    n=n+1
    insert_dict_intodb(drank_dictionary_table,n,"views","look up table for views","views",views,rank_method_code)

    downloads=run_sql("SELECT nbr_downloads FROM %s" %counts_table )
    downloads=get_percentiles_from_scores(downloads,1)
    n=n+1
    insert_dict_intodb(drank_dictionary_table,n,"downloads","look up table for downloads","downloads",downloads,rank_method_code)
    
    seens=run_sql("SELECT nbr_seens FROM %s" %counts_table )
    seens=get_percentiles_from_scores(seens,1)
    n=n+1
    insert_dict_intodb(drank_dictionary_table,n,"seens","look up table for seens","seens",seens,rank_method_code)
    
    displays=run_sql("SELECT nbr_displays FROM %s" %counts_table )    
    displays=get_percentiles_from_scores(displays,1)
    n=n+1
    insert_dict_intodb(drank_dictionary_table,n,"displays","look up table for displays ","displays",displays,rank_method_code)
    
    freshness=run_sql("SELECT freshness FROM %s" %counts_table )
    freshness=get_percentiles_from_scores(freshness,1)
    n=n+1
    insert_dict_intodb(drank_dictionary_table,n,"freshness","look up table for freshness","freshness",freshness,rank_method_code)
    
    pageviews_decayed=run_sql("SELECT nbr_pageviews_decayed FROM %s" %counts_table )
    pageviews_decayed=get_percentiles_from_scores(pageviews_decayed,1)
    n=n+1
    insert_dict_intodb(drank_dictionary_table,n,"views_decayed","look up table for pageviews decayed","pageviews_decayed",pageviews_decayed,rank_method_code)
    
    downloads_decayed=run_sql("SELECT nbr_downloads_decayed FROM %s" %counts_table )    
    downloads_decayed=get_percentiles_from_scores(downloads_decayed,1)
    n=n+1
    insert_dict_intodb(drank_dictionary_table,n,"downloads_decayed","look up table for downloads decayed","downloads_decayed",downloads_decayed,rank_method_code)
    
    seens_decayed=run_sql("SELECT nbr_seens_decayed FROM %s" %counts_table )    
    seens_decayed=get_percentiles_from_scores(seens_decayed,1)
    n=n+1
    insert_dict_intodb(drank_dictionary_table,n,"seens_decayed","look up table for seens decayed","seens_decayed",seens_decayed,rank_method_code)

    displays_decayed=run_sql("SELECT nbr_displays_decayed FROM %s" %counts_table )
    displays_decayed =get_percentiles_from_scores(displays_decayed,1)
    n=n+1
    insert_dict_intodb(drank_dictionary_table,n,"displays_decayed","look up table for displays decayed","displays_decayed",displays_decayed,rank_method_code)

    exposed_records, scores_fresh_records=get_quality_records(max_nbr_of_displays)
    print "-----------------------------------------------------------------------------------------------------------------------------"
    freshness_percentile_scores, lut_fresh =get_percentiles_from_scores(scores_fresh_records)
    print "lut{0}".format(lut_fresh)
    print "fresh{0}".format(freshness_percentile_scores)
    print "-----------------------------------------------------------------------------------------------------------------------------"
     
    n=n+1
    insert_dict_intodb(drank_dictionary_table,n,"fresh","scores for less exposed records","fresh",freshness_percentile_scores,rank_method_code)
    n=n+1
    insert_dict_intodb(drank_dictionary_table,n,"fresh","look up table for less exposed records","fresh_lut",lut_fresh,rank_method_code)
    
    wrd_sim_scores=fromDB("wrd")
    wrd_sim_percentile_scores, lut_sim=get_percentiles_from_scores(wrd_sim_scores)
    print """============================================================================================="""
    print "wrd_sim_sc:{0}".format(wrd_sim_scores)
    print "lut{0}".format(lut_sim)
    print "wrd_sim_perc:{0}".format(wrd_sim_percentile_scores)
    print "-----------------------------------------------------------------------------------------------------------------------------"
    n=n+1
    insert_dict_intodb(drank_dictionary_table,n,"wrd","word similarity scores","wrd_lut",lut_sim,rank_method_code)
    
    for id,drank_method in drank_config["drank_methods"].iteritems():        
        exposed_scores=generate_new_scores(exposed_records,drank_method,drank_dictionary_table,drank_config)
        exposed_percentile_scores, lut_exposed =get_percentiles_from_scores(exposed_scores)
        print """============================================================================================="""
        print "lut{0}".format(lut_exposed)
        print "exposed{0}".format(exposed_percentile_scores)
        

        print """============================================================================================="""
        description=drank_config[drank_method]["description"] 
        n=n+1
        insert_dict_intodb(drank_dictionary_table,n,drank_method,description,"exposed",exposed_percentile_scores,rank_method_code)
        n=n+1
        insert_dict_intodb(drank_dictionary_table,n,drank_method,"look up table for "+description,"exposed_lut",lut_exposed,rank_method_code)

def get_record_relevance_for_drank(recdict, rec_termcount, hitset, rank_limit_relevance, verbose):
    """Sorts the dictionary and returns records with a relevance higher than the given value.
    recdict - {recid: value} unsorted
    rank_limit_relevance - a value > 0 usually
    verbose - verbose value"""

    startCreate = time.time()
    global voutput
    reclist = []

    #Multiply with the number of terms of the total number of terms in the query existing in the records
    for j in recdict.keys():
        if recdict[j] > 0 and rec_termcount[j] > 1:
            recdict[j] = math.log((recdict[j] * rec_termcount[j]))
        else:
            recdict[j] = 0

    hitset -= recdict.keys()
    #gives each record a score between 0-100
    divideby = max(recdict.values())
    for (j, w) in recdict.iteritems():
#        w = int(w * 100 / divideby)
        if w >= rank_limit_relevance:
            reclist.append((j, w))

    #sort scores
    reclist.sort(lambda x, y: cmp(x[1], y[1]))

    if verbose > 0:
        voutput += "Number of records sorted: %s<br />" % len(reclist)
        voutput += "Sort time: %s<br />" % (str(time.time() - startCreate))
    return (reclist, hitset)

def sort_record_relevance_findsimilar_for_drank(recdict, rec_termcount, hitset, rank_limit_relevance, verbose):
    """Sorts the dictionary and returns records with a relevance higher than the given value.
    recdict - {recid: value} unsorted
    rank_limit_relevance - a value > 0 usually
    verbose - verbose value"""

    startCreate = time.time()
    global voutput
    reclist = []

    #Multiply with the number of terms of the total number of terms in the query existing in the records
    for j in recdict.keys():
        if recdict[j] > 0 and rec_termcount[j] > 1:
            recdict[j] = math.log((recdict[j] * rec_termcount[j]))
        else:
            recdict[j] = 0

    hitset -= recdict.keys()
    #gives each record a score between 0-100
    divideby = max(recdict.values())
    for (j, w) in recdict.iteritems():
#        w = int(w * 100 / divideby)
        if w >= rank_limit_relevance/100: ##TODO: check the range of values of w and set the appropriate treshhold
            reclist.append((j, w))

    #sort scores
    reclist.sort(lambda x, y: cmp(x[1], y[1]))

    if verbose > 0:
        voutput += "Number of records sorted: %s<br />" % len(reclist)
        voutput += "Sort time: %s<br />" % (str(time.time() - startCreate))
#    return (reclist, hitset)
    return dict(reclist) ##return dictionary

def find_similar_for_drank(drank_methods, recID, hitset, rank_limit_relevance,verbose):
    """Finding terms to use for calculating similarity. Terms are taken from the recid given, returns a list of recids's and relevance,
    input:
    rank_method_code - the code of the method, from the name field in rnkMETHOD
    recID - records to use for find similar
    hitset - a list of hits for the query found by search_engine
    rank_limit_relevance - show only records with a rank value above this
    verbose - verbose value
    output:
    reclist - a list of sorted records: [[23,34], [344,24], [1,01]]
    prefix - what to show before the rank value
    postfix - what to show after the rank value
    voutput - contains extra information, content dependent on verbose value"""

    startCreate = time.time()
    global voutput

    if verbose > 0:
        voutput += "<br />Drank: Running rank method: %s, using find_similar/word_frequency in bibrank_record_sorter<br />"
    if not rank_limit_relevance:
        rank_limit_relevance = drank_methods["relevance_parameters"]["default_min_relevance"]

    try:
        recID = int(recID)
    except Exception,e :
        return (None, "Warning: Error in record ID, please check that a number is given.", "", voutput)

    rec_terms = run_sql("""SELECT termlist FROM %sR WHERE id_bibrec=%%s""" % drank_methods["wrd"]["rnkWORD_table"][:-1],  (recID,))
    if not rec_terms:
        return (None, "Warning: Requested record does not seem to exist.", "", voutput)
    rec_terms = deserialize_via_marshal(rec_terms[0][0])

    #Get all documents using terms from the selected documents
    if len(rec_terms) == 0:
        return (None, "Warning: Record specified has no content indexed for use with this method.", "", voutput)
    else:
        terms = "%s" % rec_terms.keys()
        terms_recs = dict(run_sql("""SELECT term, hitlist FROM %s WHERE term IN (%s)""" % (drank_methods["wrd"]["rnkWORD_table"], terms[1:len(terms) - 1])))

    tf_values = {}
    #Calculate all term frequencies
    for (term, tf) in rec_terms.iteritems():
        if len(term) >= drank_methods["relevance_parameters"]["min_word_length"] and terms_recs.has_key(term) and tf[1] != 0:
            tf_values[term] =  int((1 + math.log(tf[0])) * tf[1]) #calculate term weigth
    tf_values = tf_values.items()
    tf_values.sort(lambda x, y: cmp(y[1], x[1])) #sort based on weigth

    pattern = []
    stime = time.time()
    (recdict, rec_termcount) = ({}, {})

    for (t, tf) in tf_values: #t=term, tf=term frequency
        term_recs = deserialize_via_marshal(terms_recs[t])
        if len(tf_values) <=\
         drank_methods["relevance_parameters"]["max_nr_words_lower"] or (len(term_recs) >= drank_methods["relevance_parameters"]["min_nr_words_docs"] and (((float(len(term_recs)) / float(drank_methods["wrd"]["col_size"])) <=  drank_methods["wrd"]["max_word_occurence"]) and ((float(len(term_recs)) / float(drank_methods["wrd"]["col_size"])) >= drank_methods["wrd"]["min_word_occurence"]))): #too complicated...something must be done
            pattern.append((t, drank_methods["wrd"]["rnkWORD_table"])) #list of terms used
            (recdict, rec_termcount) = get_record_relevance_for_drank((t, round(tf, 4)) , term_recs, hitset, recdict, rec_termcount, verbose, "true") #true tells the function to not calculate all unimportant terms
        if len(tf_values) > drank_methods["relevance_parameters"]["max_nr_words_lower"] and (len(pattern) ==  drank_methods["relevance_parameters"]["max_nr_words_upper"] or tf < 0):
            break

    if len(recdict) == 0 or len(pattern) == 0:
        return (None, "Could not find any similar documents, possibly because of error in ranking data.", "", voutput)
    else: #sort if we got something to sort
        (reclist, hitset) = sort_record_relevance_findsimilar_for_drank(recdict, rec_termcount, hitset, rank_limit_relevance, verbose)

    if verbose > 0:
        voutput += "<br />Number of terms: %s<br />" % run_sql("SELECT count(id) FROM %s" % drank_methods["wrd"]["rnkWORD_table"])[0][0]
        voutput += "Number of terms to use for query: %s<br />" % len(pattern)
        voutput += "Terms: %s<br />" % pattern
        voutput += "Current number of recIDs: %s<br />" % (drank_methods["wrd"]["col_size"])
        voutput += "Prepare time: %s<br />" % (str(time.time() - startCreate))
        voutput += "Total time used: %s<br />" % (str(time.time() - startCreate))
#        rank_method_stat(rank_method_code, reclist, pattern)
    
#    return (reclist[:len(reclist)], drank_methods["wrd"]["prefix"], drank_methods["wrd"]["postfix"], voutput)
    return (reclist[:len(reclist)], voutput)

        
def drank_sorter(rank_method_code,pattern,hitset,rank_limit_relevance,verbose):                
    """Ranking of records based on predetermined values.
    input:
    rank_method_code - the code of the method, from the name field in rnkMETHOD, used to get predetermined values from
    rnkMETHODDATA
    pattern - a list of words from the query
    hitset - a list of hits for the query found by search_engine
    rank_limit_relevance - show only records with a rank value above this
    verbose - verbose value
    output:
    reclist - a list of sorted records, with unsorted added to the end: [[23,34], [344,24], [1,01]]
    prefix - what to show before the rank value
    postfix - what to show after the rank value
    voutput - contains extra information, content dependent on verbose value"""
    
    global voutput
    voutput = ""
    drank_config=get_config(rank_method_code)
    
    default_method=         drank_config["drank_parameters"]["default_method"]
    description=           drank_config[default_method]["description"]
    quality_exposed = run_sql("SELECT dictionary FROM rnkDRANKLUT where dictionary_type='%s' and name='%s'" % ("exposed",default_method,))
    quality_fresh = run_sql("SELECT dictionary FROM rnkDRANKLUT where dictionary_type='%s' and name='%s'" % ("fresh","fresh",))
    if not quality_exposed:
        return (None, "Warning: Could not load ranking data for method %s." % default_method, "", voutput)

    max_recid = 0
    res = run_sql("SELECT max(id) FROM bibrec")
    if res and res[0][0]:
        max_recid = int(res[0][0])

    lwords_hitset = None
    for j in range(0, len(pattern)): #find which docs to search based on ranges..should be done in search_engine...
        if pattern[j] and pattern[j][:6] == "recid:":
            if not lwords_hitset:
                lwords_hitset = intbitset()
            lword = pattern[j][6:]
            if string.find(lword, "->") > -1:
                lword = string.split(lword, "->")
                if int(lword[0]) >= max_recid or int(lword[1]) >= max_recid + 1:
                    return (None, "Warning: Given record IDs are out of range.", "", voutput)
                for i in range(int(lword[0]), int(lword[1])):
                    lwords_hitset.add(int(i))
            elif lword < max_recid + 1:
                lwords_hitset.add(int(lword))
            else:
                return (None, "Warning: Given record IDs are out of range.", "", voutput)
    relevance_wrdsim={}

    relevance_wrdsim = find_similar_for_drank(drank_config, pattern[0][6:], hitset, rank_limit_relevance, verbose)
##    records, scores = wrdsim[0:2]
##    for k,v in enumerate(records):
##        relevance_wrdsim[v]=scores[k]

              
    relevance_wrdsim=lookup_percentile_scores(relevance_wrdsim,"wrd_lut","wrd")     
    quality_exposed = deserialize_via_marshal(quality_exposed[0][0])   
    quality_fresh = deserialize_via_marshal(quality_fresh[0][0])
    
    quality_exposed_relevance=merge_relevant_quality(relevance_wrdsim,quality_exposed)
    quality_fresh_relevance=merge_relevant_quality(relevance_wrdsim,quality_fresh)
    
    dictionary=merge_exposed_fresh_dictionary(quality_exposed_relevance,quality_fresh_relevance)
     
    
    if verbose > 0:
        voutput += "<br />Running Ranking Method: Distributed Ranking using the %s Configuration<br />" % description
        voutput += "Ranking data loaded, size of structure: %s<br />" % len(dictionary)
        voutput += "wrdsym data loaded, size of structure: %s<br />" % str(wrdsim)
        voutput += "wrd_sym data loaded, size of structure: %s<br />" % str(relevance_wrdsim)
    lrecIDs = list(hitset)

    if verbose > 0:
        voutput += "Number of records to rank: %s<br />" % len(lrecIDs)
    reclist = []
    reclist_addend = []

    if not lwords_hitset: #rank all docs, can this be speed up using something else than for loop?
        for recID in lrecIDs:
            if dictionary.has_key(recID):
                reclist.append((recID, dictionary[recID]))
                del dictionary[recID]
            else:
                reclist_addend.append((recID, 0))
    else: #rank docs in hitset, can this be speed up using something else than for loop?
        for recID in lwords_hitset:
            if dictionary.has_key(recID) and recID in hitset:
                reclist.append((recID, dictionary[recID]))
                del dictionary[recID]
            elif recID in hitset:
                reclist_addend.append((recID, 0))

    if verbose > 0:
        voutput += "Number of records ranked: %s<br />" % len(reclist)
        voutput += "Number of records not ranked: %s<br />" % len(reclist_addend)
        voutput += "Total Number of records: %s<br />" % len(reclist_addend+reclist)

    reclist.sort(lambda x, y: cmp(x[1], y[1]))
    result=reclist_addend + reclist
#    return result,voutput
    return reclist,voutput

#rank_method_code="drank"
#drank_config=get_config(rank_method_code)
##drank_sorter(rank_method_code,drank_config,pattern, hitset,"",9)
#drank_indexer("drank")
#
#print "done"