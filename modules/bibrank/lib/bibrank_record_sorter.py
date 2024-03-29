# -*- coding: utf-8 -*-
## Ranking of records using different parameters and methods on the fly.
##
## This file is part of Invenio.
## Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011 CERN.
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

__revision__ = "$Id$"

import string
import time
import math
import re
import ConfigParser
import copy

from invenio.config import \
     CFG_SITE_LANG, \
     CFG_ETCDIR
from invenio.dbquery import run_sql, deserialize_via_marshal
from invenio.errorlib import register_exception
from invenio.webpage import adderrorbox
from invenio.bibindex_engine_stemmer import stem
from invenio.bibindex_engine_stopwords import is_stopword
from invenio.bibrank_citation_searcher import get_cited_by, get_cited_by_weight
from invenio.intbitset import intbitset
from invenio.bibrank_dictionary import rnkDict



def rescale(recdict, rank_limit_relevance):
    """Rescale list so that values are between 0-100."""
    reclist = []
    divideby = max(recdict.values())
    for (j, w) in recdict.iteritems():
        w = int(w * 100 / divideby)
        if w >= rank_limit_relevance:
            reclist.append((j, w))
    return reclist

def rescale_list(reclist, rank_limit_relevance):
    """Rescale list so that values are between 0-100."""
    values_new = []
    keys, values = zip(*reclist)
    divideby = max(values)
    for item in values:
        item = int(item * 100 / divideby)
        if item >= rank_limit_relevance:
            values_new.append(item)
    reclist = zip(*(keys, values_new))
    return reclist

def compare_on_val(first, second):
    """Compare on val."""
    return cmp(second[1], first[1])

def check_term(term, col_size, term_rec, max_occ, min_occ, termlength):
    """Check if the tem is valid for use
    term - the term to check
    col_size - the number of records in database
    term_rec - the number of records which contains this term
    max_occ - max frequency of the term allowed
    min_occ - min frequence of the term allowed
    termlength - the minimum length of the terms allowed"""

    try:
        if is_stopword(term, 1) or (len(term) <= termlength) or ((float(term_rec) / float(col_size)) >= max_occ) or ((float(term_rec) / float(col_size)) <= min_occ):
            return ""
        if int(term):
            return ""
    except StandardError, e:
        pass
    return "true"

def create_rnkmethod_cache():
    """Create cache with vital information for each rank method."""

    global methods
    bibrank_meths = run_sql("SELECT name from rnkMETHOD")
    methods = {}
    global voutput
    voutput = ""

    for (rank_method_code,) in bibrank_meths:
        try:
            file = CFG_ETCDIR + "/bibrank/" + rank_method_code + ".cfg"
            config = ConfigParser.ConfigParser()
            config.readfp(open(file))
        except StandardError, e:
            pass

        cfg_function = config.get("rank_method", "function")
        if config.has_section(cfg_function):
            methods[rank_method_code] = {}
            methods[rank_method_code]["function"] = cfg_function
            methods[rank_method_code]["prefix"] = config.get(cfg_function, "relevance_number_output_prologue")
            methods[rank_method_code]["postfix"] = config.get(cfg_function, "relevance_number_output_epilogue")
            methods[rank_method_code]["chars_alphanumericseparators"] = r"[1234567890\!\"\#\$\%\&\'\(\)\*\+\,\-\.\/\:\;\<\=\>\?\@\[\\\]\^\_\`\{\|\}\~]"
        else:
            raise Exception("Error in configuration file: %s" % (CFG_ETCDIR + "/bibrank/" + rank_method_code + ".cfg"))

        i8n_names = run_sql("""SELECT ln,value from rnkMETHODNAME,rnkMETHOD where id_rnkMETHOD=rnkMETHOD.id and rnkMETHOD.name=%s""", (rank_method_code,))
        for (ln, value) in i8n_names:
            methods[rank_method_code][ln] = value

        if config.has_option(cfg_function, "table"):
            methods[rank_method_code]["rnkWORD_table"] = config.get(cfg_function, "table")
            methods[rank_method_code]["col_size"] = run_sql("SELECT count(*) FROM %sR" % methods[rank_method_code]["rnkWORD_table"][:-1])[0][0]

        if config.has_option(cfg_function, "stemming") and config.get(cfg_function, "stemming"):
            try:
                methods[rank_method_code]["stemmer"] = config.get(cfg_function, "stemming")
            except Exception,e:
                pass

        if config.has_option(cfg_function, "stopword"):
            methods[rank_method_code]["stopwords"] = config.get(cfg_function, "stopword")

        if config.has_section("find_similar"):
            methods[rank_method_code]["max_word_occurence"] = float(config.get("find_similar", "max_word_occurence"))
            methods[rank_method_code]["min_word_occurence"] = float(config.get("find_similar", "min_word_occurence"))
            methods[rank_method_code]["min_word_length"] = int(config.get("find_similar", "min_word_length"))
            methods[rank_method_code]["min_nr_words_docs"] = int(config.get("find_similar", "min_nr_words_docs"))
            methods[rank_method_code]["max_nr_words_upper"] = int(config.get("find_similar", "max_nr_words_upper"))
            methods[rank_method_code]["max_nr_words_lower"] = int(config.get("find_similar", "max_nr_words_lower"))
            methods[rank_method_code]["default_min_relevance"] = int(config.get("find_similar", "default_min_relevance"))

        if config.has_section("combine_method"):
            i = 1
            methods[rank_method_code]["combine_method"] = []
            while config.has_option("combine_method", "method%s" % i):
                methods[rank_method_code]["combine_method"].append(string.split(config.get("combine_method", "method%s" % i), ","))
                i += 1

def is_method_valid(colID, rank_method_code):
    """
    Check if RANK_METHOD_CODE method is valid for the collection given.
    If colID is None, then check for existence regardless of collection.
    """

    if colID is None:
        return run_sql("SELECT COUNT(*) FROM rnkMETHOD WHERE name=%s", (rank_method_code,))[0][0]

    enabled_colls = dict(run_sql("SELECT id_collection, score from collection_rnkMETHOD,rnkMETHOD WHERE id_rnkMETHOD=rnkMETHOD.id AND name='%s'" % rank_method_code))

    try:
        colID = int(colID)
    except TypeError:
        return 0

    if enabled_colls.has_key(colID):
        return 1
    else:
        while colID:
            colID = run_sql("SELECT id_dad FROM collection_collection WHERE id_son=%s" % colID)
            if colID and enabled_colls.has_key(colID[0][0]):
                return 1
            elif colID:
                colID = colID[0][0]
    return 0

def get_bibrank_methods(colID, ln=CFG_SITE_LANG):
    """
    Return a list of rank methods enabled for collection colID and the
    name of them in the language defined by the ln parameter.
    """

    if not globals().has_key('methods'):
        create_rnkmethod_cache()

    avail_methods = []
    for (rank_method_code, options) in methods.iteritems():
        if options.has_key("function") and is_method_valid(colID, rank_method_code):
            if options.has_key(ln):
                avail_methods.append((rank_method_code, options[ln]))
            elif options.has_key(CFG_SITE_LANG):
                avail_methods.append((rank_method_code, options[CFG_SITE_LANG]))
            else:
                avail_methods.append((rank_method_code, rank_method_code))
    return avail_methods

def ranked(_rnkdict, hitset):
    """"""
    rnkdict = rnkDict()
    rnkdict.loaddict(_rnkdict)
    rnkdict.clean()
 
    rnkdict.filter(list(hitset))
    result = (rnkdict.rank(),"(", ")", "")
    return result # rnkdict.rank()

def rank_records(rank_method_code, rank_limit_relevance, hitset_global, pattern=[], verbose=0):
    """rank_method_code, e.g. `jif' or `sbr' (word frequency vector model)
       rank_limit_relevance, e.g. `23' for `nbc' (number of citations) or `0.10' for `vec'
       hitset, search engine hits;
       pattern, search engine query or record ID (you check the type)
       verbose, verbose level
       output:
       list of records
       list of rank values
       prefix
       postfix
       verbose_output"""

    global voutput
    voutput = ""
    configcreated = ""

    starttime = time.time()
    afterfind = starttime - time.time()
    aftermap = starttime - time.time()

    try:
        hitset = copy.deepcopy(hitset_global) #we are receiving a global hitset
        if not globals().has_key('methods'):
            create_rnkmethod_cache()

        function = methods[rank_method_code]["function"]
        #we get 'citation' method correctly here
        func_object = globals().get(function)

        if func_object and pattern and pattern[0][0:6] == "recid:" and function == "word_similarity":
            result = find_similar(rank_method_code, pattern[0][6:], hitset, rank_limit_relevance, verbose)
        elif rank_method_code == "download":
            result = ranked(rank_method_code, hitset)
        elif rank_method_code == "citation":
            #we get rank_method_code correctly here. pattern[0] is the search word - not used by find_cit
            p = ""
            if pattern and pattern[0]:
                p = pattern[0][6:]
            result = find_citations(rank_method_code, p, hitset, verbose)

        elif func_object:
            result = func_object(rank_method_code, pattern, hitset, rank_limit_relevance, verbose)
        else:
            result = rank_by_method(rank_method_code, pattern, hitset, rank_limit_relevance, verbose)
    except Exception, e:
        register_exception()
        result = (None, "", adderrorbox("An error occured when trying to rank the search result "+rank_method_code, ["Unexpected error: %s<br />" % (e,)]), voutput)

    afterfind = time.time() - starttime

    if result[0] and result[1]: #split into two lists for search_engine
        results_similar_recIDs = map(lambda x: x[0], result[0])
        results_similar_relevances = map(lambda x: x[1], result[0])
        result = (results_similar_recIDs, results_similar_relevances, result[1], result[2], "%s" % configcreated + result[3])
        aftermap = time.time() - starttime;
    else:
        result = (None, None, result[1], result[2], result[3])

    if verbose > 0:
        voutput = voutput+"\nElapsed time after finding: "+str(afterfind)+"\nElapsed after mapping: "+str(aftermap)

    #add stuff from here into voutput from result
    tmp = result[4]+voutput
    result = (result[0],result[1],result[2],result[3],tmp)

    #dbg = string.join(map(str,methods[rank_method_code].items()))
    #result = (None, "", adderrorbox("Debug ",rank_method_code+" "+dbg),"",voutput);
    return result

def combine_method(rank_method_code, pattern, hitset, rank_limit_relevance,verbose):
    """combining several methods into one based on methods/percentage in config file"""

    global voutput
    result = {}
    try:
        for (method, percent) in methods[rank_method_code]["combine_method"]:
            function = methods[method]["function"]
            func_object = globals().get(function)
            percent = int(percent)

            if func_object:
                this_result = func_object(method, pattern, hitset, rank_limit_relevance, verbose)[0]
            else:
                this_result = rank_by_method(method, pattern, hitset, rank_limit_relevance, verbose)[0]

            for i in range(0, len(this_result)):
                (recID, value) = this_result[i]
                if value > 0:
                    result[recID] = result.get(recID, 0) + int((float(i) / len(this_result)) * float(percent))

        result = result.items()
        result.sort(lambda x, y: cmp(x[1], y[1]))
        return (result, "(", ")", voutput)
    except Exception, e:
        return (None, "Warning: %s method cannot be used for ranking your query." % rank_method_code, "", voutput)

def rank_by_method(rank_method_code, lwords, hitset, rank_limit_relevance,verbose):
    """Ranking of records based on predetermined values.
    input:
    rank_method_code - the code of the method, from the name field in rnkMETHOD, used to get predetermined values from
    rnkMETHODDATA
    lwords - a list of words from the query
    hitset - a list of hits for the query found by search_engine
    rank_limit_relevance - show only records with a rank value above this
    verbose - verbose value
    output:
    reclist - a list of sorted records, with unsorted added to the end: [[23,34], [344,24], [1,01]]
    prefix - what to show before the rank value
    postfix - what to show after the rank value
    voutput - contains extra information, content dependent on verbose value"""

    global voutput
    rnkdict = run_sql("SELECT relevance_data FROM rnkMETHODDATA,rnkMETHOD where rnkMETHOD.id=id_rnkMETHOD and rnkMETHOD.name='%s'" % rank_method_code)

    if not rnkdict:
        return (None, "Warning: Could not load ranking data for method %s." % rank_method_code, "", voutput)

    max_recid = 0
    res = run_sql("SELECT max(id) FROM bibrec")
    if res and res[0][0]:
        max_recid = int(res[0][0])

    lwords_hitset = None
    for j in range(0, len(lwords)): #find which docs to search based on ranges..should be done in search_engine...
        if lwords[j] and lwords[j][:6] == "recid:":
            if not lwords_hitset:
                lwords_hitset = intbitset()
            lword = lwords[j][6:]
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

    rnkdict = deserialize_via_marshal(rnkdict[0][0])
    if verbose > 0:
        voutput += "<br />Running rank method: %s, using rank_by_method function in bibrank_record_sorter<br />" % rank_method_code
        voutput += "Ranking data loaded, size of structure: %s<br />" % len(rnkdict)
    lrecIDs = list(hitset)

    if verbose > 0:
        voutput += "Number of records to rank: %s<br />" % len(lrecIDs)
    reclist = []
    reclist_addend = []

    if not lwords_hitset: #rank all docs, can this be speed up using something else than for loop?
        for recID in lrecIDs:
            if rnkdict.has_key(recID):
                reclist.append((recID, rnkdict[recID]))
                del rnkdict[recID]
            else:
                reclist_addend.append((recID, 0))
    else: #rank docs in hitset, can this be speed up using something else than for loop?
        for recID in lwords_hitset:
            if rnkdict.has_key(recID) and recID in hitset:
                reclist.append((recID, rnkdict[recID]))
                del rnkdict[recID]
            elif recID in hitset:
                reclist_addend.append((recID, 0))

    if verbose > 0:
        voutput += "Number of records ranked: %s<br />" % len(reclist)
        voutput += "Number of records not ranked: %s<br />" % len(reclist_addend)

    reclist.sort(lambda x, y: cmp(x[1], y[1]))
    return (reclist_addend + reclist, methods[rank_method_code]["prefix"], methods[rank_method_code]["postfix"], voutput)

def find_citations(rank_method_code, recID, hitset, verbose):
    """Rank by the amount of citations."""
    #calculate the cited-by values for all the members of the hitset
    #returns: ((recordid,weight),prefix,postfix,message)

    global voutput
    voutput = ""

    #If the recID is numeric, return only stuff that cites it. Otherwise return
    #stuff that cites hitset

    #try to convert to int
    recisint = True
    recidint = 0
    try:
        recidint = int(recID)
    except:
        recisint = False
    ret = []
    if recisint:
        myrecords = get_cited_by(recidint) #this is a simple list
        ret = get_cited_by_weight(myrecords)
    else:
        ret = get_cited_by_weight(hitset)
    ret.sort(lambda x,y:cmp(x[1],y[1]))      #ascending by the second member of the tuples

    if verbose > 0:
        voutput = voutput+"\nrecID "+str(recID)+" is int: "+str(recisint)+" hitset "+str(hitset)+"\n"+"find_citations retlist "+str(ret)

    #voutput = voutput + str(ret)

    if ret:
        return (ret,"(", ")", "")
    else:
        return ((),"", "", "")

def find_similar(rank_method_code, recID, hitset, rank_limit_relevance,verbose):
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
        voutput += "<br />Running rank method: %s, using find_similar/word_frequency in bibrank_record_sorter<br />" % rank_method_code
    rank_limit_relevance = methods[rank_method_code]["default_min_relevance"]

    try:
        recID = int(recID)
    except Exception,e :
        return (None, "Warning: Error in record ID, please check that a number is given.", "", voutput)

    rec_terms = run_sql("""SELECT termlist FROM %sR WHERE id_bibrec=%%s""" % methods[rank_method_code]["rnkWORD_table"][:-1],  (recID,))
    if not rec_terms:
        return (None, "Warning: Requested record does not seem to exist.", "", voutput)
    rec_terms = deserialize_via_marshal(rec_terms[0][0])

    #Get all documents using terms from the selected documents
    if len(rec_terms) == 0:
        return (None, "Warning: Record specified has no content indexed for use with this method.", "", voutput)
    else:
        terms = "%s" % rec_terms.keys()
        terms_recs = dict(run_sql("""SELECT term, hitlist FROM %s WHERE term IN (%s)""" % (methods[rank_method_code]["rnkWORD_table"], terms[1:len(terms) - 1])))

    tf_values = {}
    #Calculate all term frequencies
    for (term, tf) in rec_terms.iteritems():
        if len(term) >= methods[rank_method_code]["min_word_length"] and terms_recs.has_key(term) and tf[1] != 0:
            tf_values[term] =  int((1 + math.log(tf[0])) * tf[1]) #calculate term weigth
    tf_values = tf_values.items()
    tf_values.sort(lambda x, y: cmp(y[1], x[1])) #sort based on weigth

    lwords = []
    stime = time.time()
    (recdict, rec_termcount) = ({}, {})

    for (t, tf) in tf_values: #t=term, tf=term frequency
        term_recs = deserialize_via_marshal(terms_recs[t])
        if len(tf_values) <= methods[rank_method_code]["max_nr_words_lower"] or (len(term_recs) >= methods[rank_method_code]["min_nr_words_docs"] and (((float(len(term_recs)) / float(methods[rank_method_code]["col_size"])) <=  methods[rank_method_code]["max_word_occurence"]) and ((float(len(term_recs)) / float(methods[rank_method_code]["col_size"])) >= methods[rank_method_code]["min_word_occurence"]))): #too complicated...something must be done
            lwords.append((t, methods[rank_method_code]["rnkWORD_table"])) #list of terms used
            (recdict, rec_termcount) = calculate_record_relevance_findsimilar((t, round(tf, 4)) , term_recs, hitset, recdict, rec_termcount, verbose, "true") #true tells the function to not calculate all unimportant terms
        if len(tf_values) > methods[rank_method_code]["max_nr_words_lower"] and (len(lwords) ==  methods[rank_method_code]["max_nr_words_upper"] or tf < 0):
            break

    if len(recdict) == 0 or len(lwords) == 0:
        return (None, "Could not find any similar documents, possibly because of error in ranking data.", "", voutput)
    else: #sort if we got something to sort
        (reclist, hitset) = sort_record_relevance_findsimilar(recdict, rec_termcount, hitset, rank_limit_relevance, verbose)

    if verbose > 0:
        voutput += "<br />Number of terms: %s<br />" % run_sql("SELECT count(id) FROM %s" % methods[rank_method_code]["rnkWORD_table"])[0][0]
        voutput += "Number of terms to use for query: %s<br />" % len(lwords)
        voutput += "Terms: %s<br />" % lwords
        voutput += "Current number of recIDs: %s<br />" % (methods[rank_method_code]["col_size"])
        voutput += "Prepare time: %s<br />" % (str(time.time() - startCreate))
        voutput += "Total time used: %s<br />" % (str(time.time() - startCreate))
        rank_method_stat(rank_method_code, reclist, lwords)

    return (reclist[:len(reclist)], methods[rank_method_code]["prefix"], methods[rank_method_code]["postfix"], voutput)

def enhanced_ranking(rank_method_code, lwords, hitset, rank_limit_relevance, verbose):
    """Ranking a records containing specified words and returns a sorted list.
    input:
    rank_method_code - the code of the method, from the name field in rnkMETHOD
    lwords - a list of words from the query
    hitset - a list of hits for the query found by search_engine
    rank_limit_relevance - show only records with a rank value above this
    verbose - verbose value
    output:
    reclist - a list of sorted records: [[23, 34], [344, 24], [1, 01]]
    prefix - what to show before the rank value
    postfix - what to show after the rank value
    voutput - contains extra information, content dependent on verbose value"""

    global voutput
    startCreate = time.time()

    if verbose > 0:
        voutput += "<br />Running rank method: %s, using word_frequency function in bibrank_record_sorter<br />" % rank_method_code
 
    lwords_old = lwords
    lwords = []
    #Check terms, remove non alphanumeric characters. Use both unstemmed and stemmed version of all terms.
    for i in range(0, len(lwords_old)):
        term = string.lower(lwords_old[i])
        if not methods[rank_method_code]["stopwords"] == "True" or methods[rank_method_code]["stopwords"] and not is_stopword(term, 1):
            lwords.append((term, methods[rank_method_code]["rnkWORD_table"]))
            terms = string.split(string.lower(re.sub(methods[rank_method_code]["chars_alphanumericseparators"], ' ', term)))
            for term in terms:
                if methods[rank_method_code].has_key("stemmer"): # stem word
                    term = stem(string.replace(term, ' ', ''), methods[rank_method_code]["stemmer"])
                if lwords_old[i] != term: #add if stemmed word is different than original word
                    lwords.append((term, methods[rank_method_code]["rnkWORD_table"]))

    (recdict, rec_termcount, lrecIDs_remove) = ({}, {}, {})
    #For each term, if accepted, get a list of the records using the term
    #calculate then relevance for each term before sorting the list of records
    for (term, table) in lwords:
        term_recs = run_sql("""SELECT term, hitlist FROM %s WHERE term=%%s""" % methods[rank_method_code]["rnkWORD_table"], (term, ))
        if term_recs: #if term exists in database, use for ranking
            term_recs = deserialize_via_marshal(term_recs[0][1])
            (recdict, rec_termcount) = calculate_record_relevance((term, int(term_recs["Gi"][1])) , term_recs, hitset, recdict, rec_termcount, verbose, quick=None)
            del term_recs

    rnkdict_ = rnkDict()
    rnkdict_.put(recdict)
    rnkdict  = rnkDict()
    rnkdict.loaddict("logistic")
    rnkmerge = rnkDict()
#    rnkmerge.merge(rnkdict_, rnkdict, [0, 0.1, 0.1])
    rnklist = []
    rnklist.append(rnkdict_.getdict())
    rnklist.append(rnkdict.getdict())
    rnkmerge.octopus(rnklist, [0.5, 0.00001, 0.9999])
    recdict = rnkmerge.getdict()

    if len(recdict) == 0 or (len(lwords) == 1 and lwords[0] == ""):
        return (None, "Records not ranked. The query is not detailed enough, or not enough records found, for ranking to be possible.", "", voutput)
    else: #sort if we got something to sort
        (reclist, hitset) = sort_record_relevance(recdict, hitset, rank_limit_relevance, _rescale=0, verbose=0)

    #Add any documents not ranked to the end of the list
    if hitset:
        lrecIDs = list(hitset)                       #using 2-3mb
        reclist = zip(lrecIDs, [0] * len(lrecIDs)) + reclist      #using 6mb

    if verbose > 0:
        voutput += "<br />Current number of recIDs: %s<br />" % (methods[rank_method_code]["col_size"])
        voutput += "Number of terms: %s<br />" % run_sql("SELECT count(id) FROM %s" % methods[rank_method_code]["rnkWORD_table"])[0][0]
        voutput += "Terms: %s<br />" % lwords
        voutput += "Prepare and pre calculate time: %s<br />" % (str(time.time() - startCreate))
        voutput += "Total time used: %s<br />" % (str(time.time() - startCreate))
        rank_method_stat(reclist, lwords)

    return (reclist, methods[rank_method_code]["prefix"], methods[rank_method_code]["postfix"], voutput)


def word_similarity(rank_method_code, lwords, hitset, rank_limit_relevance, verbose):
    """Ranking a records containing specified words and returns a sorted list.
    input:
    rank_method_code - the code of the method, from the name field in rnkMETHOD
    lwords - a list of words from the query
    hitset - a list of hits for the query found by search_engine
    rank_limit_relevance - show only records with a rank value above this
    verbose - verbose value
    output:
    reclist - a list of sorted records: [[23,34], [344,24], [1,01]]
    prefix - what to show before the rank value
    postfix - what to show after the rank value
    voutput - contains extra information, content dependent on verbose value"""

    global voutput
    startCreate = time.time()

    if verbose > 0:
        voutput += "<br />Running rank method: %s, using word_frequency function in bibrank_record_sorter<br />" % rank_method_code

    lwords_old = lwords
    lwords = []
    #Check terms, remove non alphanumeric characters. Use both unstemmed and stemmed version of all terms.
    for i in range(0, len(lwords_old)):
        term = string.lower(lwords_old[i])
        if not methods[rank_method_code]["stopwords"] == "True" or methods[rank_method_code]["stopwords"] and not is_stopword(term, 1):
            lwords.append((term, methods[rank_method_code]["rnkWORD_table"]))
            terms = string.split(string.lower(re.sub(methods[rank_method_code]["chars_alphanumericseparators"], ' ', term)))
            for term in terms:
                if methods[rank_method_code].has_key("stemmer"): # stem word
                    term = stem(string.replace(term, ' ', ''), methods[rank_method_code]["stemmer"])
                if lwords_old[i] != term: #add if stemmed word is different than original word
                    lwords.append((term, methods[rank_method_code]["rnkWORD_table"]))

    (recdict, rec_termcount, lrecIDs_remove) = ({}, {}, {})
    #For each term, if accepted, get a list of the records using the term
    #calculate then relevance for each term before sorting the list of records
    for (term, table) in lwords:
        term_recs = run_sql("""SELECT term, hitlist FROM %s WHERE term=%%s""" % methods[rank_method_code]["rnkWORD_table"], (term,))
        if term_recs: #if term exists in database, use for ranking
            term_recs = deserialize_via_marshal(term_recs[0][1])
            (recdict, rec_termcount) = calculate_record_relevance((term, int(term_recs["Gi"][1])) , term_recs, hitset, recdict, rec_termcount, verbose, quick=None)
            del term_recs

    if len(recdict) == 0 or (len(lwords) == 1 and lwords[0] == ""):
        return (None, "Records not ranked. The query is not detailed enough, or not enough records found, for ranking to be possible.", "", voutput)
    else: #sort if we got something to sort
        (reclist, hitset) = sort_record_relevance(recdict, hitset, rank_limit_relevance, _rescale=0, verbose=0)

    #Add any documents not ranked to the end of the list
    if hitset:
        lrecIDs = list(hitset)                       #using 2-3mb
        reclist = zip(lrecIDs, [0] * len(lrecIDs)) + reclist      #using 6mb

    if verbose > 0:
        voutput += "<br />Current number of recIDs: %s<br />" % (methods[rank_method_code]["col_size"])
        voutput += "Number of terms: %s<br />" % run_sql("SELECT count(id) FROM %s" % methods[rank_method_code]["rnkWORD_table"])[0][0]
        voutput += "Terms: %s<br />" % lwords
        voutput += "Prepare and pre calculate time: %s<br />" % (str(time.time() - startCreate))
        voutput += "Total time used: %s<br />" % (str(time.time() - startCreate))
        rank_method_stat(rank_method_code, reclist, lwords)

    return (reclist, methods[rank_method_code]["prefix"], methods[rank_method_code]["postfix"], voutput)

def calculate_record_relevance(term, invidx, hitset, recdict, rec_termcount, verbose, quick=None):
    """Calculating the relevance of the documents based on the input, calculates only one word
    term - (term, query term factor) the term and its importance in the overall search
    invidx - {recid: tf, Gi: norm value} The Gi value is used as a idf value
    hitset - a hitset with records that are allowed to be ranked
    recdict - contains currently ranked records, is returned with new values
    rec_termcount - {recid: count} the number of terms in this record that matches the query
    verbose - verbose value
    quick - if quick=yes only terms with a positive qtf is used, to limit the number of records to sort"""


    (t, qtf) = term
    if invidx.has_key("Gi"):#Gi = weigth for this term, created by bibrank_word_indexer
        Gi = invidx["Gi"][1]
        del invidx["Gi"]
    else: #if not existing, bibrank should be run with -R
        return (recdict, rec_termcount)

    if not quick or (qtf >= 0 or (qtf < 0 and len(recdict) == 0)):
        #Only accept records existing in the hitset received from the search engine
        for (j, tf) in invidx.iteritems():
            if j in hitset:#only include docs found by search_engine based on query
                try: #calculates rank value
                    recdict[j] = recdict.get(j, 0) + int(math.log(tf[0] * Gi * tf[1] * qtf))
                except:
                    return (recdict, rec_termcount)
                rec_termcount[j] = rec_termcount.get(j, 0) + 1 #number of terms from query in document
    elif quick: #much used term, do not include all records, only use already existing ones
        for (j, tf) in recdict.iteritems(): #i.e: if doc contains important term, also count unimportant
            if invidx.has_key(j):
                tf = invidx[j]
                recdict[j] = recdict.get(j, 0) + int(math.log(tf[0] * Gi * tf[1] * qtf))
                rec_termcount[j] = rec_termcount.get(j, 0) + 1 #number of terms from query in document

    return (recdict, rec_termcount)

def calculate_record_relevance_findsimilar(term, invidx, hitset, recdict, rec_termcount, verbose, quick=None):
    """Calculating the relevance of the documents based on the input, calculates only one word
    term - (term, query term factor) the term and its importance in the overall search
    invidx - {recid: tf, Gi: norm value} The Gi value is used as a idf value
    hitset - a hitset with records that are allowed to be ranked
    recdict - contains currently ranked records, is returned with new values
    rec_termcount - {recid: count} the number of terms in this record that matches the query
    verbose - verbose value
    quick - if quick=yes only terms with a positive qtf is used, to limit the number of records to sort"""


    (t, qtf) = term
    if invidx.has_key("Gi"): #Gi = weigth for this term, created by bibrank_word_indexer
        Gi = invidx["Gi"][1]
        del invidx["Gi"]
    else: #if not existing, bibrank should be run with -R
        return (recdict, rec_termcount)

    if not quick or (qtf >= 0 or (qtf < 0 and len(recdict) == 0)):
        #Only accept records existing in the hitset received from the search engine
        for (j, tf) in invidx.iteritems():
            if j in hitset: #only include docs found by search_engine based on query
                #calculate rank value
                recdict[j] = recdict.get(j, 0) + int((1 + math.log(tf[0])) * Gi * tf[1] * qtf)
                rec_termcount[j] = rec_termcount.get(j, 0) + 1 #number of terms from query in document
    elif quick: #much used term, do not include all records, only use already existing ones
        for (j, tf) in recdict.iteritems(): #i.e: if doc contains important term, also count unimportant
            if invidx.has_key(j):
                tf = invidx[j]
                recdict[j] = recdict[j] + int((1 + math.log(tf[0])) * Gi * tf[1] * qtf)
                rec_termcount[j] = rec_termcount.get(j, 0) + 1 #number of terms from query in document

    return (recdict, rec_termcount)

def sort_record_relevance(recdict, hitset, rank_limit_relevance, _rescale=0, verbose=0):
    """Sorts the dictionary and returns records with a relevance higher than the given value.
    recdict - {recid: value} unsorted
    rank_limit_relevance - a value > 0 usually
    verbose - verbose value"""

    startCreate = time.time()
    global voutput
    reclist = []

    #remove all ranked documents so that unranked can be added to the end
    hitset -= recdict.keys()

    #gives each record a score between 0-100
    divideby = max(recdict.values())
    for (j, w) in recdict.iteritems():
#        w = int(w * 100 / divideby)
#        if w >= rank_limit_relevance:
        reclist.append((j, w))
        
    if _rescale:
        reclist = rescale(recdict, rank_limit_relevance)

    #sort scores
    reclist.sort(lambda x, y: cmp(x[1], y[1]))

    if verbose > 0:
        voutput += "Number of records sorted: %s<br />" % len(reclist)
        voutput += "Sort time: %s<br />" % (str(time.time() - startCreate))
    return (reclist, hitset)

def sort_record_relevance_findsimilar(recdict, rec_termcount, hitset, rank_limit_relevance, verbose):
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
        w = int(w * 100 / divideby)
        if w >= rank_limit_relevance:
            reclist.append((j, w))

    #sort scores
    reclist.sort(lambda x, y: cmp(x[1], y[1]))

    if verbose > 0:
        voutput += "Number of records sorted: %s<br />" % len(reclist)
        voutput += "Sort time: %s<br />" % (str(time.time() - startCreate))
    return (reclist, hitset)

def rank_method_stat(rank_method_code, reclist, lwords):
    """Shows some statistics about the searchresult.
    rank_method_code - name field from rnkMETHOD
    reclist - a list of sorted and ranked records
    lwords - the words in the query"""

    global voutput
    if len(reclist) > 20:
        j = 20
    else:
        j = len(reclist)

    voutput += "<br />Rank statistics:<br />"
    for i in range(1, j + 1):
        voutput += "%s,Recid:%s,Score:%s<br />" % (i,reclist[len(reclist) - i][0],reclist[len(reclist) - i][1])
        for (term, table) in lwords:
            term_recs = run_sql("""SELECT hitlist FROM %s WHERE term=%%s""" % table, (term,))
            if term_recs:
                term_recs = deserialize_via_marshal(term_recs[0][0])
                if term_recs.has_key(reclist[len(reclist) - i][0]):
                    voutput += "%s-%s / " % (term, term_recs[reclist[len(reclist) - i][0]])
        voutput += "<br />"

    voutput += "<br />Score variation:<br />"
    count = {}
    for i in range(0, len(reclist)):
        count[reclist[i][1]] = count.get(reclist[i][1], 0) + 1
    i = 100
    while i >= 0:
        if count.has_key(i):
            voutput += "%s-%s<br />" % (i, count[i])
        i -= 1

try:
    import psyco
    psyco.bind(find_similar)
    psyco.bind(rank_by_method)
    psyco.bind(calculate_record_relevance)
    psyco.bind(word_similarity)
    psyco.bind(enhanced_ranking)
    psyco.bind(sort_record_relevance)
except StandardError, e:
    pass

