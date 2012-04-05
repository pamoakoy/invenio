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

"""This file defines class used for the management of DRANK ranking tables. 
It includes rank normalization, rescaling, filtering and aggregation/merging."""
import sys
from math import exp
from invenio.dbquery import deserialize_via_marshal, serialize_via_marshal, run_sql
from invenio.bibrank_downloads_indexer import uniq
from scipy.stats.mstats import scoreatpercentile
from scipy.stats.kde import gaussian_kde
from numpy import inf
from Numeric import array 
  
class rnkDict:
    """Class to manage ranking tables."""

    def __init__(self, rnkdict={}):
        self.content = rnkdict

    def strip(self):
        """Remove entries that contain nil scores. Inverse of fill."""
        _content = {}
        keys     = filter(lambda x: self.content[x] > 0, self.content.keys())
        values   = map(lambda x: self.content[x], keys)
        content  = zip(keys, values)
        for item in content:
            _content[item[0]] = item[1] 
        self.content = _content
        return

    def fill(self):
        """Fill ranking table so that it contains all keys including those with nil scores. Inverse of strip."""
        self.sanitize()
        return

    def sanitize(self, nil=0):
        """Sanitize the ranking table and make sure all keys are present."""
        keys = self.content.keys()
        for item in filter(lambda x: x not in keys, range(max(keys))):
            self.content[item] = nil
        return

    def put(self, dictionary):
        """Insert content given in the dictionary into the table. Inverse of get."""
        self.content = dictionary
        return
  
    def getdict(self):
        """Get the content of the table, return a dictionary. Inverse of put."""
        return self.content

    def read_from_file(self, filename, separator="|"):
        """Read ranking table from a CSV file."""
        filehandle = open(filename, "r")
        for line in filehandle.readlines():
            (key, value) = line.split(separator)
            try:
                self.content[int(key)] = int(value)
            except ValueError as err:
                sys.stderr.write("Warning: %s\n" % err)
        return
    
    def lookup_in_lut(self,lut,reclist):
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
                        looked_up_scores[int(key)]= self.get_interpolated_value(score, freq_values[k], freq_values[k+1], lut[freq_values[k]], lut[freq_values[k+1]])
                        break
    #    print "record_pscore:{0}".format(record_pscore)
        self.content=looked_up_scores
        return
    
    def get_interpolated_value(self,x, x1, x2, y1, y2):
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
    
    def loadlut(self, dictname):
        """Load the ranking table into memory. Inverse of savedict."""
#        identifier = run_sql("SELECT id from rnkMETHOD where name=\"%s\"" % dictname)
        res = run_sql("SELECT drank_lut FROM rnkDRANKLUT WHERE drank_name=\"%s\"" %dictname)
        
        if res:
            self.content = deserialize_via_marshal(res[0][0])
        else:
            self.content = {}
        return
    
    def loaddict(self, dictname):
        """Load the ranking table into memory. Inverse of savedict."""
        identifier = run_sql("SELECT id from rnkMETHOD where name=\"%s\"" % dictname)
        res = run_sql("SELECT relevance_data FROM rnkMETHODDATA WHERE id_rnkMETHOD=\"%s\"" %  identifier[0][0])
        if res:
            self.content = deserialize_via_marshal(res[0][0])
        else:
            self.content = {}
        return
       
    def savedict(self, dictname):
        """Save ranking table into the database (rewrite). Inverse of loaddict."""
        mid = run_sql("SELECT id from rnkMETHOD where name=\"%s\"" % dictname)
        self.delete_from_db(dictname)
        serdata = serialize_via_marshal(self.content)
        midstr  = str(mid[0][0])
        run_sql("INSERT INTO rnkMETHODDATA(id_rnkMETHOD, relevance_data) VALUES (%s,%s)", (midstr, serdata,))
#        run_sql("UPDATE rnkMETHOD SET last_updated=%s WHERE name=%s", (date, rank_method_code))
        return

    def write(self, mode=1):
        """Print the ranking table to the standard output."""
        if mode:
            for key in self.content.keys():
                sys.stdout.write("%s : %s\n" % (key, self.content[key]))
        return    

    def length(self):
        """Get length of the ranking table."""
        return len(self.content.keys())

    def rescale(self):
        """Rescale the ranking table so that scores appear on scale in [0, 1]."""
        _list = self.content.values()
        xmin = min(_list)
        xmax = max(_list)
        span = xmax - xmin
        return map(lambda x: (x-xmin)/span, _list)
  
    def percrescale(self):
        """Rescale the ranking table to integers so that scores appear to be on scale in [0, 100]."""
        _content = {}
        xmin = min(self.content.values())
        xmax = max(self.content.values())
        span = xmax - xmin
        res  = map(lambda y, x: (y, int(100*(x-xmin)/span)), self.content.keys(), self.content.values())
        for item in res:
            _content[item[0]] = item[1]
        self.content = _content
        return

    def rank(self):
        """Rank. Return list of ranked key/value pairs."""
        return map(lambda x:list(x), sorted(self.content.items(), key=lambda (k, v):(v, k)))

    def octopus(self, ldict, lweight):
        """Octopus merge. Merge two and more ranking tables in one go."""
        allkeys  = []
        _allkeys = []
        _c_dict  = {}
        for i in range(len(ldict)-1):
            for key in ldict[i].keys():
                allkeys.append(key)
        _allkeys = uniq(allkeys)
        counter = 0
        for key in _allkeys:
            counter += 1
            cumulated = 0
            for i in range(len(ldict)-1):
                dico = ldict[i]
                if dico.has_key(key):
                    cumulated += float(dico[key]) * lweight[i+1]
            _c_dict[key] = 1/(1+exp(1)**(-1*cumulated))
        self.content = _c_dict
        return

    def merge(self, a_dict, b_dict, weight=[0, 0, 0]):
        """Merge two ranking tables."""
#        c_list = map(lambda x,y: 1/(1+exp(1)**(-1*(w[0]+w[1]*x+w[2]*y))), a_list, b_list)
        _a_dict = a_dict.getdict()
        _b_dict = b_dict.getdict()
        _c_dict = {}
  
        for item in _a_dict.keys():
            if _b_dict.has_key(item):
                _c_dict[item] = 1/(1+exp(1)**(-1*(weight[0]+weight[1]*_a_dict[item]+weight[2]*_b_dict[item])))
                del _b_dict[item]
            else:
                _c_dict[item] = 1/(1+exp(1)**(-1*(weight[0]+weight[1]*_a_dict[item])))
        for item in _b_dict.keys():
            _c_dict[item] = 1/(1+exp(1)**(-1*(weight[0]+weight[2]*_b_dict[item])))
        self.content = _c_dict
        return 
    
    
    def multiply(self, a_dict, b_dict):
        """Create a new table multiplying values from two other ranking tables.""" 
        _a_dict = a_dict.getdict()
        _b_dict = b_dict.getdict()
        _c_dict = {}

        for item in _a_dict.keys():
            if _b_dict.has_key(item):
                _c_dict[item] = _a_dict[item] * _b_dict[item]
                del _b_dict[item]
            else:
                _c_dict[item] = 0
        for item in _b_dict.keys():
            _c_dict[item] = 0
        self.content = _c_dict
        return        

    def normalize(self, outlier=1):
        """Normalize scores in the ranking table w/kde. Remove outliers first."""
        scorelist = self.content.values()
        if outlier:
            outlier_above = scoreatpercentile(scorelist, 100-outlier)
            outlier_below = scoreatpercentile(scorelist, outlier)
            scorelist_new = []
            for item in scorelist:
                if float(item) >= float(outlier_below) and float(item) <= float(outlier_above):
                    scorelist_new.append(item)
            scorelist = scorelist_new
        density_estimate = gaussian_kde(array(scorelist))
        conversiontable  = {}
        for item in self.content.items():
            self.content[item[0]] = density_estimate.integrate_box_1d(-inf, item[1])
        return

    def normalize_take_displays_into_account(self, outlier=1):
        """Normalize scores in the ranking table w/kde. Remove outliers first. Take displays into account."""
        _dico = rnkDict()
        _dico.loaddict("displays")
        mydico = _dico.getdict()
        scorelist = []
        for key in self.content.keys():
            if mydico.has_key(key):
                for iter in range(int(mydico[key])):
                    scorelist.append(self.content[key])
            else:
                scorelist.append(self.content[key])
        if outlier:
            outlier_above = scoreatpercentile(scorelist, 100-outlier)
            outlier_below = scoreatpercentile(scorelist, outlier)
            scorelist_new = []
            for item in scorelist:
                if float(item) >= float(outlier_below) and float(item) <= float(outlier_above):
                    scorelist_new.append(item)
            scorelist = scorelist_new
        density_estimate = gaussian_kde(array(scorelist))
        conversiontable  = {}
        for item in self.content.items():
            self.content[item[0]] = density_estimate.integrate_box_1d(-inf, item[1])
        return

    def _remove_ties(self):
        """Arbitrary tie removal."""
        temp = {}
        for item in self.content.values():
            temp[item] = item
        return

    def delete_from_db(self, dictname):
        """Delete the ranking table from database."""
        identifier = run_sql("SELECT id from rnkMETHOD where name=\"%s\"" % dictname)
        run_sql("DELETE FROM rnkMETHODDATA WHERE id_rnkMETHOD=\"%s\"" % identifier[0][0])
        return
 
    def gethirsch(self):
        """Get hirsch index of the ranking table.""" 
        hirsch  = 0
        counter = 0
        ranked  = self.rank()
        ranked.reverse()
        for item in ranked:
            counter += 1
            if item[1] >= counter:
                hirsch += 1
        return hirsch

    def filter(self, keylist):
        """Filter the ranking table with a given list of keys."""
        _content = {}
        for key in keylist:
            if self.content.has_key(key):
                _content[key] = self.content[key]
        self.content = _content
        return

    def clean(self):
        """Clean the ranking table."""
#        """Remove empty spaces from strings and make integer."""
        _content = {}
        for key in self.content.keys():
            _value   = int(self.content[key].lstrip().rstrip())
            _key     = int(key.lstrip().rstrip())
            _content[_key] = _value
        self.content = _content
        return
