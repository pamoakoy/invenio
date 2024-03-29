# -*- coding: utf-8 -*-

## This file is part of Invenio.
## Copyright (C) 2008, 2010, 2011 CERN.
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

# pylint: disable=C0301

"""Invenio Search Engine query parsers."""

import re
import string
from datetime import datetime

try:
    import dateutil.parser
    GOT_DATEUTIL = True
except ImportError:
    # Ok, no date parsing is possible, but continue anyway,
    # since this package is only recommended, not mandatory.
    GOT_DATEUTIL = False

from invenio.bibindex_engine_tokenizer import BibIndexFuzzyNameTokenizer as FNT
from invenio.logicutils import to_cnf


NameScanner = FNT()


class InvenioWebSearchMismatchedParensError(Exception):
    """Exception for parse errors caused by mismatched parentheses."""
    def __init__(self, message):
        """Initialization."""
        self.message = message
    def __str__(self):
        """String representation."""
        return repr(self.message)


class SearchQueryParenthesisedParser(object):
    """Search query parser that handles arbitrarily-nested parentheses

    Parameters:
    * substitution_dict: a dictionary mapping strings to other strings.  By
      default, maps 'and', 'or' and 'not' to '+', '|', and '-'.  Dictionary
      values will be treated as valid operators for output.
    """

    def __init__(self, substitution_dict = {'and': '+', 'or': '|', 'not': '-'}):
        self.substitution_dict = substitution_dict
        self.specials = set(['(', ')', '+', '|', '-', '+ -'])
        self.__tl_idx = 0
        self.__tl_len = 0

    # I think my names are both concise and clear
    # pylint: disable=C0103
    def _invenio_to_python_logical(self, q):
        """Translate the + and - in invenio query strings into & and ~."""
        p = q
        p = re.sub('\+ -', '&~', p)
        p = re.sub('\+', '&', p)
        p = re.sub('-', '~', p)
        p = re.sub(' ~', ' & ~', p)
        return p

    def _python_logical_to_invenio(self, q):
        """Translate the & and ~ in logical expression strings into + and -."""
        p = q
        p = re.sub('\& ~', '-', p)
        p = re.sub('~', '-', p)
        p = re.sub('\&', '+', p)
        return p
    # pylint: enable=C0103

    def parse_query(self, query):
        """Make query into something suitable for search_engine.

        This is the main entry point of the class.

        Given an expression of the form:
        "expr1 or expr2 (expr3 not (expr4 or expr5))"
        produces annoted list output suitable for consumption by search_engine,
        of the form:
        ['+', 'expr1', '|', 'expr2', '+', 'expr3 - expr4 | expr5']

        parse_query() is a wrapper for self.tokenize() and self.parse().
        """
        toklist = self.tokenize(query)
        depth, balanced, dummy_d0_p = self.nesting_depth_and_balance(toklist)
        if not balanced:
            raise SyntaxError("Mismatched parentheses in "+str(toklist))
        toklist, var_subs = self.substitute_variables(toklist)
        if depth > 1:
            toklist = self.tokenize(self.logically_reduce(toklist))
        return self.parse(toklist, var_subs)

    def substitute_variables(self, toklist):
        """Given a token list, return a copy of token list in which all free
        variables are bound with boolean variable names of the form 'pN'.
        Additionally, all the substitutable logical operators are exchanged
        for their symbolic form and implicit ands are made explicit

        e.g., ((author:'ellis, j' and title:quark) or author:stevens jones)
        becomes:
              ((p0 + p1) | p2 + p3)
        with the substitution table:
        {'p0': "author:'ellis, j'", 'p1': "title:quark",
         'p2': "author:stevens", 'p3': "jones" }

        Return value is the substituted token list and a copy of the
        substitution table.
        """
        def labels():
            i = 0
            while True:
                yield 'p'+str(i)
                i += 1

        def filter_front_ands(toklist):
            """Filter out extra logical connectives and whitespace from the front."""
            while toklist[0] == '+' or toklist[0] == '|' or toklist[0] == '':
                toklist = toklist[1:]
            return toklist

        var_subs = {}
        labeler = labels()
        new_toklist = ['']
        cannot_be_anded = self.specials.difference((')',))
        for token in toklist:
            token = token.lower()
            if token in self.substitution_dict:
                if token == 'not' and new_toklist[-1] == '+':
                    new_toklist[-1] = '-'
                else:
                    new_toklist.append(self.substitution_dict[token])
            elif token == '(':
                if new_toklist[-1] not in self.specials:
                    new_toklist.append('+')
                new_toklist.append(token)
            elif token not in self.specials:
                # apparently generators are hard for pylint to figure out
                # Turns off msg about labeler not having a 'next' method
                # pylint: disable=E1101
                label = labeler.next()
                # pylint: enable=E1101
                var_subs[label] = token
                if new_toklist[-1] not in cannot_be_anded:
                    new_toklist.append('+')
                new_toklist.append(label)
            else:
                if token == '-' and new_toklist[-1] == '+':
                    new_toklist[-1] = '-'
                else:
                    new_toklist.append(token)
        return filter_front_ands(new_toklist), var_subs

    def nesting_depth_and_balance(self, token_list):
        """Checks that parentheses are balanced and counts how deep they nest"""
        depth = 0
        maxdepth = 0
        depth0_pairs = 0
        good_depth = True
        for i in range(len(token_list)):
            token = token_list[i]
            if token == '(':
                if depth == 0:
                    depth0_pairs += 1
                depth += 1
                if depth > maxdepth:
                    maxdepth += 1
            elif token == ')':
                depth -= 1
            if depth == -1:        # can only happen with unmatched )
                good_depth = False # so force depth check to fail
                depth = 0          # but keep maxdepth in good range
        return maxdepth, depth == 0 and good_depth, depth0_pairs

    def logically_reduce(self, token_list):
        """Return token_list in conjunctive normal form as a string.

        CNF has the property that there will only ever be one level of
        parenthetical nesting, and all distributable operators (such as
        the not in -(p | q) will be fully distributed (as -p + -q).
        """

        maxdepth, dummy_balanced, d0_p = self.nesting_depth_and_balance(token_list)
        s = ' '.join(token_list)
        s = self._invenio_to_python_logical(s)
        last_maxdepth = 0
        while maxdepth != last_maxdepth:             # XXX: sometimes NaryExpr doesn't
            try:                                     # fully flatten Expr; but it usually
                s = str(to_cnf(s))                   # does in 2 passes FIXME: diagnose
            except SyntaxError:
                raise SyntaxError(str(s)+" couldn't be converted to a logic expression.")
            last_maxdepth = maxdepth
            maxdepth, dummy_balanced, d0_p = self.nesting_depth_and_balance(self.tokenize(s))
        if d0_p == 1 and s[0] == '(' and s[-1] == ')': # s can come back with extra parens
            s = s[1:-1]
        s = self._python_logical_to_invenio(s)
        return s

    def tokenize(self, query):
        """Given a query string, return a list of tokens from that string.

        * Isolates meaningful punctuation: ( ) + | -
        * Keeps single- and doubl-quoted strings together without interpretation.
        * Splits everything else on whitespace.

        i.e.:
        "expr1|expr2 (expr3-(expr4 or expr5))"
        becomes:
        ['expr1', '|', 'expr2', '(', 'expr3', '-', '(', 'expr4', 'or', 'expr5', ')', ')']
        """
        ###
        # Invariants:
        # * Query is never modified
        # * In every loop iteration, querytokens grows to the right
        # * The only return point is at the bottom of the function, and the only
        #   return value is querytokens
        ###

        def get_tokens(s):
            """
            Given string s, return a list of s's tokens.

            Adds space around special punctuation, then splits on whitespace.
            """
            s = ' '+s
            s = s.replace('->', '####DATE###RANGE##OP#') # XXX: Save '->'
            s = re.sub('(?P<outside>[a-zA-Z0-9_,]+)\((?P<inside>[a-zA-Z0-9_,]*)\)',
                       '#####\g<outside>####PAREN###\g<inside>##PAREN#', s) # XXX: Save U(1) and SL(2,Z)
            for char in self.specials:
                if char == '-':
                    s = s.replace(' -', ' '+char+' ')
                    s = s.replace(')-', ') '+char+' ')
                    s = s.replace('-(', ' '+char+' (')
                else:
                    s = s.replace(char, ' '+char+' ')
            s = re.sub('#####(?P<outside>[a-zA-Z0-9_,]+)####PAREN###(?P<inside>[a-zA-Z0-9_,]*)##PAREN#',
                       '\g<outside>(\g<inside>)', s) # XXX: Restore U(1) and SL(2,Z)
            s = s.replace('####DATE###RANGE##OP#', '->') # XXX: Restore '->'
            return s.split()

        querytokens = []
        current_position = 0
        re_quotepairs = re.compile(r'[^\\](".*?[^\\]")|[^\\](\'.*?[^\\]\')')

        for match in re_quotepairs.finditer(query):
            # special case for when our regexp captures a single char before
            # the quotes.  This is faster (in development time and code)
            # than a much more complicated and complete regexp or using an
            # FSM for quote balancing.  XXX: But is there a better way?
            match_start = match.start()
            quoted_region = match.group(0).strip()
            if quoted_region[0] not in "'\"":
                match_start += 1            # capture char 0 in unquoted below
                quoted_region = quoted_region[1:]

            # clean the content after the previous quotes and before current quotes
            unquoted = query[current_position : match_start]
            querytokens.extend(get_tokens(unquoted))

            # XXX: sometimes we end up with labeled quoted regions split off from their
            # labels, e.g., 'title:"compton scattering"' becomes
            # ['title:', '"compton scattering"'] rather than ['title:"compton scattering"']
            # This corrects for that.
            if querytokens[-1][-1] == ':':
                querytokens[-1] += quoted_region
            else:
                # add our newly tokenized content to the token list
                querytokens.extend([quoted_region])

            # move current position to the end of the tokenized content
            current_position = match.end()

        # get tokens from the last appearance of quotes until the query end
        unquoted = query[current_position : len(query)]
        querytokens.extend(get_tokens(unquoted))

        return querytokens

    def parse(self, token_list, variable_substitution_dict=None):
        """Make token_list consumable by search_engine.

        Turns a list of tokens and a variable mapping into a grouped list
        of subexpressions in the format suitable for use by search_engine,
        e.g.:
        ['+', 'searchterm', '-', 'searchterm to exclude', '|', 'another term']

        Incidentally, this works recursively so parens can cause arbitrarily
        deep nestings.  But since the search_engine doesn't know about nested
        structures, we need to flatten the input structure first.
        """
        ###
        # Invariants:
        # * Token list is never modified
        # * Balanced parens remain balanced; unbalanced parens are an error
        # * Individual tokens may only be exchanged for items in the variable
        #   substitution dict; otherwise they pass through unmolested
        # * Return value is built up mostly as a stack
        ###

        op_symbols = self.substitution_dict.values()
        self.__tl_idx = 0
        self.__tl_len = len(token_list)

        def inner_parse(token_list, open_parens=False):

            if open_parens:
                parsed_values = []
            else:
                parsed_values = ['+']

            i = 0
            while i < len(token_list):
                token = token_list[i]
                if i > 0 and parsed_values[-1] not in op_symbols:
                    parsed_values.append('+')
                if token == '(':
                    offset = self.__tl_len - len(token_list)
                    inner_value = inner_parse(token_list[i+1:], True)
                    inner_value = ' '.join(inner_value)
                    parsed_values.append(inner_value)
                    self.__tl_idx += 1
                    i = self.__tl_idx - offset
                elif token == ')':
                    if parsed_values[-1] in op_symbols:
                        parsed_values = parsed_values[:-1]
                    if parsed_values[0] == '+' and parsed_values[1] in op_symbols:
                        parsed_values = parsed_values[1:]
                    return parsed_values
                elif token in op_symbols:
                    if len(parsed_values) > 0:
                        parsed_values[-1] = token
                    else:
                        parsed_values = [token]
                else:
                    if variable_substitution_dict != None and token in variable_substitution_dict:
                        token = variable_substitution_dict[token]
                    parsed_values.append(token)
                i += 1
                self.__tl_idx += 1

            # If we have an extra start symbol, remove the default one
            if parsed_values[1] in op_symbols:
                parsed_values = parsed_values[1:]
            return parsed_values

        return inner_parse(token_list, False)


class SpiresToInvenioSyntaxConverter:
    """Converts queries defined with SPIRES search syntax into queries
    that use Invenio search syntax.
    """

    # Constants defining fields
    _DATE_ADDED_FIELD = 'datecreated:'
    _DATE_UPDATED_FIELD = 'datemodified:'
    _DATE_FIELD = 'year:'

    _A_TAG = 'author:'
    _EA_TAG = 'exactauthor:'


    # Dictionary containing the matches between SPIRES keywords
    # and their corresponding Invenio keywords or fields
    # SPIRES keyword : Invenio keyword or field
    _SPIRES_TO_INVENIO_KEYWORDS_MATCHINGS = {
        # affiliation
        'affiliation' : 'affiliation:',
        'affil' : 'affiliation:',
        'aff' : 'affiliation:',
        'af' : 'affiliation:',
        'institution' : 'affiliation:',
        'inst' : 'affiliation:',
        # any field
        'any' : 'anyfield:',
        # bulletin
        'bb' : 'reportnumber:',
        'bbn' : 'reportnumber:',
        'bull' : 'reportnumber:',
        'bulletin-bd' : 'reportnumber:',
        'bulletin-bd-no' : 'reportnumber:',
        'eprint' : 'reportnumber:',
        # citation / reference
        'c' : 'reference:',
        'citation' : 'reference:',
        'cited' : 'reference:',
        'jour-vol-page' : 'reference:',
        'jvp' : 'reference:',
        # collaboration
        'collaboration' : 'collaboration:',
        'collab-name' : 'collaboration:',
        'cn' : 'collaboration:',
        # conference number
        'conf-number' : '111__g:',
        'cnum' : '111__g:',
        # country
        'cc' : '044__a:',
        'country' : '044__a:',
        # date
        'date': _DATE_FIELD,
        'd': _DATE_FIELD,
        # date added
        'date-added': _DATE_ADDED_FIELD,
        'dadd': _DATE_ADDED_FIELD,
        'da': _DATE_ADDED_FIELD,
        # date updated
        'date-updated': _DATE_UPDATED_FIELD,
        'dupd': _DATE_UPDATED_FIELD,
        'du': _DATE_UPDATED_FIELD,
        # first author
        'fa' : '100__a:',
        'first-author' : '100__a:',
        # author
        'a' : 'author:',
        'au' : 'author:',
        'author' : 'author:',
        'name' : 'author:',
        # exact author
        # this is not a real keyword match. It is pseudo keyword that
        # will be replaced later with author search
        'ea' : 'exactauthor:',
        'exact-author' : 'exactauthor:',
        # experiment
        'exp' : 'experiment:',
        'experiment' : 'experiment:',
        'expno' : 'experiment:',
        'sd' : 'experiment:',
        'se' : 'experiment:',
        # journal
        'journal' : 'journal:',
        'j' : 'journal:',
        'published_in' : 'journal:',
        'spicite' : 'journal:',
        'vol' : 'journal:',
        # journal page
        'journal-page' : '773__c:',
        'jp' : '773__c:',
        # journal year
        'journal-year' : '773__y:',
        'jy' : '773__y:',
        # key
        'key' : '970__a:',
        'irn' : '970__a:',
        'record' : '970__a:',
        'document' : '970__a:',
        'documents' : '970__a:',
        # keywords
        'k' : 'keyword:',
        'keywords' : 'keyword:',
        'kw' : 'keyword:',
        # note
        'note' : '500__a:',
        'n' : '500__a:',
        # old title
        'old-title' : '246__a:',
        'old-t' : '246__a:',
        'ex-ti' : '246__a:',
        'et' : '246__a:',
        # ppf subject
        'ppf-subject' : '650__a:',
        'ps' : '650__a:',
        'scl' : '650__a:',
        'status' : '650__a:',
        # report number
        'r' : 'reportnumber:',
        'rn' : 'reportnumber:',
        'rept' : 'reportnumber:',
        'report' : 'reportnumber:',
        'report-num' : 'reportnumber:',
        # title
        't' : 'title:',
        'ti' : 'title:',
        'title' : 'title:',
        'with-language' : 'title:',
        # fulltext
        'fulltext' : 'fulltext:',
        'ft' : 'fulltext:',
        # topic
        'topic' : '695__a:',
        'tp' : '695__a:',
        'hep-topic' : '695__a:',
        'desy-keyword' : '695__a:',
        'dk' : '695__a:',

        # second-order operators
        'refersto' : 'refersto:',
        'refs': 'refersto:',
        'citedby' : 'citedby:',

        # replace all the keywords without match with empty string
        # this will remove the noise from the unknown keywrds in the search
        # and will in all fields for the words following the keywords


        # category
        'arx' : '037__c:',
        'category' : '037__c:',
        # primarch
        'parx' : '037__c:',
        'primarch' : '037__c:',
        # texkey
        'texkey' : '035__z:',
        # type code
        'tc' : '690C_a:',
        'ty' : '690C_a:',
        'type' : '690C_a:',
        'type-code' : '690C_a',
        # field code
        'f' : '65017a:',
        'fc' : '65017a:',
        'field' : '65017a:',
        'field-code' : '65017a:',
        # coden
        'bc' : '',
        'browse-only-indx' : '',
        'coden' : '',
        'journal-coden' : '',
        # energy
        'e' : '',
        'energy' : '',
        'energyrange-code' : '',
        # exact experiment number
        'ee' : '',
        'exact-exp' : '',
        'exact-expno' : '',
        # hidden note
        'hidden-note' : '',
        'hn' : '',
        # ppf
        'ppf' : '',
        'ppflist' : '',
        # slac topics
        'ppfa' : '',
        'slac-topics' : '',
        'special-topics' : '',
        'stp' : '',
        # test index
        'test' : '',
        'testindex' : '',
        }

    def __init__(self):
        """Initialize the state of the converter"""
        self._months = {}
        self._month_name_to_month_number = {}
        self._init_months()
        self._compile_regular_expressions()

    def _compile_regular_expressions(self):
        """Compiles some of the regular expressions that are used in the class
        for higher performance."""

        # regular expression that matches the contents in single and double quotes
        # taking in mind if they are escaped.
        self._re_quotes_match = re.compile('[^\\\\](".*?[^\\\\]")|[^\\\\](\'.*?[^\\\\]\')')

        # for matching cases where kw needs distributing
        self._re_distribute_keywords = re.compile(r'\b(?P<keyword>\S*:)(?P<content>.+?)\s*(?P<combination>and not | and | or | not )\s*(?P<last_content>[^:]*?)(?= and not | and | or | not |$)', re.IGNORECASE)

        # regular expression that matches author patterns
        self._re_author_match = re.compile(r'\bauthor:\s*(?P<name>.+?)\s*(?= and not | and | or | not |$)', re.IGNORECASE)

        # regular expression that matches exact author patterns
        # the group defined in this regular expression is used in method
        # _convert_spires_exact_author_search_to_invenio_author_search(...)
        # in case of changes correct also the code in this method
        self._re_exact_author_match = re.compile(r'\bexactauthor:(?P<author_name>[^\'\"].*?[^\'\"]\b)(?= and not | and | or | not |$)', re.IGNORECASE)

        # regular expression that matches search term, its content (words that
        # are searched) and the operator preceding the term. In case that the
        # names of the groups defined in the expression are changed, the
        # chagned should be reflected in the code that use it.
        self._re_search_term_pattern_match = re.compile(r'\b(?P<combine_operator>find|and|or|not)\s+(?P<search_term>title:|keyword:|fulltext:)(?P<search_content>.+?)(?= and not | and | or | not |$)', re.IGNORECASE)

        # regular expression matching date after pattern
        self._re_date_after_match = re.compile(r'\b(?P<searchop>d|date|dupd|dadd|da|date-added|du|date-updated)\b\s*(after|>)\s*(?P<search_content>.+?)(?= and not | and | or | not |$)', re.IGNORECASE)

        # regular expression matching date after pattern
        self._re_date_before_match = re.compile(r'\b(?P<searchop>d|date|dupd|dadd|da|date-added|du|date-updated)\b\s*(before|<)\s*(?P<search_content>.+?)(?= and not | and | or | not |$)', re.IGNORECASE)

        # regular expression that matches date searches which have been
        # keyword-substituted
        #self._re_keysubbed_date_expr = re.compile(r'\b(?P<term>(' + self._DATE_ADDED_FIELD + ')|(' + self._DATE_UPDATED_FIELD + ')|(' + self._DATE_FIELD + '))\s*(?P<content>.+)(?= and not | and | or | not |$)', re.IGNORECASE)
        self._re_keysubbed_date_expr = re.compile(r'\b(?P<term>(' + self._DATE_ADDED_FIELD + ')|(' + self._DATE_UPDATED_FIELD + ')|(' + self._DATE_FIELD + '))(?P<content>.+?)(?= and not | and | or | not |$)', re.IGNORECASE)

        # for finding (and changing) a variety of different SPIRES search keywords
        self._re_spires_find_keyword = re.compile('^(?P<find>f|fin|find)\s+(?P<query>.*)$', re.IGNORECASE)

        # patterns for subbing out spaces within quotes temporarily
        self._re_pattern_single_quotes = re.compile("'(.*?)'")
        self._re_pattern_double_quotes = re.compile("\"(.*?)\"")
        self._re_pattern_regexp_quotes = re.compile("\/(.*?)\/")
        self._re_pattern_space = re.compile("__SPACE__")
        self._re_pattern_IRN_search = re.compile(r'970__a:(?P<irn>\d+)')

    def is_applicable(self, query):
        """Is this converter applicable to this query?

        Returns True IFF the query starts with SPIRES' 'find' keyword or some
        acceptable variation thereof."""
        if self._re_spires_find_keyword.match(query.lower()):
            return True
        else:
            return False

    def convert_query(self, query):
        """Convert SPIRES syntax queries to Invenio syntax.

        Do nothing to queries not in SPIRES syntax."""

        # SPIRES syntax allows searches with 'find' or 'fin'.
        if self.is_applicable(query):

            # Everywhere else make the assumption that all and only queries
            # starting with 'find' are SPIRES queries.  Turn fin into find.
            query = self._re_spires_find_keyword.sub(lambda m: 'find '+m.group('query'), query)

            # these calls are before keywords replacement becuase when keywords
            # are replaced, date keyword is replaced by specific field search
            # and the DATE keyword is not match in DATE BEFORE or DATE AFTER
            query = self._convert_spires_date_before_to_invenio_span_query(query)
            query = self._convert_spires_date_after_to_invenio_span_query(query)

            # call to _replace_spires_keywords_with_invenio_keywords should be at the
            # beginning because the next methods use the result of the replacement
            query = self._replace_spires_keywords_with_invenio_keywords(query)
            query = self._distribute_keywords_across_combinations(query)

            query = self._convert_dates(query)
            query = self._convert_irns_to_spires_irns(query)
            query = self._convert_spires_author_search_to_invenio_author_search(query)
            query = self._convert_spires_exact_author_search_to_invenio_author_search(query)
            query = self._convert_spires_truncation_to_invenio_truncation(query)
            query = self._expand_search_patterns(query)

            # remove FIND in the beginning of the query as it is not necessary in Invenio
            query = query[4:]
            query = query.strip()

        return query

    def _init_months(self):
        """Defines a dictionary matching the name
        of the month with its corresponding number"""

        # this dictionary is used when generating match patterns for months
        self._months = {'jan':'01', 'january':'01',
                         'feb':'02', 'february':'02',
                         'mar':'03', 'march':'03',
                         'apr':'04', 'april':'04',
                         'may':'05', 'may':'05',
                         'jun':'06', 'june':'06',
                         'jul':'07', 'july':'07',
                         'aug':'08', 'august':'08',
                         'sep':'09', 'september':'09',
                         'oct':'10', 'october':'10',
                         'nov':'11', 'november':'11',
                         'dec':'12', 'december':'12'}
        # this dictionary is used to transform name of the month
        # to a number used in the date format. By this reason it
        # contains also the numbers itself to simplify the conversion
        self._month_name_to_month_number = {'1':'01', '01':'01',
                                            '2':'02', '02':'02',
                                            '3':'03', '03':'03',
                                            '4':'04', '04':'04',
                                            '5':'05', '05':'05',
                                            '6':'06', '06':'06',
                                            '7':'07', '07':'07',
                                            '8':'08', '08':'08',
                                            '9':'09', '09':'09',
                                            '10':'10',
                                            '11':'11',
                                            '12':'12',}
        # combine it with months in order to cover all the cases
        self._month_name_to_month_number.update(self._months)

    def _get_month_names_match(self):
        """Retruns part of a patter that matches month in a date"""

        months_match = ''
        for month_name in self._months.keys():
            months_match = months_match + month_name + '|'

        months_match = r'\b(' + months_match[0:-1] + r')\b'

        return months_match

    def _convert_dates(self, query):
        """Tries to find dates in query and make them look like ISO-8601."""

        def mangle_with_dateutils(query):
            DEFAULT = datetime(datetime.today().year, 1, 1)
            result = ''
            position = 0
            for match in self._re_keysubbed_date_expr.finditer(query):
                result += query[position : match.start()]

                isodates = []
                dates = match.group('content').split('->') # Warning: generalizing but should only ever be 2 items
                for datestamp in dates:
                    if datestamp != None:
                        if re.match('[0-9]{1,4}$', datestamp):
                            isodates.append(datestamp)
                        else:
                            try:
                                dtobj = dateutil.parser.parse(datestamp, default=DEFAULT)
                                if dtobj.day == 1:
                                    isodates.append("%d-%02d" % (dtobj.year, dtobj.month))
                                else:
                                    isodates.append("%d-%02d-%02d" % (dtobj.year, dtobj.month, dtobj.day))
                            except ValueError:
                                isodates.append(datestamp)

                daterange = '->'.join(isodates)
                result += match.group('term') + daterange
                position = match.end()
            result += query[position : ]
            return result

        if GOT_DATEUTIL:
            query = mangle_with_dateutils(query)
        # else do nothing with the dates
        return query

    def _convert_irns_to_spires_irns(self, query):
        """Prefix IRN numbers with SPIRES- so they match the INSPIRE format."""
        def create_replacement_pattern(match):
            """method used for replacement with regular expression"""
            return '970__a:SPIRES-' + match.group('irn')
        query = self._re_pattern_IRN_search.sub(create_replacement_pattern, query)
        return query

    def _convert_spires_date_after_to_invenio_span_query(self, query):
        """Converts date after SPIRES search term into invenio span query"""

        def create_replacement_pattern(match):
            """method used for replacement with regular expression"""
            return match.group('searchop') + ' ' + match.group('search_content') + '->9999'

        query = self._re_date_after_match.sub(create_replacement_pattern, query)

        return query

    def _convert_spires_date_before_to_invenio_span_query(self, query):
        """Converts date before SPIRES search term into invenio span query"""

        # method used for replacement with regular expression
        def create_replacement_pattern(match):
            return match.group('searchop') + ' ' + '0->' + match.group('search_content')

        query = self._re_date_before_match.sub(create_replacement_pattern, query)

        return query

    def _expand_search_patterns(self, query):
        """Expands search queries.

        If a search term is followed by several words e.g.
        author: ellis or title:THESE THREE WORDS it is expanded to
        author:ellis or (title: THESE and title:THREE...)

        Not all the search terms are expanded this way, but only title:,
        keyword:, and fulltext: - see _re_search_term_pattern_match for details.
        """

        def create_replacements(term, content):
            result = ''
            content = content.strip()

            # replace spaces within quotes by __SPACE__ temporarily:
            content = self._re_pattern_single_quotes.sub(lambda x: "'"+string.replace(x.group(1), ' ', '__SPACE__')+"'", content)
            content = self._re_pattern_double_quotes.sub(lambda x: "\""+string.replace(x.group(1), ' ', '__SPACE__')+"\"", content)
            content = self._re_pattern_regexp_quotes.sub(lambda x: "/"+string.replace(x.group(1), ' ', '__SPACE__')+"/", content)

            words = content.split()
            if len(words) > 1:
                result = '(' + term + words[0]
                for word in words[1:]:
                    result += ' and ' + term + word
                result += ')'
            else:
                result = term + words[0]
            # replace back __SPACE__ by spaces:
            result = self._re_pattern_space.sub(" ", result)
            return result.strip()

        result = ''
        current_position = 0
        for match in self._re_search_term_pattern_match.finditer(query):
            result += query[current_position : match.start()]
            result += ' ' + match.group('combine_operator') + ' '
            result += create_replacements(match.group('search_term'), match.group('search_content'))
            current_position = match.end()
        result += query[current_position : len(query)]
        return result.strip()

    def _convert_spires_truncation_to_invenio_truncation(self, query):
        """Replace SPIRES truncation symbol # with invenio trancation symbol *"""
        return query.replace('#', '*')

    def _convert_spires_exact_author_search_to_invenio_author_search(self, query):
        """Converts SPIRES search patterns for exact author into search pattern
        for invenio"""

        # method used for replacement with regular expression
        def create_replacement_pattern(match):
            # the regular expression where this group name is defined is in
            # the method _compile_regular_expressions()
            return self._EA_TAG + '"' + match.group('author_name') + '"'

        query = self._re_exact_author_match.sub(create_replacement_pattern, query)

        return query

    def _convert_spires_author_search_to_invenio_author_search(self, query):
        """Converts SPIRES search patterns for authors to search patterns in invenio
        that give similar results to the spires search.
        """

        # result of the replacement
        result = ''
        current_position = 0
        for match in self._re_author_match.finditer(query):
            result += query[current_position : match.start() ]
            scanned_name = NameScanner.scan(match.group('name'))
            author_atoms = self._create_author_search_pattern_from_fuzzy_name_dict(scanned_name)
            if author_atoms.find(' ') == -1:
                result += author_atoms + ' '
            else:
                result += '(' + author_atoms + ') '
            current_position = match.end()
        result += query[current_position : len(query)]
        return result

    def _create_author_search_pattern_from_fuzzy_name_dict(self, fuzzy_name):
        """Creates an invenio search pattern for an author from a fuzzy name dict"""

        author_name = ''
        author_middle_name = ''
        author_surname = ''
        if len(fuzzy_name['nonlastnames']) > 0:
            author_name = fuzzy_name['nonlastnames'][0]
        if len(fuzzy_name['nonlastnames']) == 2:
            author_middle_name = fuzzy_name['nonlastnames'][1]
        if len(fuzzy_name['nonlastnames']) > 2:
            author_middle_name = ' '.join(fuzzy_name['nonlastnames'][1:])
        author_surname = ' '.join(fuzzy_name['lastnames'])

        NAME_IS_INITIAL = (len(author_name) == 1)
        NAME_IS_NOT_INITIAL = not NAME_IS_INITIAL

        # we expect to have at least surname
        if author_surname == '' or author_surname == None:
            return ''

        # ellis ---> "author:ellis"
        #if author_name == '' or author_name == None:
        if not author_name:
            return self._A_TAG + author_surname

        # ellis, j ---> "ellis, j*"
        if NAME_IS_INITIAL and not author_middle_name:
            return self._A_TAG + '"' + author_surname + ', ' + author_name + '*"'

        # if there is middle name we expect to have also name and surname
        # ellis, j. r. ---> ellis, j* r*
        # j r ellis ---> ellis, j* r*
        # ellis, john r. ---> ellis, j* r* or ellis, j. r. or ellis, jo. r.
        # ellis, john r. ---> author:ellis, j* r* or exactauthor:ellis, j r or exactauthor:ellis jo r
        if author_middle_name:
            search_pattern = self._A_TAG + '"' + author_surname + ', ' + author_name + '*' + ' ' + author_middle_name.replace(" ","* ") + '*"'
            if NAME_IS_NOT_INITIAL:
                for i in range(1, len(author_name)):
                    search_pattern += ' or ' + self._EA_TAG + "\"%s, %s %s\"" % (author_surname, author_name[0:i], author_middle_name)
            return search_pattern

        # ellis, jacqueline ---> "ellis, jacqueline" or "ellis, j.*" or "ellis, j" or "ellis, ja.*" or "ellis, ja" or "ellis, jacqueline *, ellis, j *"
        # in case we don't use SPIRES data, the ending dot is ommited.
        search_pattern = self._A_TAG + '"' + author_surname + ', ' + author_name + '*"'
        search_pattern += " or " + self._EA_TAG + "\"%s, %s *\"" % (author_surname, author_name[0])
        if NAME_IS_NOT_INITIAL:
            for i in range(1,len(author_name)):
                search_pattern += ' or ' + self._EA_TAG + "\"%s, %s\"" % (author_surname, author_name[0:i])

        return search_pattern


    def _replace_spires_keywords_with_invenio_keywords(self, query):
        """Replaces SPIRES keywords that have directly
        corresponding Invenio keywords

        Replacements are done only in content that is not in quotes."""

        # result of the replacement
        result = ""
        current_position = 0

        for match in self._re_quotes_match.finditer(query):
            # clean the content after the previous quotes and before current quotes
            cleanable_content = query[current_position : match.start()]
            cleanable_content = self._replace_all_spires_keywords_in_string(cleanable_content)

            # get the content in the quotes (group one matches double
            # quotes, group 2 singles)
            if match.group(1):
                quoted_content = match.group(1)
            elif match.group(2):
                quoted_content = match.group(2)

            # append the processed content to the result
            result = result + cleanable_content + quoted_content

            # move current position at the end of the processed content
            current_position = match.end()

        # clean the content from the last appearance of quotes till the end of the query
        cleanable_content = query[current_position : len(query)]
        cleanable_content = self._replace_all_spires_keywords_in_string(cleanable_content)
        result = result + cleanable_content

        return result

    def _replace_all_spires_keywords_in_string(self, query):
        """Replaces all SPIRES keywords in the string with their
        corresponding Invenio keywords"""

        for spires_keyword, invenio_keyword in self._SPIRES_TO_INVENIO_KEYWORDS_MATCHINGS.iteritems():
            query = self._replace_keyword(query, spires_keyword, \
                                          invenio_keyword)

        return query

    def _replace_keyword(self, query, old_keyword, new_keyword):
        """Replaces old keyword in the query with a new keyword"""
        # perform case insensitive replacement with regular expression

        regex_string = r'\b(?P<operator>(find|and|or|not)\b[\s\(]*)' + \
                       old_keyword + r'(?P<end>[\s\(]+|$)'
        regular_expression = re.compile(regex_string, re.IGNORECASE)
        result = regular_expression.sub(r'\g<operator>' + new_keyword + r'\g<end>', query)
        result = re.sub(':\s+', ':', result)
        return result

    def _distribute_keywords_across_combinations(self, query):
        # method used for replacement with regular expression

        def create_replacement_pattern(match):
            # the regular expression where this group name is defined is in
            # the method _compile_regular_expressions()
            return match.group('keyword') + match.group('content') + \
                   ' ' +  match.group('combination') + ' ' + match.group('keyword') + \
                   match.group('last_content')

        query = self._re_distribute_keywords.sub(create_replacement_pattern, query)
        return query
