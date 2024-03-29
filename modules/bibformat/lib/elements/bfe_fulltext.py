# -*- coding: utf-8 -*-
##
## This file is part of Invenio.
## Copyright (C) 2006, 2007, 2008, 2009, 2010, 2011 CERN.
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
"""BibFormat element - Prints a links to fulltext
"""
__revision__ = "$Id$"

import re
from invenio.bibdocfile import BibRecDocs, file_strip_ext
from invenio.messages import gettext_set_language
from invenio.config import CFG_SITE_URL, CFG_CERN_SITE
from cgi import escape
from urlparse import urlparse
from os.path import basename
import urllib

cern_arxiv_categories = ["astro-ph", "chao-dyn", "cond-mat", "gr-qc",
                         "hep-ex", "hep-lat", "hep-ph", "hep-th", "math-ph",
                         "math", "nucl-ex", "nucl-th", "out", "physics",
                         "quant-ph", "q-alg", "cs", "adap-org", "comp-gas",
                         "chem-ph", "cs", "math", "neuro-sys", "patt-sol",
                         "solv-int", "acc-phys", "alg-geom", "ao-sci",
                         "atom-ph", "cmp-lg", "dg-ga", "funct-an", "mtrl-th",
                         "plasm-ph", "q-alg", "supr-con"]

def format_element(bfo, style, separator='; ', show_icons='no', focus_on_main_file='no'):
    """
    This is the default format for formatting fulltext links.

    When possible, it returns only the main file(s) (+ link to
    additional files if needed). If no distinction is made at
    submission time between main and additional files, returns
    all the files

    @param separator: the separator between urls.
    @param style: CSS class of the link
    @param show_icons: if 'yes', print icons for fulltexts
    @param focus_on_main_file: if 'yes' and a doctype 'Main' is found,
    prominently display this doctype. In that case other doctypes are
    summarized with a link to the Files tab, named "Additional files"
    """
    _ = gettext_set_language(bfo.lang)

    out = ''

    # Retrieve files
    (parsed_urls, old_versions, additionals) = get_files(bfo, \
                                                         distinguish_main_and_additional_files=focus_on_main_file.lower() == 'yes')

    main_urls = parsed_urls['main_urls']
    others_urls = parsed_urls['others_urls']
    if parsed_urls.has_key('cern_urls'):
        cern_urls = parsed_urls['cern_urls']

    # Prepare style and icon
    if style != "":
        style = 'class="'+style+'"'

    if show_icons.lower() == 'yes':
        file_icon = '<img style="border:none" src="%s/img/file-icon-text-12x16.gif" alt="%s"/>' % (CFG_SITE_URL, _("Download fulltext"))
    else:
        file_icon = ''

    # Build urls list.
    # Escape special chars for <a> tag value.

    additional_str = ''
    if additionals:
        additional_str = ' <small>(<a '+style+' href="'+CFG_SITE_URL+'/record/'+str(bfo.recID)+'/files/">%s</a>)</small>' % _("additional files")

    versions_str = ''
    #if old_versions:
        #versions_str = ' <small>(<a '+style+' href="'+CFG_SITE_URL+'/record/'+str(bfo.recID)+'/files/">%s</a>)</small>' % _("older versions")

    if main_urls:
        last_name = ""
        main_urls_keys = sort_alphanumerically(main_urls.keys())
        for descr in main_urls_keys:
            urls = main_urls[descr]
            if re.match(r'^\d+\s', descr) and urls[0][2] == 'png':
                # FIXME: we have probably hit a Plot (as link
                # description looks like '0001 This is Caption'), so
                # do not take it.  This test is not ideal, we should
                # rather study doc type, and base ourselves on
                # Main/Additional/Plot etc.
                continue
            out += "<strong>%s:</strong> " % descr
            url_list = []
            ## FIXME: This is so ugly!
            urls_dict = {}
            for url, name, format in urls:
                urls_dict[url] = (name, format)
            urls_dict_keys = sort_alphanumerically(urls_dict.keys())
            for url in urls_dict_keys:
                name, format = urls_dict[url]
                if not name == last_name and len(main_urls) > 1:
                    print_name = "<em>%s</em> - " % name
                else:
                    print_name = ""
                last_name = name
                url_list.append(print_name + '<a '+style+' href="'+escape(url)+'">'+ \
                                file_icon + format.upper()+'</a>')
            out += separator.join(url_list) + additional_str + versions_str + '<br />'

    if CFG_CERN_SITE and cern_urls:
        link_word = len(cern_urls) == 1 and _('%(x_sitename)s link') or _('%(x_sitename)s links')
        out += '<strong>%s</strong>: ' % (link_word % {'x_sitename': 'CERN'})
        url_list = []
        for url, descr in cern_urls:
            url_list.append('<a '+style+' href="'+escape(url)+ '''" onclick="try{var r;if(typeof XMLHttpRequest=='undefined'){r=new ActiveXObject('Msxml2.XMLHTTP.3.0');}else{r=new XMLHttpRequest();}r.open('GET', \''''+CFG_SITE_URL+"/record/"+str(bfo.recID)+'''/extlink/?url=' + escape(this.href), false);r.send(null);}catch(e){}return true;">'''+ \
                            file_icon + escape(str(descr))+'</a>')
        out += separator.join(url_list)

    if others_urls:
        external_link = len(others_urls) == 1 and _('external link') or _('external links')
        out += '<strong>%s</strong>: ' % external_link.capitalize()
        url_list = []
        for url, descr in others_urls:
            url_list.append('<a '+style+' href="'+escape(url)+ '''" onclick="try{var r;if(typeof XMLHttpRequest=='undefined'){r=new ActiveXObject('Msxml2.XMLHTTP.3.0');}else{r=new XMLHttpRequest();}r.open('GET', \''''+CFG_SITE_URL+"/record/"+str(bfo.recID)+'''/extlink/?url=' + escape(this.href), false);r.send(null);}catch(e){}return true;">'''+ \
                            file_icon + escape(str(descr))+'</a>')
        out += separator.join(url_list) + '<br />'

    if out.endswith('<br />'):
        out = out[:-len('<br />')]

    # When exported to text (eg. in WebAlert emails) we do not want to
    # display the link to the fulltext:
    if out:
        out = '<!--START_NOT_FOR_TEXT-->' + out + '<!--END_NOT_FOR_TEXT-->'

    return out

def escape_values(bfo):
    """
    Called by BibFormat in order to check if output of this element
    should be escaped.
    """
    return 0

def get_files(bfo, distinguish_main_and_additional_files=True):
    """
    Returns the files available for the given record.
    Returned structure is a tuple (parsed_urls, old_versions, additionals):
     - parsed_urls: contains categorized URLS (see details below)
     - old_versions: set to True if we can have access to old versions
     - additionals: set to True if we have other documents than the 'main' document

    'parsed_urls' is a dictionary in the form:
    {'main_urls' : {'Main'      : [('http://CFG_SITE_URL/record/1/files/aFile.pdf', 'aFile', 'PDF'),
                                   ('http://CFG_SITE_URL/record/1/files/aFile.gif', 'aFile', 'GIF')],
                    'Additional': [('http://CFG_SITE_URL/record/1/files/bFile.pdf', 'bFile', 'PDF')]},

     'other_urls': [('http://externalurl.com/aFile.pdf', 'Fulltext'),      # url(8564_u), description(8564_z/y)
                    ('http://externalurl.com/bFile.pdf', 'Fulltext')],

     'cern_urls' : [('http://cern.ch/aFile.pdf', 'Fulltext'),              # url(8564_u), description(8564_z/y)
                    ('http://cern.ch/bFile.pdf', 'Fulltext')],
    }

    Some notes about returned structure:
    - key 'cern_urls' is only available on CERN site
    - keys in main_url dictionaries are defined by the BibDoc.
    - older versions are not part of the parsed urls
    - returns only main files when possible, that is when doctypes
      make a distinction between 'Main' files and other
      files. Otherwise returns all the files as main. This is only
      enabled if distinguish_main_and_additional_files is set to True
    """
    _ = gettext_set_language(bfo.lang)

    urls = bfo.fields("8564_")
    bibarchive = BibRecDocs(bfo.recID)

    old_versions = False # We can provide link to older files. Will be
                         # set to True if older files are found.
    additionals = False  # We have additional files. Will be set to
                         # True if additional files are found.

    # Prepare object to return
    parsed_urls = {'main_urls':{},    # Urls hosted by Invenio (bibdocs)
                  'others_urls':[]    # External urls
                  }
    if CFG_CERN_SITE:
        parsed_urls['cern_urls'] = [] # cern.ch urls

    # Doctypes can of any type, but when there is one file marked as
    # 'Main', we consider that there is a distinction between "main"
    # and "additional" files. Otherwise they will all be considered
    # equally as main files
    distinct_main_and_additional_files = False
    if len(bibarchive.list_bibdocs(doctype='Main')) > 0 and \
           distinguish_main_and_additional_files:
        distinct_main_and_additional_files = True

    # Parse URLs
    for complete_url in urls:
        if complete_url.has_key('u'):
            url = complete_url['u']
            (dummy, host, path, dummy, params, dummy) = urlparse(url)
            filename = urllib.unquote(basename(path))
            name = file_strip_ext(filename)
            format = filename[len(name):]
            if format.startswith('.'):
                format = format[1:]

            descr = ''
            if complete_url.has_key('y'):
                descr = complete_url['y']
            if not url.startswith(CFG_SITE_URL): # Not a bibdoc?
                if not descr: # For not bibdoc let's have a description
                    # Display the URL in full:
                    descr = url
                if CFG_CERN_SITE and 'cern.ch' in host and \
                       ('/setlink?' in url or \
                        'cms' in host or \
                        'documents.cern.ch' in url or \
                        'doc.cern.ch' in url or \
                        'preprints.cern.ch' in url):
                    url_params_dict = dict([part.split('=') for part in params.split('&') \
                                            if len(part) == 2])
                    if url_params_dict.has_key('categ') and \
                           (url_params_dict['categ'].split('.', 1)[0] in cern_arxiv_categories) and \
                           url_params_dict.has_key('id'):
                        # Old arXiv links, used to be handled by
                        # setlink. Provide direct links to arXiv
                        for file_format, label in [('pdf', "PDF")]:#,
                            #('ps', "PS"),
                            #('e-print', "Source (generally TeX or LaTeX)"),
                            #('abs', "Abstract")]:
                            url = "http://arxiv.org/%(format)s/%(category)s/%(id)s" % \
                                  {'format': file_format,
                                   'category': url_params_dict['categ'],
                                   'id': url_params_dict['id']}
                            parsed_urls['others_urls'].append((url, "%s/%s %s" % \
                                                               (url_params_dict['categ'],
                                                                url_params_dict['id'],
                                                                label)))
                else:
                    parsed_urls['others_urls'].append((url, descr)) # external url
            else: # It's a bibdoc!
                assigned = False
                for doc in bibarchive.list_bibdocs():
                    if int(doc.get_latest_version()) > 1:
                        old_versions = True
                    if True in [f.fullname.startswith(filename) \
                                for f in doc.list_all_files()]:
                        assigned = True
                        #doc.getIcon()
                        if not doc.doctype == 'Main' and \
                               distinct_main_and_additional_files == True:
                            # In that case we record that there are
                            # additional files, but don't add them to
                            # returned structure.
                            additionals = True
                        else:
                            if not descr:
                                descr = _('Fulltext')
                            if not parsed_urls['main_urls'].has_key(descr):
                                parsed_urls['main_urls'][descr] = []
                            parsed_urls['main_urls'][descr].append((url, name, format))
                if not assigned: # Url is not a bibdoc :-S
                    if not descr:
                        descr = filename
                    parsed_urls['others_urls'].append((url, descr)) # Let's put it in a general other url

    return (parsed_urls, old_versions, additionals)

_RE_SPLIT = re.compile(r"\d+|\D+")
def sort_alphanumerically(elements):
    elements = [([not token.isdigit() and token or int(token) for token in _RE_SPLIT.findall(element)], element) for element in elements]
    elements.sort()
    return [element[1] for element in elements]
