## -*- mode: python; coding: utf-8; -*-
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
## 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.
""" Calculation of rnkUSAGEDATA and
    rnkDOWNLOADS/rnkEXTLINKS/rnkPAGEVIEWS ranks for rnkPAGEVIEWS and rnkDOWNLOADS """
import cPickle
import zlib
import datetime
import os, sys, cgi

from invenio.dbquery import run_sql
from invenio.bibrank_record_sorter import rank_records
from invenio.intbitset import intbitset as HitSet
from invenio.shellutils import run_shell_command, escape_shell_arg
from invenio.config import CFG_LOGDIR, CFG_PATH_MYSQL, CFG_DATABASE_HOST, \
     CFG_DATABASE_USER, CFG_DATABASE_PASS, CFG_DATABASE_NAME
from invenio.bibtask import task_init, write_message, task_set_option, \
        task_get_option, task_update_progress, task_update_status, \
        task_get_task_param, task_sleep_now_if_required
from invenio.config import CFG_BIBRANK_DRANK_POSTPROCESSING_TIME_WINDOW
from invenio.config import CFG_SITE_URL, CFG_SITE_SECURE_URL, CFG_CACHEDIR

# Update of the seens will happen from the "from_date" time up to the
# current time
current_time = datetime.datetime.now()
last_run_time_file = os.path.join(CFG_CACHEDIR, 'queryanalyzer_runtime')

def get_last_runtime():
    """Return time when script was last time launched (as datetime)."""
    last_run_file = open(last_run_time_file)
    last_runtime = cPickle.load(last_run_file)
    last_run_file.close()
    return last_runtime

if os.path.isfile(last_run_time_file):
    last_runtime = get_last_runtime()
else:
    # Never run this tool? Start analyzing from last year...
    last_runtime = current_time - datetime.timedelta(days = 365)
    last_run_file = open(last_run_time_file, 'w')
    cPickle.dump(last_runtime, last_run_file)
    last_run_file.close()

def _dbdump_elaborate_submit_param(key, value, dummyopts, dummyargs):
    """
    Elaborate task submission parameter.  See bibtask's
    task_submit_elaborate_specific_parameter_fnc for help.
    """
    if key in ('-r', '--rnkUD-dump'):
        try:
            task_set_option('rnkUD-dump', value)
        except ValueError:
            raise StandardError, "ERROR: Input '%s' is not valid." % \
                  value
    elif key in ('-o', '--output'):
        if os.path.isdir(value):
            task_set_option('output', value)
        else:
            raise StandardError, "ERROR: Output '%s' is not a directory." % \
                  value
    elif key in ('-f', '--type'):
        try:
            task_set_option('type', value)
        except ValueError:
            raise StandardError, "ERROR: Input '%s' is not valid." % \
                  value
    elif key in ('-b', '--bibrec'):
        try:
            task_set_option('bibrec', bool(value))
        except ValueError:
            raise StandardError, "ERROR: Input '%s' is not valid." % \
                  value
    elif key in ('-q', '--user-query'):
        try:
            task_set_option('user-query', bool(value))
        except ValueError:
            raise StandardError, "ERROR: Input '%s' is not valid." % \
                  value
    elif key in ('-w', '--wrd-sim'):
        try:
            task_set_option('wrd-sim', bool(value))
        except ValueError:
            raise StandardError, "ERROR: Input '%s' is not valid." % \
                  value
    elif key == '--site-url':
        try:
            task_set_option('site-url', value)
        except ValueError:
            raise StandardError, "ERROR: Input '%s' is not valid." % \
                  value
    elif key == '--until':
        try:
            until_date = datetime.datetime.strptime(value, '%Y-%m-%d')
            task_set_option('until', until_date)
        except ValueError:
            raise StandardError, "ERROR: Input '%s' is not valid." % \
                  value
    else:
        return False
    return True

def dump_bibrec(dump_folder, bibrec_dump, dump_type="sql"):
    """Dumping bibrec"""

    dump_bibrec_folder = dump_folder + 'bibrec'
    if bibrec_dump and dump_type == 'sql':
        if not os.path.isdir(dump_bibrec_folder):
            os.makedirs(dump_bibrec_folder)
        cmd = CFG_PATH_MYSQL + 'dump'
        if not os.path.exists(cmd):
            write_message("ERROR: cannot find %s." % cmd, stream=sys.stderr)
            task_update_status("ERROR")
            sys.exit(1)

        dump_file = '/drank_dump_bibrec_' + current_time.strftime("%Y%m%d%H%M%S") + '.sql'

        write_message("... writing %s" % dump_bibrec_folder + dump_file)

        cmd += " --skip-opt --add-drop-table --add-locks --create-options " \
               " --quick --extended-insert --set-charset --disable-keys " \
               " --host=%s --user=%s --password=%s %s bibrec" % \
               (escape_shell_arg(CFG_DATABASE_HOST),
                escape_shell_arg(CFG_DATABASE_USER),
                escape_shell_arg(CFG_DATABASE_PASS),
                escape_shell_arg(CFG_DATABASE_NAME))

        dummy1, dummy2, dummy3 = run_shell_command(cmd, None, dump_bibrec_folder + dump_file)
        if dummy1:
            write_message("ERROR: mysqldump exit code is %s." % repr(dummy1),
                          stream=sys.stderr)
            task_update_status("ERROR")
            sys.exit(1)
        if dummy2:
            write_message("ERROR: mysqldump stdout is %s." % repr(dummy1),
                          stream=sys.stderr)
            task_update_status("ERROR")
            sys.exit(1)
        if dummy3:
            write_message("ERROR: mysqldump stderr is %s." % repr(dummy1),
                          stream=sys.stderr)
            task_update_status("ERROR")
            sys.exit(1)
    elif bibrec_dump and dump_type == 'csv':
        if not os.path.isdir(dump_bibrec_folder):
            os.makedirs(dump_bibrec_folder)

        dump_file = '/drank_dump_bibrec_' + current_time.strftime("%Y%m%d%H%M%S") + '.csv'

        write_message("... writing %s" % dump_bibrec_folder + dump_file)

        run_sql("""SELECT * INTO OUTFILE %s"""\
                """ FIELDS TERMINATED BY ',' OPTIONALLY ENCLOSED BY '"'"""\
                """ LINES TERMINATED BY '\n'"""\
                """ FROM bibrec""",\
                 ('/tmp'+dump_file,)) #'

        dummy1, dummy2, dummy3 = run_shell_command("cp /tmp/%s %s", (dump_file, dump_bibrec_folder))
        if dummy1:
            write_message("ERROR: mysqldump exit code is %s." % repr(dummy1),
                          stream=sys.stderr)
            task_update_status("ERROR")
            sys.exit(1)
        if dummy2:
            write_message("ERROR: mysqldump stdout is %s." % repr(dummy1),
                          stream=sys.stderr)
            task_update_status("ERROR")
            sys.exit(1)
        if dummy3:
            write_message("ERROR: mysqldump stderr is %s." % repr(dummy1),
                          stream=sys.stderr)
            task_update_status("ERROR")
            sys.exit(1)

def dump_rnkusagedata(dump_folder, full_dump, dump_type='sql'):
    """
    Dump rnkUSAGEDATA table into SQL file called FILENAME living in
    DIRNAME. If full_dump is specified we dump whole table, otherwise
    we dump only those rows which were modified from the last time
    script has been run.
    """
    dump_rnkUSAGEDATA_folder = dump_folder + 'rnkUSAGEDATA/'
    if not os.path.isdir(dump_rnkUSAGEDATA_folder):
        os.makedirs(dump_rnkUSAGEDATA_folder)
        full_dump_folder = True
    else:
        full_dump_folder = False

    if dump_type == 'sql':
        cmd = CFG_PATH_MYSQL + 'dump'
        if not os.path.exists(cmd):
            write_message("ERROR: cannot find %s." % cmd, stream=sys.stderr)
            task_update_status("ERROR")
            sys.exit(1)

        dump_file = 'drank_dump_rnkUSAGEDATA_' + current_time.strftime("%Y%m%d%H%M%S") + '.sql'

        if full_dump=='1' or full_dump_folder:
            write_message("... writing FULL DUMP %s" % dump_rnkUSAGEDATA_folder + dump_file)

            cmd += " --skip-opt --add-drop-table --add-locks --create-options " \
                   " --quick --extended-insert --set-charset --disable-keys " \
                   " --host=%s --user=%s --password=%s %s rnkUSAGEDATA" % \
                   (escape_shell_arg(CFG_DATABASE_HOST),
                    escape_shell_arg(CFG_DATABASE_USER),
                    escape_shell_arg(CFG_DATABASE_PASS),
                    escape_shell_arg(CFG_DATABASE_NAME))
        else:
            write_message("... writing INCREMENTAL DUMP %s" % dump_rnkUSAGEDATA_folder + dump_file)

            cmd += " --skip-opt --add-drop-table --add-locks --create-options " \
                   " --quick --extended-insert --set-charset --disable-keys " \
                   " --host=%s --user=%s --password=%s %s rnkUSAGEDATA --where='time_stamp>%s'" % \
                   (escape_shell_arg(CFG_DATABASE_HOST),
                    escape_shell_arg(CFG_DATABASE_USER),
                    escape_shell_arg(CFG_DATABASE_PASS),
                    escape_shell_arg(CFG_DATABASE_NAME),
                    'current_time')

        dummy1, dummy2, dummy3 = run_shell_command(cmd, None, dump_rnkUSAGEDATA_folder + dump_file)
        if dummy1:
            write_message("ERROR: mysqldump exit code is %s." % repr(dummy1),
                          stream=sys.stderr)
            task_update_status("ERROR")
            sys.exit(1)
        if dummy2:
            write_message("ERROR: mysqldump stdout is %s." % repr(dummy1),
                          stream=sys.stderr)
            task_update_status("ERROR")
            sys.exit(1)
        if dummy3:
            write_message("ERROR: mysqldump stderr is %s." % repr(dummy1),
                          stream=sys.stderr)
            task_update_status("ERROR")
            sys.exit(1)
    elif dump_type == 'csv':

        dump_file = '/drank_dump_rnkUSAGEDATA_' + current_time.strftime("%Y%m%d%H%M%S") + '.csv'

        if full_dump=='1' or full_dump_folder:
            write_message("... writing FULL DUMP %s" % '/tmp' + dump_file)

            run_sql("""SELECT * INTO OUTFILE %s"""\
                    """ FIELDS TERMINATED BY ',' OPTIONALLY ENCLOSED BY '"'"""\
                    """ LINES TERMINATED BY '\n'"""\
                    """ FROM rnkUSAGEDATA""",\
                     ('/tmp' + dump_file,))
        else:
            write_message("... writing INCREMENTAL DUMP %s" % '/tmp' + dump_file)
            run_sql("""SELECT * INTO OUTFILE %s"""\
                    """ FIELDS TERMINATED BY ',' OPTIONALLY ENCLOSED BY '"'"""\
                    """ LINES TERMINATED BY '\n'"""\
                    """ FROM rnkUSAGEDATA WHERE time_stamp>%s""",\
                     ('/tmp' + dump_file, current_time))

        dummy1, dummy2, dummy3 = run_shell_command("cp /tmp/%s %s", (dump_file, dump_rnkUSAGEDATA_folder))
        if dummy1:
            write_message("ERROR: mysqldump exit code is %s." % repr(dummy1),
                          stream=sys.stderr)
            task_update_status("ERROR")
            sys.exit(1)
        if dummy2:
            write_message("ERROR: mysqldump stdout is %s." % repr(dummy1),
                          stream=sys.stderr)
            task_update_status("ERROR")
            sys.exit(1)
        if dummy3:
            write_message("ERROR: mysqldump stderr is %s." % repr(dummy1),
                          stream=sys.stderr)
            task_update_status("ERROR")
            sys.exit(1)

def update_display_counts(reclist_in_collection):
    """number of displayes is updated here"""
    if reclist_in_collection:
        data = [(row, 1, str(current_time)) for row in reclist_in_collection]

        run_sql("INSERT INTO rnkUSAGEDATA"\
            " (id_bibrec, nbr_displays, time_stamp)"\
            " VALUES %s"\
            " ON DUPLICATE KEY UPDATE nbr_displays=nbr_displays+1, time_stamp='%s'" % (str(data)[1:-1], str(current_time)))

def update_seens_counts(action_recids, reclist_in_collection):
    '''number of seens is updated here'''
    records_position = {}
    if action_recids and reclist_in_collection:
        for recid in action_recids:
            if recid in reclist_in_collection:
                reclist_printed = reclist_in_collection[:]
                record_position = reclist_printed.index(recid)
                records_position[recid] = record_position
    else:
        reclist_printed = []
    if records_position:
        last_clicked_record_id = max(records_position, key=records_position.get)
        last_clicked_position = records_position[last_clicked_record_id] + 1
    else:
        last_clicked_record_id = 0
        last_clicked_position = 0
        reclist_printed = []

    reclist_seen = reclist_printed[:last_clicked_position]

    if len(reclist_seen) == 1:
        run_sql("UPDATE rnkUSAGEDATA"\
                " SET nbr_seens=nbr_seens+1, time_stamp=%s"\
                " WHERE id_bibrec=%s",\
                (current_time, reclist_seen[0]))
    elif len(reclist_seen)>1:
        run_sql("UPDATE rnkUSAGEDATA"\
                " SET nbr_seens=nbr_seens+1, time_stamp=%s"\
                " WHERE id_bibrec IN %s",\
                (current_time, tuple(reclist_seen)))

def get_viewed_recids(ip_address, user_id, referrer, window_time, query_time):
    "Returns record ids which were detail viewed if any, otherwise returns empty list."
    viewed_records = []
    viewed_recids = run_sql("""SELECT id_bibrec
                               FROM rnkPAGEVIEWS
                               WHERE client_host=%s
                               AND view_time BETWEEN %s AND %s
                               AND referer=%s
                               AND id_user=%s""",
                             (ip_address, query_time, window_time, referrer, user_id))
    return [row[0] for row in viewed_recids]

def get_extdownloaded_recids(viewed_records, ip_address, user_id, referrer, window_time, query_time):
    "gets record ids which were downloaded"
    extdownloaded_recids = []
    if not viewed_records:
        db_args = (referrer, query_time, window_time, ip_address, user_id)
        extdownloaded_recids = run_sql("SELECT id_bibrec"\
                                    " FROM rnkEXTLINKS"\
                                    " WHERE referer=%s"\
                                    " AND download_time BETWEEN %s and %s"\
                                    " AND client_host=%s"\
                                    " AND id_user=%s",\
                                     db_args)
    else:
        if len(viewed_records)==1:
            db_args = (viewed_records[0], referrer, query_time, window_time, ip_address, user_id)
            extdownloaded_recids = run_sql("SELECT id_bibrec"\
                                        " FROM rnkEXTLINKS"\
                                        " WHERE (id_bibrec=%s OR referer=%s)"\
                                        " AND download_time BETWEEN %s and %s"\
                                        " AND client_host=%s"\
                                        " AND id_user=%s",\
                                         db_args)
        elif len(viewed_records)>1:
            db_args = (tuple(viewed_records), referrer, query_time, window_time, ip_address, user_id)
            extdownloaded_recids = run_sql("SELECT id_bibrec"\
                                        " FROM rnkEXTLINKS"\
                                        " WHERE (id_bibrec IN %s OR referer=%s)"\
                                        " AND download_time BETWEEN %s and %s"\
                                        " AND client_host=%s"\
                                        " AND id_user=%s",\
                                         db_args)
    return [row[0] for row in extdownloaded_recids]

def get_downloaded_recids(viewed_records, ip_address, user_id, referrer, window_time, query_time):
    "gets record ids which were downloaded"
    downloaded_recids = []
    if not viewed_records:
        db_args = (referrer, query_time, window_time, ip_address, user_id)
        downloaded_recids = run_sql("SELECT id_bibrec"\
                                    " FROM rnkDOWNLOADS"\
                                    " WHERE referer=%s"\
                                    " AND download_time BETWEEN %s AND %s"\
                                    " AND client_host=%s"\
                                    " AND id_user=%s",\
                                     db_args)
    else:
        if len(viewed_records)==1:
            db_args = (viewed_records[0], referrer, query_time, window_time, ip_address, user_id)
            downloaded_recids = run_sql("SELECT id_bibrec"\
                                        " FROM rnkDOWNLOADS"\
                                        " WHERE (id_bibrec=%s OR referer=%s)"\
                                        " AND download_time BETWEEN %s AND %s"\
                                        " AND client_host=%s"\
                                        " AND id_user=%s",\
                                         db_args)
        elif len(viewed_records)>1:
            db_args = (tuple(viewed_records), referrer, user_id, ip_address, query_time, window_time)
            downloaded_recids = run_sql("SELECT id_bibrec"\
                                        " FROM rnkDOWNLOADS"\
                                        " WHERE (id_bibrec IN %s OR referer=%s)"\
                                        " AND id_user=%s"\
                                        " AND client_host=%s"\
                                        " AND download_time BETWEEN %s AND %s",\
                                         db_args)
    return [row[0] for row in downloaded_recids]

def update_rnkPAGEVIEWS(viewed_recids, user_id, ip_address, reclist, referrer, jrec, query_date, query_time_window, query_id):
    """Updates detail viewed document ranks."""
    if reclist and viewed_recids:
        for viewed_recid in viewed_recids:
            if viewed_recid in reclist:
                if jrec:
                    view_position = reclist.index(viewed_recid)+jrec
                else:
                    view_position = reclist.index(viewed_recid)+1
                run_sql("UPDATE rnkPAGEVIEWS"\
                        " SET display_position=%s, view_delay=TIMEDIFF(view_time, %s), id_query=%s"\
                        " WHERE id_bibrec=%s"\
                        " AND view_time BETWEEN %s AND %s"\
                        " AND id_user=%s"\
                        " AND referer=%s"\
                        " AND client_host=%s",\
                        (view_position, query_date, query_id, viewed_recid, query_date, query_time_window, user_id, referrer, ip_address))

def update_rnkDOWNLOADS(downloaded_recids, user_id, ip_address, reclist, referrer, jrec, query_date, query_time_window, query_id, site_url):
    """Updates downloaded document ranks."""
    if downloaded_recids and reclist:
        for downloaded_recid in downloaded_recids:
            if downloaded_recid in reclist:
                if jrec:
                    download_position = reclist.index(downloaded_recid)+jrec
                else:
                    download_position = reclist.index(downloaded_recid)+1

                run_sql("UPDATE rnkDOWNLOADS"\
                        " SET display_position=%s, download_delay=TIMEDIFF(download_time, %s), id_query=%s"\
                        " WHERE id_bibrec=%s"\
                        " AND download_time BETWEEN %s AND %s"\
                        " AND id_user=%s"\
                        " AND client_host=%s"\
                        " AND (referer=%s OR referer LIKE %s)",\
                         (download_position, query_date, query_id, downloaded_recid, query_date, query_time_window, user_id, ip_address, referrer, site_url+'/record/'+str(downloaded_recid)+'%'))

def update_rnkEXTLINKS(extdownloaded_recids, user_id, ip_address, reclist, referrer, jrec, query_date, query_time_window, query_id, site_url):
    """Updates downloaded document ranks."""
    if extdownloaded_recids and reclist:
        for extdownloaded_recid in extdownloaded_recids:
            if extdownloaded_recid in reclist:
                if jrec:
                    download_position = reclist.index(extdownloaded_recid)+jrec
                else:
                    download_position = reclist.index(extdownloaded_recid)+1

                run_sql("UPDATE rnkEXTLINKS"\
                        " SET display_position=%s, download_delay=TIMEDIFF(download_time, %s), id_query=%s"\
                        " WHERE id_bibrec=%s"\
                        " AND download_time BETWEEN %s AND %s"\
                        " AND id_user=%s"\
                        " AND client_host=%s"\
                        " AND (referer=%s OR referer LIKE %s)",\
                         (download_position, query_date, query_id, extdownloaded_recid, query_date, query_time_window, user_id, ip_address, referrer, site_url+'/record/'+str(extdownloaded_recid)+'%'))

def getWrdSimilarity(reclist, patterns, query_id, dump_folder):
    """This function inserts word similarity for each entry of reclist in external file. """
    hitset = HitSet(reclist)
    rank_values = rank_records("wrd", 0, hitset, patterns, 0)
    records, scores = rank_values[0:2]

    wrd_sim_data = repr(query_id) + '---' + repr(records) + '---' + repr(scores) + '\n'

    dump_wrd_sim_folder = dump_folder + 'wrd_sim'
    print dump_wrd_sim_folder

    dump_wrd_sim_file = 'drank_wrd_sim_dump_' + current_time.strftime("%Y%m%d%H%M%S") + '.csv'

    if not os.path.isdir(dump_wrd_sim_folder):
        os.makedirs(dump_wrd_sim_folder)

    if os.path.isfile(dump_wrd_sim_folder + dump_wrd_sim_file):
        fd = open(dump_wrd_sim_folder + dump_wrd_sim_file,'a')
        fd.write(wrd_sim_data)
        fd.close()
    else:
        fd = open(dump_wrd_sim_folder + dump_wrd_sim_file,'w')
        write_message("... writing %s" % dump_wrd_sim_folder+dump_wrd_sim_file)
        fd.write(wrd_sim_data)
        fd.close()

def dump_user_query(dump_folder, user_id, ip_address, query_id, date, pattern, reclist):
    """Dumping user_query information"""
    user_query_data = repr(user_id) + '---' + repr(ip_address) + '---' +  repr(query_id) + '---' + repr(date) + '---' + repr(pattern) + '---' + repr(reclist) + '\n'

    user_query_dump_folder = dump_folder + 'user_query'

    user_query_file = '/drank_user_query_dump_' + current_time.strftime("%Y%m%d%H%M%S") + '.csv'

    if not os.path.isdir(user_query_dump_folder):
        os.makedirs(user_query_dump_folder)

    if os.path.isfile(user_query_dump_folder + user_query_file):
        fd = open(user_query_dump_folder + user_query_file,'a')
        fd.write(user_query_data)
        fd.close()
    else:
        fd = open(user_query_dump_folder + user_query_file,'w')
        write_message("... writing %s" % user_query_dump_folder + user_query_file)
        fd.write(user_query_data)
        fd.close()

def get_jrec(args):
    "Returns jrec if any, otherwise None."
    if args.has_key('jrec'):
        try:
            jrec = int(args['jrec'][0])
        except:
            jrec = None
        return jrec
    else:
        return None

def unify(sequence):
    "Takes list as input and returns unified list"
    if sequence:
        yada = set(sequence)
        return list(yada)

def update_download_counts(from_date, until_date):
    """
    Update download counts in the rnkUSAGEDATA table
    """
    if until_date:
        data_downloads = run_sql("""SELECT id_bibrec, CONCAT(COUNT(*)), %s
                                    FROM  rnkDOWNLOADS
                                    WHERE id_query > 0
                                    AND   download_time>=%s
                                    AND   download_time<=%s
                                    GROUP BY id_bibrec""", (current_time, from_date, until_date))

        data_extdownloads = run_sql("""SELECT id_bibrec, CONCAT(COUNT(*)), %s
                                       FROM  rnkEXTLINKS
                                       WHERE id_query > 0
                                       AND   download_time>=%s
                                       AND   download_time<=%s
                                       GROUP BY id_bibrec""", (current_time, from_date, until_date))
    else:
        data_downloads = run_sql("""SELECT id_bibrec, CONCAT(COUNT(*)), %s
                                    FROM  rnkDOWNLOADS
                                    WHERE id_query > 0
                                    AND   download_time>=%s
                                    GROUP BY id_bibrec""", (current_time, from_date))

        data_extdownloads = run_sql("""SELECT id_bibrec, CONCAT(COUNT(*)), %s
                                       FROM  rnkEXTLINKS
                                       WHERE id_query > 0
                                       AND   download_time>=%s
                                       GROUP BY id_bibrec""", (current_time, from_date))

    data = data_downloads + data_extdownloads

    if data:
        for x in range(len(data)/100000+1):
            run_sql("INSERT INTO rnkUSAGEDATA"\
                    " (id_bibrec, nbr_downloads, time_stamp)"\
                    " VALUES %s"\
                    " ON DUPLICATE KEY UPDATE"\
                    " nbr_downloads=nbr_downloads+VALUES(nbr_downloads), time_stamp=VALUES(time_stamp)" % \
                    (str(data[x*100000:(x+1)*100000])[1:-1],))

def update_pageview_counts(from_date, until_date):
    """
    Update detailed view counts in the rnkUSAGEDATA table. Update is
    done in chuncks, one chunk contains 100 000 entries.
    """
    if until_date:
        data = run_sql("""SELECT id_bibrec, CONCAT(COUNT(*)), %s
                            FROM  rnkPAGEVIEWS
                            WHERE id_query > 0
                            AND   view_time>=%s
                            AND   view_time<=%s
                            GROUP BY id_bibrec""", (current_time, from_date, until_date))
    else:
        data = run_sql("""SELECT id_bibrec, CONCAT(COUNT(*)), %s
                            FROM  rnkPAGEVIEWS
                            WHERE id_query > 0
                            AND   view_time>=%s
                            GROUP BY id_bibrec""", (current_time, from_date))

    if data:
        for x in range(len(data)/100000+1):
            run_sql("INSERT INTO rnkUSAGEDATA"\
                    " (id_bibrec, nbr_pageviews, time_stamp)"\
                    " VALUES %s"\
                    " ON DUPLICATE KEY UPDATE"\
                    " nbr_pageviews=nbr_pageviews+VALUES(nbr_pageviews), time_stamp=VALUES(time_stamp)" % \
                    (str(data[x*100000:(x+1)*100000])[1:-1],))

def get_drank_counts():
    """Core function for updating rnkUSAGEDATA counts and dumping necessary information."""
#    startime = datetime.datetime.now()

    # Getting some parameters
    dump_folder = task_get_option('output', CFG_LOGDIR + '/drank_dumps/')
    rnkUD_dump_type = task_get_option('rnkUD-dump', 0)
    dump_type = task_get_option('type', 'sql')
    bibrec_dump = task_get_option('bibrec', bool(False))
    user_query_dump = task_get_option('user-query', bool(False))
    wrd_sim_dump = task_get_option('wrd-sim', bool(False))
    verbose = task_get_option('verbose', 0)

    site_url = task_get_option('site-url', CFG_SITE_URL)
    until_date = task_get_option('until', None)

    task_update_progress("Starting ...")
    write_message("Starting ...")

    # The following restrict analysis to queries run from this website,
    # and ignore all other queries (such as those coming from an external
    # website redirecting to this search, or all other searches done through
    # the search APIs) in order to reduce the amount of meaningless queries
    # to analyze (we run into memory trouble when removing this constraint:
    # further analysis as to if 'referer != null' would be sufficient should
    # be carried on)
    if until_date:
        res = run_sql("""SELECT id_user, id_query, hostname, date, reclist, referer, client_host
                         FROM user_query
                         WHERE date>=%s
                         AND date<=%s
                         AND referer LIKE %s""",
                      (last_runtime, until_date, site_url + '%'))
    else:
        res = run_sql("""SELECT id_user, id_query, hostname, date, reclist, referer, client_host
                         FROM user_query
                         WHERE date>=%s
                         AND referer LIKE %s""",
                      (last_runtime, site_url + '%'))
    write_message("Found %i user queries to process" % len(res), verbose=5)
    user_ids = [row[0] for row in res]
    query_ids = [row[1] for row in res]
    #hostnames = [row[2] for row in res]
    dates = [row[3] for row in res]
    reclists = [cPickle.loads(zlib.decompress(row[4])) for row in res]
    #referers = [row[5] for row in res]
    ip_addresses = [row[6] for row in res]

    # We are going to update tables. Store on disk the last runtime of
    # the script before.
    last_run_file = open(last_run_time_file, 'w')
    if until_date:
        # If we specify date range, we don't want to update with
        # current time, but specified date (only if specifed date is
        # newer)
        if until_date > last_runtime:
            cPickle.dump(until_date, last_run_file)
    else:
        cPickle.dump(current_time, last_run_file)
    last_run_file.close()

    write_message("Processing Task 1/4: user_query analysis")
    task_update_progress("Processing Task 1/4: user_query analysis")

    for i, query_id in enumerate(query_ids):
        query_args = run_sql("SELECT urlargs FROM query WHERE id=%s",
                              (query_id,))[0][0]
        merged_recids = []
        referrer = site_url + '/search?' + query_args
        args = cgi.parse_qs(query_args)
        jrec = get_jrec(args)
        query_time_window = dates[i] + datetime.timedelta(minutes = CFG_BIBRANK_DRANK_POSTPROCESSING_TIME_WINDOW)

        viewed_recids = get_viewed_recids(ip_addresses[i], user_ids[i], referrer, query_time_window, dates[i])
        extdownloaded_recids = get_extdownloaded_recids(viewed_recids, ip_addresses[i], user_ids[i], referrer, query_time_window, dates[i])
        downloaded_recids = get_downloaded_recids(viewed_recids, ip_addresses[i], user_ids[i], referrer, query_time_window, dates[i])

        action_recids = unify(downloaded_recids + viewed_recids + extdownloaded_recids)

        task_update_progress("Processing Task 1/4: user_query entry %d out of %d" % (i+1, len(query_ids)))
#        write_message("Processing Task 1/4: user_query entry %d out of %d." % (i+1, len(query_ids)))

        for record_in_list in reclists[i]:
            # Name of the collection
            #collection_name = record_in_list[0]
            # Records list in the collection
            reclist_in_collection = record_in_list[1]
            # Extra 20 record ids, which were not displayed
            extra_reclist_in_collection = record_in_list[2]

            # Rank and delay from query of viewed record is updated here
            update_rnkPAGEVIEWS(viewed_recids, user_ids[i], ip_addresses[i],
                                reclist_in_collection, referrer, jrec, dates[i],
                                query_time_window, query_ids[i])

            # Rank and delay from query of downloaded records is updated here
            update_rnkDOWNLOADS(downloaded_recids, user_ids[i], ip_addresses[i],
                                reclist_in_collection, referrer, jrec, dates[i],
                                query_time_window, query_ids[i], site_url)
            update_rnkEXTLINKS(extdownloaded_recids, user_ids[i], ip_addresses[i],
                               reclist_in_collection, referrer, jrec, dates[i],
                               query_time_window, query_ids[i], site_url)

            # Here we run update on displays and seens counts
            update_display_counts(reclist_in_collection)
            update_seens_counts(action_recids, reclist_in_collection)

            # Here we collect record ids and patterns for calculating wrd similarity
            pattern = []
            if args.has_key('p') and reclist_in_collection:
                pattern_string = cgi.parse_qs(query_args)['p'][0]
                pattern = pattern_string.split(' ')
                merged_recids.extend(reclist_in_collection)
                merged_recids.extend(extra_reclist_in_collection)

        # Here we call wrd similarity and dump user_query table
        if pattern:
            if wrd_sim_dump:
                getWrdSimilarity(merged_recids, pattern, query_id, dump_folder)
            if user_query_dump:
                dump_user_query(dump_folder, user_ids[i], ip_addresses[i],
                                query_id, dates[i], pattern, reclists[i])
        task_sleep_now_if_required()

    write_message("Done Task 1/4: Processed user_query entries")
    task_update_progress("Done Task 1/4: Processed user_query entries")

    task_update_progress("Starting Task 2/4: Updating download counts")
    write_message("Starting Task 2/4: Updating download counts")
    task_sleep_now_if_required()
    update_download_counts(last_runtime, until_date)
    write_message("Done Task 2/4: Updating download counts")

    task_update_progress("Starting Task 3/4: Updating page view counts")
    write_message("Starting Task 3/4: Updating page view counts")
    task_sleep_now_if_required()
    update_pageview_counts(last_runtime, until_date)
    write_message("Done Task 3/4: Updating page view counts")

    task_update_progress("Starting Task 4/4: Dumping ...")
    write_message("Starting Task 4/4: Dumping ...")
    task_sleep_now_if_required(can_stop_too=True)
    if rnkUD_dump_type:
        dump_rnkusagedata(dump_folder, rnkUD_dump_type, dump_type)
    dump_bibrec(dump_folder, bibrec_dump, dump_type)
    write_message("Done Task 4/4: Dumping.")
    write_message("Done")
    task_update_progress("Done")

#    endtime = datetime.datetime.now()
#    total_time = endtime-startime
#    print "TOTAL RUNNING TIME:", total_time
    return 1

def main():
    """main function which is constructing bibtask."""

    task_init(authorization_action='run_drank_dump',
              authorization_msg="Query Analyzer",
              help_specific_usage="""\
  -o, --output=DIR      Output directory [default=%s]
  -r, --rnkUD-dump=FD   Should we dump full rnkUSAGEDATA table(1), only changes(2) from last dump or not to dump it(0) [default=0]
  -f, --type=FT         Dump output format (csv,sql) [default=sql]
  -b, --bibrec=BR       Should we dump bibrec table (True) or not (False) [default=false]
  -q, --user-query=UQ   Should we dump user_query data (True) or not (False) [default=false]
  -w, --wrd-sim=WS      Should we dump word similarity data (True) or not (False) [default=false]
      --site-url=URL    Site URL, if analyzing imported logs [default=%s]
      --until=DATE      Analyze logs until given DATE (format YYYY-MM-DD)

Examples:
Run periodically (for eg. weekly basis) and dump last changes from rnkUSAGEDATA:
   $ queryanalyzer --rnkUD-dump=2  -u admin -s 1d -L "Sunday 01:00-05:00"
""" % (CFG_LOGDIR + '/drank_dumps', CFG_SITE_URL),
              specific_params=("r:o:f:b:q:w:",
                               ["rnkUD-dump=", "output=", "type=", "bibrec=",
                                "user-query=", "wrd-sim=", "site-url=", "until="]),
              task_submit_elaborate_specific_parameter_fnc = _dbdump_elaborate_submit_param,
              task_run_fnc = get_drank_counts)

if __name__ == '__main__':
    main()
