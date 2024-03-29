## This file is part of Invenio.
## Copyright (C) 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010 CERN.
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

import unittest

from invenio import bibrank_dictionary
from invenio.testutils import make_test_suite, run_test_suite

class TestRankingDictionaryOperations(unittest.TestCase):
    """Test ranking dictionary operations."""

    def test_strip(self):
        """Test strip function"""
        self.assertEqual({}, bibrank_dictionary.strip({"0":0}))

TEST_SUITE = make_test_suite(TestRankingDictionaryOperations)

if __name__ == "__main__":
    run_test_suite(TEST_SUITE)
