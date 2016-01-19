#!/usr/apps/python2.7.10/bin/python2.7

""" This script gets Chombo source code if necessary, then does a build if
necessary.  It requires Python 2.7 for argparse.

Note:
Normally the first line in this file would be:
  #!/bin/env python
but that gives us python 2.6.
This line gets python 2.7.x:
  #!/bin/env python2.7
but then PyCharm can't find the 2.7 site-packages.
"""

import argparse
import os
import pprint

import support.local_logging as local_logging
import support.runner as runner

VERSION_STRING = '1.3'

# Declaring the logger at the highest level because it must be available when
# the first decorator is parsed:
_logger = local_logging.Logger(__file__)
_Auto_Log = local_logging.getLocalAutoLog(_logger)
_logger.set_debug_on()
_logger.debug('START')


class Builder(object):
    """ Builds COGENT
    """

    @_Auto_Log('Initialize Builder', debug_only=True)
    def __init__(self):
        self.logger = _logger
        self.runner = runner.Runner()
        # <project>/<repo>/scripts:
        self.scripts_dir = os.path.dirname(os.path.realpath(__file__))
        # <project>/<repo>:
        self.repo_dir = os.path.dirname(self.scripts_dir)
        # <project>/<repo>/hypre-2.9.0b:
        self.hypre_dir = os.path.join(self.repo_dir, 'hypre-2.9.0b')
        # <project>/<repo>/exec:
        self.make_dir = os.path.join(self.repo_dir, 'exec')
        # Project dir contains other repos:
        # <project>:
        self.project_dir = os.path.dirname(self.repo_dir)
        # <project>/Chombo:
        self.chombo_dir = os.path.join(self.project_dir, 'Chombo')
        self.hypre_config = 'cab-mvapich2-gnu-2.1-dbg'
        self._define_args()

    @_Auto_Log('Run Builder')
    def run(self):
        self._process_args()
        self.logger.debug('Self is:\n"%s"' % pprint.pformat(self.__dict__))
        self._get_chombo_source()
        self._adjust_source()
        self._do_build()

    def _define_args(self):
        parser = argparse.ArgumentParser(
            description='This script checks out Chombo if necessary, and then '
                        'does a build if necessary.   ===> ',
            epilog='You only need to provide enough of each argument name to '
                   'make it unique - e.g. --i and --input_file are equivalent.',
        )
        parser.add_argument(
            '-v', '--version', action='version',
            version='%(prog)s ' + VERSION_STRING)
        # Add the rest of the arguments in alphabetical order, so they come out
        # that way in help:
        parser.add_argument(
            '--debug', action='store_true',
            dest='debug',
            help='Turn on debugging output [default: %(default)s]')
        parser.add_argument(
            '--effort_only', action='store_true',
            dest='effort_only',
            help="Log what steps would be taken, but don't actually do them. [default: %(default)s]")
        parser.add_argument(
            '--skip_build',
            '--no_build',
            action='store_true',
            dest='skip_build',
            help='Skip the build [default: %(default)s]')
        parser.add_argument(
            '--skip_get_chombo_source',
            '--skip_checkouts', '--skip_clones', '--skip_fetches', '--skip_gets',
            '--skip_pulls', '--skip_updates',
            '--no_get_chombo_source',
            '--no_checkouts', '--no_clones', '--no_fetches', '--no_gets',
            '--no_pulls', '--no_updates',
            action='store_true',
            dest='skip_checkouts',
            help='Skip getting the latest version of Chombo source [default: %(default)s]')
        parser.set_defaults(
            debug=False,
            effort_only=False,
            skip_build=False,
            skip_checkouts=False
        )
        self.parser = parser

    def _process_args(self):
        self.args = self.parser.parse_args()
        self.logger.set_debug(self.args.debug)
        self.runner.setEffortOnly(self.args.effort_only)

    @_Auto_Log('Getting Chombo source')
    def _get_chombo_source(self):
        if self.args.skip_checkouts:
            self._log_skip('Getting Chombo source')
        else:
            # Creates a dir called ./Chombo:
            self.runner.callOrLog(
                callArgs=(
                    'svn',
                    # '--non-interactive',
                    # '--trust-server-cert',
                    #'--username', 'jcompton',
                    'checkout',
                    'https://anag-repo.lbl.gov/svn/Chombo/trunk', 'Chombo'),
                dir=self.project_dir)

    @_Auto_Log('Putting our Make.defs.local in Chombo')
    def _adjust_source(self):
        """ Adds Make.defs.local to Chombo source.
        """
        self.runner.execOrLog(
            command='import shutil; shutil.copy(src="%s", dst="%s")' %
            (os.path.join(self.repo_dir, 'exec', 'Make.defs.local'),
             os.path.join(self.chombo_dir, 'lib', 'mk', 'Make.defs.local')),
            locals=locals(),
            globals=globals(),
            doReraise=True)

    @_Auto_Log('Making COGENT')
    def _do_build(self):
        if self.args.skip_build:
            self._log_skip('Making COGENT')
        else:
            self._do_hypre()
            self.runner.callOrLog(
                callArgs=('do_make',),
                dir=self.make_dir)

    @_Auto_Log('Configuring and installing hypre')
    def _do_hypre(self):
        self.runner.callOrLog(
            callArgs=('./doconfig-%s' % self.hypre_config,),
            dir=self.hypre_dir)
        self.runner.callOrLog(
            callArgs=('./doinstall',),
            dir=self.hypre_dir)
        # Below is done by doinstall:
        hypre_loc_path = os.path.join(self.hypre_dir, 'hypre_loc')
        if os.path.exists(hypre_loc_path):
            self.runner.execOrLog(
                command='os.remove("%s")' % hypre_loc_path,
                locals=locals(),
                globals=globals(),
                doReraise=True)
        self.runner.execOrLog(
            command='os.symlink("%s", "%s")' %
                    (os.path.join('lib', self.hypre_config),
                     hypre_loc_path),
            locals=locals(),
            globals=globals(),
            doReraise=True)

    # Logging utility methods:

    def _log_skip(self, section):
        self.logger.info('SKIPPING %s' % section)

if __name__ == '__main__':
    Builder().run()

_logger.debug('END')
