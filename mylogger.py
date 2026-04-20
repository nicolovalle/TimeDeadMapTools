import inspect
import sys
from datetime import datetime

INFO = 'INFO '
IMPORTANT = 'IMPO '
WARNING = 'WARN '
ERROR = 'ERROR'
FATAL = 'FATAL'
DEBUG = 'DEBUG'
WARN = WARNING
IMPO = IMPORTANT

SEVERITY_LEVELS = {
    DEBUG: 10,
    INFO: 20,
    WARNING: 30,
    ERROR: 40,
    FATAL: 50,
    IMPORTANT: 100
}


def is_higher_severity(s1, s2):
    return SEVERITY_LEVELS[s1] > SEVERITY_LEVELS[s2]
def is_atleast_severity(s1, s2):
    return SEVERITY_LEVELS[s1] >= SEVERITY_LEVELS[s2]

# TODO: import newly introduced severity levels
class Logger:
    def __init__(self, logfile):
        self.logfile = logfile
        self.highlight_keyword = False
        self.print_on_terminal = True # IMPORTANT are printer on terminal even if False
        self.verbosity = 1 # -1: suppress all, 0: only IMPORTANT, 1: no DEBUG, 2: also DEBUG

    def set_highlight_keyword(self, hl: bool):
        self.highlight_keyword = hl

    def set_print_on_terminal(self, tr: bool):
        self.print_on_terminal = tr

    def set_verbosity(self, vb: int):
        self.verbosity = vb

    def cropstr(self, s, n):
        n = max(n, 4)
        return s.ljust(n) if len(s) <= n else s[:n-2] + '..' 

    def log(self, severity, *message):

        if self.verbosity < 2 and severity == DEBUG:
            return
        if self.verbosity < 1 and severity != IMPORTANT:
            return
        if self.verbosity < 0:
            return

        filename = str(inspect.stack()[1][1]).split('/')[-1]
        filename = self.cropstr(filename, len(sys.argv[0]) - 2)
        funcname = str(inspect.stack()[2][3])
        funcname = self.cropstr(funcname, 10)

        tt = datetime.now()
        tstamp = "%s%s%s-%s:%s:%s" % (str(tt.year)[-2:], str(tt.month).zfill(2), str(tt.day).zfill(2), str(tt.hour).zfill(2), str(tt.minute).zfill(2), str(tt.second).zfill(2))
        coldict = {INFO: '92', WARNING: '38;5;214', ERROR: '91', FATAL: '97;41', DEBUG: '0', IMPORTANT: '30;104'}

        formatted_message_parts = []
        for part in message:
            words = str(part).split()
            formatted_words = []
            for word in words:
                if self.highlight_keyword and word in coldict.keys():
                    formatted_words.append("\033[%sm%s\033[0m" % (coldict[word], word))
                else:
                    formatted_words.append(word)

            formatted_message_parts.append(' '.join(formatted_words))
        formatted_message = ' '.join(formatted_message_parts)
                
        if self.print_on_terminal or severity == IMPORTANT:
            print("\033[%sm[%s][%s][%s]\033[0m %s" % (coldict[severity], tstamp, severity, filename + ':' + funcname, formatted_message))

        writestring = ' '.join(map(str, message))
        with open(self.logfile, 'a') as f:
            f.write("[%s][%s][%s] %s\n" % (tstamp, severity, filename + ':' + funcname, writestring))

