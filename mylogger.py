import inspect
import sys
from datetime import datetime

INFO = 'INFO '
WARNING = 'WARN '
ERROR = 'ERROR'
FATAL = 'FATAL'
DEBUG = 'DEBUG'

class Logger:
    def __init__(self, logfile):
        self.logfile = logfile
        self.highlight_keyword = False

    def set_highlight_keyword(self, hl: bool):
        self.highlight_keyword = hl

    def cropstr(self, s, n):
        n = max(n, 4)
        return s.ljust(n) if len(s) <= n else s[:n-2] + '..'

    def log(self, severity, *message):
        filename = str(inspect.stack()[1][1]).split('/')[-1]
        filename = self.cropstr(filename, len(sys.argv[0]) - 2)
        funcname = str(inspect.stack()[2][3])
        funcname = self.cropstr(funcname, 10)

        tt = datetime.now()
        tstamp = "%s-%s-%s-%s:%s:%s" % (str(tt.year)[-2:], str(tt.month).zfill(2), str(tt.day).zfill(2), str(tt.hour).zfill(2), str(tt.minute).zfill(2), str(tt.second).zfill(2))
        coldict = {INFO: '92', WARNING: '38;5;214', ERROR: '91', FATAL: '97;41', DEBUG: '0'}

        formatted_message_parts = []
        for part in message:
            words = str(part).split()
            formatted_words = []
            for word in words:
                if self.highlight_keyword and word in [INFO,WARNING,ERROR,FATAL,DEBUG]:
                    formatted_words.append("\033[%sm%s\033[0m" % (coldict[word], word))
                else:
                    formatted_words.append(word)

            formatted_message_parts.append(' '.join(formatted_words))
        formatted_message = ' '.join(formatted_message_parts)
                
        
        print("\033[%sm[%s][%s][%s]\033[0m %s" % (coldict[severity], tstamp, severity, filename + ':' + funcname, formatted_message))

        writestring = ' '.join(map(str, message))
        with open(self.logfile, 'a') as f:
            f.write("[%s][%s][%s] %s\n" % (tstamp, severity, filename + ':' + funcname, writestring))

