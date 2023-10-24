#!/usr/bin/env python3
'''
Logging module thanks to Brian Le from Pythium
'''
import logging
import re
import sys
sys.tracebacklimit = None
BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)

COLORS = {
    'WARNING': YELLOW,
    'INFO': GREEN,
    'DEBUG': BLUE,
    'CRITICAL': MAGENTA,
    'ERROR': RED
}

#These are the sequences need to get colored ouput
RESET_SEQ = "\033[0m"
COLOR_SEQ = "\033[1;%dm"
BOLD_SEQ = "\033[1m"
BLINK_SEQ = "\033[4m"

def formatter_message(message, use_color = True):
    if use_color:
        message = message.replace("$RESET", RESET_SEQ).replace("$BOLD", BOLD_SEQ)
    else:
        message = message.replace("$RESET", "").replace("$BOLD", "")
    return message

class ColoredFormatter(logging.Formatter):
    def __init__(self, msg, use_color = True):
        logging.Formatter.__init__(self, msg)
        self.use_color = use_color

    def format(self, record):
        levelname = record.levelname
        if self.use_color and levelname in COLORS:
            levelname_color = COLOR_SEQ % (30 + COLORS[levelname]) + levelname + RESET_SEQ
            record.levelname = levelname_color
            message = record.msg
            record.msg = ''
            if record.levelno == 10:
                record.msg = COLOR_SEQ % (30 + COLORS[levelname]) + record.msg
            if record.levelno > 20:
                if record.levelno > 30:
                    record.msg += COLOR_SEQ % (30 + COLORS[levelname])
                record.msg += COLOR_SEQ % (30 + COLORS[levelname])
            record.msg += COLOR_SEQ % (30 + COLORS[levelname]) + message
            if not record.levelno == 20:
                filemsg = '$RESET\t$RESET($BOLD{fl}:{fn}$RESET:{ln})'
                filemsg = filemsg.format(fl=record.filename, fn=record.funcName, ln=record.lineno)
                filemsg = filemsg.replace('$RESET', RESET_SEQ).replace('$BOLD', BOLD_SEQ)
                record.msg += filemsg
            record.msg += RESET_SEQ
        return logging.Formatter.format(self, record)

class ShutdownHandler(logging.StreamHandler):
    def emit(self, record):
        if record.levelno >= logging.ERROR:
            sys.exit(1)

class ColoredLogger(logging.Logger):
    def __init__(self, name=False):
        if name:
            self.FORMAT = "[$BOLD%(name)-s:%(levelname)-s$RESET]  %(message)s"
        else:
            self.FORMAT = "[$BOLD%(module)-s:%(levelname)-s$RESET]  %(message)s"
        self.COLOR_FORMAT = formatter_message(self.FORMAT, True)
        logging.Logger.__init__(self, name, logging.INFO)                

        color_formatter = ColoredFormatter(self.COLOR_FORMAT)

        console = logging.StreamHandler()
        console.setFormatter(color_formatter)

        self.addHandler(console)
        self.addHandler(ShutdownHandler())
        return

logging.setLoggerClass(ColoredLogger)

if __name__ == '__main__':
    logger = ColoredLogger()
    logger.debug('More verbose')
    logger.info('Info here')
    logger.warning('DANGER DANGER!!')
    logger.error('SOMEONE FUCKED UP')
    logger.critical('WHY DO YOU DO THIS????')
