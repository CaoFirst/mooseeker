'''
Autor: caoyh
Date: 2022-11-20 11:24:41
LastEditTime: 2022-11-20 11:31:21
'''
import logging
from logging import handlers

class Logger(object):
    level_relations = {
        'debug': logging.DEBUG,
        'info' : logging.INFO,
        'warning': logging.WARNING,
        'error': logging.ERROR,
        'critical': logging.CRITICAL
    }

    def __init__(self,filename, level='info',
                 fmt='%(asctime)s-%(pathname)s[line:%(lineno)d]-%(levelname)s : %(message)s'
                ):
        #create a logger
        self.logger = logging.getLogger()
        self.logger.setLevel(self.level_relations.get(level)) #Logger_level 总开关？？？
        format_str = logging.Formatter(fmt) #set log format

        # create a handler to input
        ch = logging.StreamHandler()
        ch.setLevel(self.level_relations.get(level))
        ch.setFormatter(format_str)

        #create a handler to filer
        fh = logging.FileHandler(filename=filename, mode='a')
        fh.setLevel(self.level_relations.get(level))
        fh.setFormatter(format_str)

        self.logger.addHandler(ch)
        self.logger.addHandler(fh)

def main():
    logger_file = 'cache/log/example.log'
    log = Logger(logger_file, level='debug')
    log.logger.debug('debug')
    log.logger.info('info')
    log.logger.warning('警告')
    log.logger.error('报错')
    log.logger.critical('严重')
    log.logger.critical('************************************************')

if __name__=='__main__':
    main()




















