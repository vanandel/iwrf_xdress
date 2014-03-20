package = 'iwrf'
packagedir = 'iwrf'

from xdress.utils import DEFAULT_PLUGINS
from xdress.utils import apiname

plugins = list(DEFAULT_PLUGINS) + ['xdress.autoall',]

#xdress.autoall.clang_includes = ('../include',)
classes = [apiname('iwrf_ts_processing', 
                   ('../iwrf_lib/iwrf_functions.c'),
                   incfiles='../include/iwrf_functions.h'),
apiname('iwrf_packet_info', 
                   ('../iwrf_lib/iwrf_functions.c'),
                   incfiles='../include/iwrf_data.h')]

functions = [apiname('iwrf_ts_processing_init',('../iwrf_lib/iwrf_functions.c'),
             incfiles= '../include/iwrf_functions.h')]

