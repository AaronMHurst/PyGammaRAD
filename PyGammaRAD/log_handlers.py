import logging

class FieldReporter(object):
    __doc__="""Class of methods to filter log records and handle specific log 
    levels differently."""
    log_levels = {0:"NOTSET", 10:"DEBUG", 20:"INFO",
                  30:"WARNING", 40:"ERROR", 50:"CRITICAL"}
    def show_only_critical(record):
        return record.levelname == FieldReporter.log_levels[50]
    def show_only_errors(record):
        return record.levelname == FieldReporter.log_levels[40]
    def show_only_warnings(record):
        return record.levelname == FieldReporter.log_levels[30]
    def show_only_info(record):
        return record.levelname == FieldReporter.log_levels[20]
    def show_only_debug(record):
        return record.levelname == FieldReporter.log_levels[10]
    def show_only_notset(record):
        return record.levelname == FieldReporter.log_levels[0]

class LogLevelContext:
    __doc__="""Class to handle temporary changes to the logging level for a 
    specific action or running a block of code.  The level is reverted to its 
    original state after completion of the action."""

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    
    #def __init__(self, logger, level):
    def __init__(self, level):
        self.logger = logger
        self.new_level = level
        self.original_level = None

    def __enter__(self):
        self.original_level = self.logger.level
        self.logger.setLevel(self.new_level)
        return self.logger

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.logger.setLevel(self.original_level)
        
#logging.basicConfig(level=logging.DEBUG)
#logging.basicConfig(format="%(asctime)s:%(levelname)s:%(name)s:%(message)s")
#logging.basicConfig(style="{")
#logging.basicConfig(datefmt="%Y-%m-%d %H:%M")

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter(
    "{asctime} - {name} - {levelname}: \n File: {filename} - line no.: {lineno} - function: {funcName} -\n Message:'{message}'",
    style="{",
    datefmt="%Y-%m-%d %H:%M",
)

def write_logs(level=logging.DEBUG):
    file_handler = logging.FileHandler("app_PyGammaRAD.log", mode="a", encoding="utf-8")
    file_handler.setFormatter(formatter)
    file_handler.setLevel(level)
    logger.addHandler(file_handler)
    print(f"Logging will be written to file: {file_handler}")

console_handler = logging.StreamHandler()
console_handler.setFormatter(formatter)
#console_handler.setLevel(logging.INFO)
console_handler.setLevel(logging.DEBUG)
#console_handler.addFilter(FieldReporter.show_only_errors)
#console_handler.addFilter(FieldReporter.show_only_warnings)
logger.addHandler(console_handler)


