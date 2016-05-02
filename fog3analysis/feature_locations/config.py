import ConfigParser
import re
import os

def config(
    filepath="/groups/Kimble/Aman Prasad/redo_fbf/analysis/src/config.ini"):
    Config = ConfigParser.ConfigParser()
    Config.read(filepath)
    lib = ConfigSectionMap('library', Config)
    return lib


def check_exists(lib):
    for k in lib:
        if re.search('\.\w+', lib[k]) and re.search('/', lib[k]):
            li = k + ":"
            if os.path.exists(lib[k]):
                li += " ...exists (%s)\n" % lib[k]
                print li,
            else:
                li += " ...does not exist (%s)\n" % lib[k]
                print li,
                raw_input('Press key to continue anyway...')
                
 
def ConfigSectionMap(section, Config):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1
