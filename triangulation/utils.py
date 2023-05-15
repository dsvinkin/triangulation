import os
import glob

def search_file(file_name):

    lst_files = glob.glob(file_name)

    if len(lst_files) > 1:
        print("Found duplicate {:s} files: {:s}".format(file_name, str(lst_files)))
        return None

    if len(lst_files) == 0:
        #print("No {:s} files found".format(file_name) )
        return None
   
    #print("Found single {:s} for the template {:s}".format(lst_files[0], file_name))
    return lst_files[0]