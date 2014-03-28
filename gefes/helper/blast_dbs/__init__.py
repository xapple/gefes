# Built-in modules #


import json
import os
mod_path = os.path.dirname(os.path.realpath(__file__))
with open(mod_path + "/" + "blast_dbs.json") as j:
    blast_dbs = json.load(j)

blast_header = ["subject","identity","length","mismatches","gaps","start","end","s_start","s_end","e_value","bit"]
