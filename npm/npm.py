import json
import os

with open("/home/fusion809/OBS/home:fusion809:arch_extra/arch-wiki-man/package/package.json") as json_file:
    json_data = json.load(json_file)
    deps = json_data["dependencies"]
    for key, value in deps.items():
        print(key)
        os.system("cpobsn" + " " + key)
