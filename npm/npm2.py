import json
from sys import argv
print(argv[1])
from subprocess import call

for i in argv[1]:
    with open("/home/fusion809/OBS/home:fusion809:arch_extra/nodejs-" + argv[1] + "/src/package/package.json") as json_file:
        json_data = json.load(json_file)
        deps = json_data["dependencies"]
        for key, value in deps.items():
            print(key)
            call(["cpobsn", key])
