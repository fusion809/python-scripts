import json
from sys import argv
print(argv[1])
from subprocess import call
from subprocess import Popen, PIPE

with open("/home/fusion809/OBS/home:fusion809:arch_extra/nodejs-" + argv[1] + "/src/package/package.json") as json_file:
    json_data = json.load(json_file)
    deps = json_data["dependencies"]
    LEN=len(deps)
    print(LEN)
    i=0
    DEP=list()
    print(DEP)
    for key, value in deps.items():
        print(key)
        DEP.append(key)
        i = i+1
        print(i)
        #call(["cpobsn", key, argv[1]])
    p=Popen(["depends"] + DEP, stdout=PIPE)
    output, err = p.communicate(b"input data that is passed to subprocess' stdin")
    rc = p.returncode
    print(output)
    call(["moddepends", "nodejs-" + argv[1], output])
