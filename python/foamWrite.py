import json, sys
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile # Integration to OpenFOAM


#Output list from z0 txt file
def write_list(jparams):
    z0_str = []
    for z0 in jparams["z0_list"]:
        z0_str.append(str(z0))
    fin_list = []
    with open(jparams["input"], "r") as f:
        print("-----Reading z0 values from z0.txt file-----")
        for r in f:
            r_stripped = r.rstrip("\n")
            r_dstripped = r_stripped.rstrip("\r")
            if r_dstripped in z0_str:
                fin_list.append(float(r_dstripped))
    return fin_list

# For integration to openfoam (need to add option for writing to another file
def write_to_foam(list_z0, foamLocation, patch_label):
    f = ParsedParameterFile(foamLocation)

    for boundary in f["boundaryField"]:
        if boundary in patch_label:
            f["boundaryField"][boundary]["z0"] = "nonuniform List<Scalar>\n", list_z0
    f.writeFile()

def main():
    # read parameters from the file 'setz0Dict.json'
    try:
        jparams = json.load(open('foamWrite.json'))
    except:
        print("ERROR: something is wrong with the setZ0Dict.json file.")
        sys.exit


    z0_list = write_list(jparams)
    if jparams["writetoBC"] == "on":
        write_to_foam(z0_list, jparams["foamLocation"],jparams["patch_name"] )
    else:
        print("The output is not written to file.")


if __name__ == '__main__':
    main()