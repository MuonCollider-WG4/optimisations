import glob
import subprocess

def is_number(a_string):
    try:
        float(a_string)
        return True
    except Exception:
        return False

def sort_key(fname):
    print(fname)
    a_list = [float(x) for x in reversed(fname.split("_")) if is_number(x)]
    return a_list

def main():
    file_list = glob.glob("2022-04-08_optics-scan/*.png")
    file_list = sorted(file_list, key=sort_key)
    file_list = ["mf://"+x for x in file_list]
    output = subprocess.check_output(["mencoder"]+file_list+
                            ["-mf", "w=800:h=600:fps=1:type=png",
                            "-ovc", "lavc",
                            "-lavcopts", "vcodec=msmpeg4:vbitrate=2000:mbd=2:trell",
                            "-oac", "copy",
                            "-o", "optics.avi"])

if __name__ == "__main__":
    main()
