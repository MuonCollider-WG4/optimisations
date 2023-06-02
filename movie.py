import os
import glob
import subprocess

def is_number(a_string):
    try:
        float(a_string)
        return True
    except Exception:
        return False

def sort_key(fname):
    #print(fname)
    a_list = [float(x) for x in fname.split("_") if is_number(x)]
    a_list = [x for x in reversed(a_list)]
    return a_list

def main():
    for folder in glob.glob("optics-scan_v16"):
        here = os.getcwd()
        os.chdir(folder)
        file_list = glob.glob("*.png")
        file_list = sorted(file_list, key=sort_key)
        file_list_mencoder = ["mf://"+x for x in file_list]
        print(file_list)
        output = subprocess.check_output(["mencoder"]+file_list_mencoder+
                                ["-mf", "w=800:h=600:fps=1:type=png",
                                "-ovc", "lavc",
                                "-lavcopts", "vcodec=msmpeg4:vbitrate=2000:mbd=2:trell",
                                "-oac", "copy",
                                "-o", "optics_mencoder.avi"])
        file_list_ffmpeg = []
        for f in file_list:
            file_list_ffmpeg += ["-i", f]
        """
        output = subprocess.check_output(["ffmpeg"]+file_list_ffmpeg+
                                ["-r", "1", "optics_ffmpeg.avi"])
        """
        os.chdir(here)

if __name__ == "__main__":
    main()
