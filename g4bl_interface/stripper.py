import operator
import os


class Stripper:
    def __init__(self):
        self.n_skip_lines = 3
        self.word = 2
        self.value = "5"
        self.operator = operator.eq
        self.fin = None
        self.fout = None
        self.overwrite = True

    def strip(self, filename):
        with open(filename) as self.fin:
            with open(filename+".tmp", "w") as self.fout:
                self.skip_lines()
                for line in self.fin.readlines():
                    line_out = self.strip_line(line)
                    self.fout.write(line_out)
        if self.overwrite:
            print(f"Stripped {filename} and overwriting")
            os.rename(filename+".tmp", filename)


    def skip_lines(self):
        for i in range(self.n_skip_lines):
            line = self.fin.readline()
            self.fout.write(line)

    def strip_line(self, line):
        words = line.split()
        if self.operator(words[self.word], self.value): # operator is true
            return ""
        else:
            return line

    @classmethod
    def G4BLStripper(cls):
        stripper = Stripper()
        stripper.n_skip_lines = 2
        stripper.word = 7
        stripper.value = "2212"
        return stripper

    @classmethod
    def For009Stripper(cls):
        stripper = Stripper()
        stripper.n_skip_lines = 3
        stripper.word = 2
        stripper.value = "5"
        return stripper

