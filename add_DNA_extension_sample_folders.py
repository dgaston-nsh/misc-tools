#!/usr/bin/env python

# Standard packages
import os
import sys
import glob
import shutil

if __name__ == "__main__":
    rootdir = os.getcwd()
    for file in os.listdir(rootdir):
        d = os.path.join(rootdir, file)
        if os.path.isdir(d):
            new_dir = f"{d}-DNA"
            sys.stdout.write(f"Renaming Directory {d} to {new_dir}\n")
            os.rename(d, new_dir)
