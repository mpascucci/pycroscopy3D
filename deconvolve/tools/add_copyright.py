#!/usr/bin/env python3

#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#  
#   Copyright (C) 2018-2020
#    Laboratory of Systems Biology, Department of Cybernetics,
#    School of Science, Tallinn University of Technology
#   This file is part of project: IOCBIO Deconvolve



import os


def main(args):

    excludes = ['.svn' , 'thirdparty', '__init__.py']
    python_types = ['.py', '.pyx']
    c_types = ['.c', '.cpp', '.h', '.hpp']
    file_types = python_types + c_types

    py_cp_right = '\n'
    c_cp_right = '\n/*\n'
    with open(args.copyright, 'r') as f:
        for line in f.readlines():
            py_cp_right += '#  ' + line
            c_cp_right += ' *  ' + line
    py_cp_right += '\n'
    c_cp_right += ' */\n\n'
    
    print(py_cp_right)
    print(c_cp_right)

    py_cp_right = [i+'\n' for i in py_cp_right.split('\n')]
    c_cp_right = [i+'\n' for i in c_cp_right.split('\n')]

    print(py_cp_right)
    print(c_cp_right)

    def proceed(s):
        for ex in excludes:
            if ex in s:
                return False
        return True

    def check_filetype(s):
        for t in file_types:
            if s.endswith(t):
                return True
        return False

    def get_type(s):
        for t in python_types:
            if s.endswith(t): return 'py_type'
        for t in c_types:
            if s.endswith(t): return 'c_type'
        return None

    all_files = []
    for root, dirs, files in os.walk(args.basedir):
        files = [os.path.join(root, f) for f in files]
        for fn in files:
            if proceed(fn) and check_filetype(fn):
                all_files.append(fn)

    for fn in all_files:
        print(fn)
        ftype = get_type(fn)

        # original file
        with open(fn, 'r') as f:
            original = f.readlines()

        if ftype is 'py_type':
            with open(fn, 'w') as f:
                st = 0
                if original[st].startswith('#!'):
                    f.writelines(original[st])
                    st = 1
                f.writelines(py_cp_right)
                f.writelines(original[st:])

        elif ftype is 'c_type':
            with open(fn, 'w') as f:
                f.writelines(c_cp_right)
                f.writelines(original)

        else:
            print('Unknown type', fn)



if __name__ == '__main__':
    import argparse, sys

    parser = argparse.ArgumentParser(description='Script for appending source code with copyright statement')
    parser.add_argument('basedir', type=str, help='Base directory ...')
    parser.add_argument('copyright', type=str, help='The copyright text file')
    # parser.add_argument('-e', '--exclude', type=str, help='Defines which folders/files are not appended with copyright')

    args = parser.parse_args()
    main(args)
