#!/usr/bin/env python 
# coding: utf-8

'''
path: untils > createDir
objection:create new dir

创建文件夹的小工具
'''

import os
import argparse

def get_parser():
    parser = argparse.ArgumentParser(description="若指定的文件夹不存在，则创建")
    parser.add_argument('work_dir', metavar='WORK_DIR', nargs = 1, type = str, help='指定文件夹')
    parser.add_argument('new_dir', metavar='NEW_DIR', nargs = '+', type = str, help='指定文件夹')
    return parser

# split evenly into tmp dir
def create_dir(workDir, newDir):
    try:
        if isinstance(newDir, str):
            newDir = [os.path.join(workDir,newDir)]
        else:
            newDir = [ os.path.join(workDir,dir) for dir in newDir]
        for dir in newDir:
            if not os.path.exists(dir):
                os.makedirs(dir)
                print("{} 文件夹创建成功".format(dir))
            else:
                print("{} 文件夹已存在".format(dir))
    except Exception as e:
        print(e)

def main():
    args = vars(get_parser().parse_args())  # 返回对象的属性-属性值的字典对象  return a key-value dict
    work_dir = args['work_dir'][0]
    new_dir = args['new_dir'][0]
    create_dir(work_dir, new_dir)

if __name__ == '__main__':
    main()