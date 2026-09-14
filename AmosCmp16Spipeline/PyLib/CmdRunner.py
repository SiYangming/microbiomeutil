#!/usr/bin/env python

import os, sys

class CmdRunnerException(Exception):

    def __init__(self, cmd, ret):
        self.cmd = cmd
        self.ret = ret
        sys.stderr.write("CMD: " + cmd + " died with ret(" + str(ret) + ")")

    def __repr__(self):
        return("CMD: " + self.cmd + " died with ret(" + str(self.ret) + ")")

    def __str__(self):
        return(repr(self))


def process_cmd(cmd):
    cmd_pipe = os.popen(cmd)
    output = cmd_pipe.readlines()
    ret = cmd_pipe.close()
    if (ret):
        ret >> 8 # shift out to get the exit value
        raise CmdRunnerException(cmd, ret)

    output = "".join(output)
    return(output)

