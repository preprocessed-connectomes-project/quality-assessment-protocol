#!/usr/local/bin/miniconda/envs/qap_fix_test_2/bin/python
# -*- coding: utf-8 -*-
# @Author: oesteban
# @Date:   2015-10-16 13:36:40
# @Last Modified by:   oesteban
# @Last Modified time: 2015-10-16 14:21:04

from qap.cli import QAProtocolCLI

if __name__ == "__main__":
    obj = QAProtocolCLI()
    obj.run()
