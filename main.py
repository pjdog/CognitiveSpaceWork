# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 06:49:02 2022

@author: paulj
"""
#*****************************************************************************#
#                         Main for Sim                                        #
#*****************************************************************************#

import utils.sim_config

def main():
    con_win = utils.sim_config.config_self
    con_win.draw_window()
    

if __name__ == "__main__":
    main()