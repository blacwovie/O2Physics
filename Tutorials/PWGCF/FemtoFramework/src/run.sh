#!/bin/bash

# in case the central server is not reachable, set SE yourself
# https://alimonitor.cern.ch/stats?page=SE/table
# export alien_CLOSE_SE=ALICE::UPB::EOS

# make sure you can connect to central server
# alien.py

# number of task you want to run, run task 1 by default
Task="${1:-1}"
#这行代码定义了一个名为 Task 的变量，用于指定要运行的任务编号。
#如果用户在运行脚本时没有提供参数，则默认使用任务 1。
#${1:-1} 的意思是：如果脚本运行时第一个参数（$1）不存在，则使用默认值 1。

# config file to use, use config.json by default
ConfigFile="${2:-config.json}"
#用于指定要使用的配置文件。如果没有提供第二个参数（$2），则默认使用 config.json 作为配置文件

# options passed to each workflow
Options=("-b" "--configuration" "json://${ConfigFile}")

# command to be executed
Command="o2-analysistutorial-cf-femtodream-tutorial-${Task} ${Options[*]}"
#这是根据 Task 变量生成的任务名称。${Task} 是任务编号，
#如 o2-analysistutorial-cf-femtodream-tutorial-1 代表第 1 个任务。
# print comand before executing it
echo "$Command"
eval "$Command"

exit 0
