#!/bin/bash

# 定义一个函数来执行指令并报告结果
execute_command() {
    cmd="$1"
    echo "Executing: $cmd"
    eval $cmd
    if [ $? -eq 0 ]; then
        echo "Command executed successfully: $cmd"
    else
        echo "Error executing command: $cmd"
        exit 1
    fi
}

# 执行指令1
execute_command "./partitioner -t ../testcase/case01 -s ./out1 > ./debuginfo1"

# 执行指令2
execute_command "./partitioner -t ../testcase/case02 -s ./out2 > ./debuginfo2"

# 执行指令3
execute_command "./partitioner -t ../testcase/case03 -s ./out3 > ./debuginfo3"

# 执行指令4
execute_command "./partitioner -t ../testcase/case04 -s ./out4 > ./debuginfo4"