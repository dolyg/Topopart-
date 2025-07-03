#!/bin/bash

# ����һ��������ִ��ָ�������
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

# ִ��ָ��1
execute_command "./partitioner -t ../testcase/case01 -s ./out1 > ./debuginfo1"

# ִ��ָ��2
execute_command "./partitioner -t ../testcase/case02 -s ./out2 > ./debuginfo2"

# ִ��ָ��3
execute_command "./partitioner -t ../testcase/case03 -s ./out3 > ./debuginfo3"

# ִ��ָ��4
execute_command "./partitioner -t ../testcase/case04 -s ./out4 > ./debuginfo4"