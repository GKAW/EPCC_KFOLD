#!/bin/bash

find ../src/ -name "*f90" | xargs grep -n "\!\!\!"
