#!/bin/bash

autoreconf 
automake --add-missing 
touch NEWS AUTHORS ChangeLog
autoreconf 

