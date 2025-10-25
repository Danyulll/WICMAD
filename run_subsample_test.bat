@echo off
cd /d "C:\Users\danie\Desktop\WICMAD"
R --vanilla --slave -e "source('scripts/arrowhead/test_arrowhead_subsample.R')"
pause
