@echo off
cd /d "C:\Users\danie\Desktop\WICMAD"
R --vanilla --slave -e "source('scripts/arrowhead/profile_arrowhead_analysis.R')"
pause
