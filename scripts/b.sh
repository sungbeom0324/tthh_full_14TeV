targetpath="./skimmed/"
mkdir -p $targetpath

root -l -b -q skimmer/ana_os2l.C'("tthh", "'$targetpath'")' &> ./$targetpath/log/log_b_tthh &


