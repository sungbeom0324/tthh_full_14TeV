targetpath="./skimmed/"
mkdir -p $targetpath

root -l -b -q skimmer/ana_os2l.C'("tthh", "'$targetpath'")' &> ./$targetpath/log/log_tthh_osdl &
root -l -b -q skimmer/ana_os2l.C'("tth", "'$targetpath'")' &> ./$targetpath/log/log_tth_osdl &
root -l -b -q skimmer/ana_os2l.C'("ttbbh", "'$targetpath'")' &> ./$targetpath/log/log_ttbbh_osdl &
root -l -b -q skimmer/ana_os2l.C'("ttzh", "'$targetpath'")' &> ./$targetpath/log/log_ttzh_osdl &
root -l -b -q skimmer/ana_os2l.C'("ttvv", "'$targetpath'")' &> ./$targetpath/log/log_ttvv_osdl &
root -l -b -q skimmer/ana_os2l.C'("ttbbv", "'$targetpath'")' &> ./$targetpath/log/log_ttbbv_osdl &
root -l -b -q skimmer/ana_os2l.C'("ttbbbb", "'$targetpath'")' &> ./$targetpath/log/log_ttbbbb_osdl &
root -l -b -q skimmer/ana_os2l.C'("ttbb", "'$targetpath'")' &> ./$targetpath/log/log_ttbb_osdl &
root -l -b -q skimmer/ana_os2l.C'("tttt", "'$targetpath'")' &> ./$targetpath/log/log_tttt_osdl &


