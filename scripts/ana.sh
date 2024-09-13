targetpath="./skimmed/"
mkdir -p $targetpath

#root -l -b -q skimmer/ana.C'("tthh", "'$targetpath'")' &> ./$targetpath/log/log_ana_tthh &
#root -l -b -q skimmer/ana.C'("tth", "'$targetpath'")' &> ./$targetpath/log/log_ana_tthbb &
#root -l -b -q skimmer/ana.C'("ttbbh", "'$targetpath'")' &> ./$targetpath/log/log_ana_ttbbh &
#root -l -b -q skimmer/ana.C'("ttzh", "'$targetpath'")' &> ./$targetpath/log/log_ana_ttzh &
#root -l -b -q skimmer/ana.C'("ttvv", "'$targetpath'")' &> ./$targetpath/log/log_ana_ttvv &
#root -l -b -q skimmer/ana.C'("ttbbv", "'$targetpath'")' &> ./$targetpath/log/log_ana_ttbbv &
#root -l -b -q skimmer/ana.C'("ttbbbb", "'$targetpath'")' &> ./$targetpath/log/log_ana_ttbbbb &
#root -l -b -q skimmer/ana.C'("ttbb", "'$targetpath'")' &> ./$targetpath/log/log_ana_ttbb &
#root -l -b -q skimmer/ana.C'("tttt", "'$targetpath'")' &> ./$targetpath/log/log_ana_tttt &

root -l -b -q skimmer/ana.C'("tthh", "'$targetpath'")' &> ./$targetpath/log/log_ana_tthh &

