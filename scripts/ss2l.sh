targetpath="./skimmed/"
mkdir -p $targetpath

root -l -b -q skimmer/ana_ss2l.C'("tthh", "'$targetpath'")' &> ./$targetpath/log/log_ss2l_tthh &
root -l -b -q skimmer/ana_ss2l.C'("tth", "'$targetpath'")' &> ./$targetpath/log/log_ss2l_tth &
root -l -b -q skimmer/ana_ss2l.C'("ttbbh", "'$targetpath'")' &> ./$targetpath/log/log_ss2l_ttbbh &
root -l -b -q skimmer/ana_ss2l.C'("ttzh", "'$targetpath'")' &> ./$targetpath/log/log_ss2l_ttzh &
root -l -b -q skimmer/ana_ss2l.C'("ttvv", "'$targetpath'")' &> ./$targetpath/log/log_ss2l_ttvv &
root -l -b -q skimmer/ana_ss2l.C'("ttbbv", "'$targetpath'")' &> ./$targetpath/log/log_ss2l_ttbbv &
root -l -b -q skimmer/ana_ss2l.C'("ttbbbb", "'$targetpath'")' &> ./$targetpath/log/log_ss2l_ttbbbb &
root -l -b -q skimmer/ana_ss2l.C'("ttbb", "'$targetpath'")' &> ./$targetpath/log/log_ss2l_ttbb &
root -l -b -q skimmer/ana_ss2l.C'("tttt", "'$targetpath'")' &> ./$targetpath/log/log_ss2l_tttt &

