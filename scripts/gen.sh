targetpath="./skimmed/"
mkdir -p $targetpath

root -l -b -q skimmed.gen.C'("tthh", "'$targetpath'")' &> ./$targetpath/log/log_gen_tthh &
root -l -b -q skimmed.gen.C'("tth", "'$targetpath'")' &> ./$targetpath/log/log_gen_tth &
root -l -b -q skimmed.gen.C'("ttbbh", "'$targetpath'")' &> ./$targetpath/log/log_gen_ttbbh &
root -l -b -q skimmed.gen.C'("ttzh", "'$targetpath'")' &> ./$targetpath/log/log_gen_ttzh &
root -l -b -q skimmed.gen.C'("ttvv", "'$targetpath'")' &> ./$targetpath/log/log_gen_ttvv &
root -l -b -q skimmed.gen.C'("ttbbv", "'$targetpath'")' &> ./$targetpath/log/log_gen_ttbbv &
root -l -b -q skimmed.gen.C'("ttbb", "'$targetpath'")' &> ./$targetpath/log/log_gen_ttbb &
root -l -b -q skimmed.gen.C'("ttbbbb", "'$targetpath'")' &> ./$targetpath/log/log_gen_ttbbbb &
root -l -b -q skimmed.gen.C'("tttt", "'$targetpath'")' &> ./$targetpath/log/log_gen_tttt &

