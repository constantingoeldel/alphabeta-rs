# WT
cargo run --release --bin metaprofile -- \
-m /mnt/nas/zhilin/others/MA3_new_total_original_methylome/ \
-g ../methylome/gbM_gene_anotation_extract_Arabidopsis.bed \
-o /mnt/extStorage/workingDir/constantin_not_owned_by_postgres/windows/wt_ci \
--name wildtype \
-s 1 -w 5 
alphabeta \
-n ../methylome/nodelist.txt \
-e ../methylome/edgelist.txt \

# # # WT
# # cargo run --release --bin metaprofile -- \
# # -m ../methylome/within_gbM_genes/ \
# # -g ../methylome/gbM_gene_anotation_extract_Arabidopsis.bed \
# # -o ../windows/testing \
# # --name wildtype \
# # -s 1 -w 5 
# # alphabeta \
# # -e ../methylome/edgelist.txt \
# # -n ../methylome/nodelist.txt \



# CMT3
cargo run --release --bin metaprofile -- \
-m /mnt/nas/zhilin/others/constantin-sergio/CMT3/total_original_methylome \
-g ../methylome/gbM_gene_anotation_extract_Arabidopsis.bed \
-o /mnt/extStorage/workingDir/constantin_not_owned_by_postgres/windows/cmt3_ci \
--name cmt3 \
-s 1 -w 5 
alphabeta \
-e /home/constantin/methylome/cmt3_edgelist.txt \
-n /home/constantin/methylome/cmt3_nodelist.txt \

suv 4/5/6
cargo run --release --bin metaprofile -- \
-m /mnt/nas/zhilin/others/constantin-sergio/SUV456/total_original_methylome \
-g ../methylome/gbM_gene_anotation_extract_Arabidopsis.bed \
-o /mnt/extStorage/workingDir/constantin_not_owned_by_postgres/windows/suv_ci \
--name suv \
-s 1 -w 5 
alphabeta \ 
-e /mnt/nas/zhilin/others/constantin-sergio/SUV456/SUV456_edgelist.txt \
-n /mnt/nas/zhilin/others/constantin-sergio/SUV456/SUV456_nodelist.txt \

# ros 
#cargo run --release -- \
#-m /mnt/nas/zhilin/others/constantin-sergio/biostress-data \
#-g ../../methylome/gbM_gene_anotation_extract_Arabidopsis.bed \
#-o ../../windows \
#--name ros \
#-s 1 -w 5 
#alphabeta \
#-e /mnt/nas/zhilin/others/constantin-sergio/biostress-data/edgelist_ros_mock.tsv \
#-n /mnt/nas/zhilin/others/constantin-sergio/biostress-data/nodelist_ros_mock.tsv \

# # nrpe 
# cargo run --release -- \
# -m /mnt/nas/zhilin/others/constantin-sergio/biostress-data \
# -g ../../methylome/gbM_gene_anotation_extract_Arabidopsis.bed \
# -o ../../windows \
# --name nrpe \
# -s 1 -w 5 
# -n /mnt/nas/zhilin/others/constantin-sergio/biostress-data/nodelist_nrpe_mock.tsv \

# -e /mnt/nas/zhilin/others/constantin-sergio/biostress-data/edgelist_nrpe_mock.tsv \

# # # col 
#  cargo run --release -- \
#  -m /mnt/nas/zhilin/others/constantin-sergio/biostress-data \
#  -g ../../methylome/gbM_gene_anotation_extract_Arabidopsis.bed \
#  -o /mnt/extStorage/workingDir/constantin_not_owned_by_postgres/windows/windows_col \
#  --name col \
#  -s 1 -w 5 \
# alphabeta
#  -e /mnt/nas/zhilin/others/constantin-sergio/biostress-data/edgelist_col_mock.tsv \
#  -n ../../methylome/nodelist_col_mock.tsv \

# mods
# cargo run --release -- \
# -m ../../modifications/bed \
# -g ../../methylome/gbM_gene_anotation_extract_Arabidopsis.bed \
# -s 1 -w 1 \
# -o ../../windows \
# --force
# --cutoff-gene-length


# h2az
# cargo run --release -- \
# -m ../../h2az \
# -g ../../methylome/gbM_gene_anotation_extract_Arabidopsis.bed \
# -s 1 -w 5 \
# -o ../../windows \
# --name h2az \
# # --cutoff-gene-length

# # chromatin states
# cargo run --release -- \
# -m ../../chr_states \
# -g ../../methylome/gbM_gene_anotation_extract_Arabidopsis.bed \
# -s 1 -w 1 \
# -o ../../windows \

# # chromatin states within red CS && gbM genes
# cargo run  -- \
# -m ../../chr_states \
# -g ../../methylome/redCS_SPMRs_in_gbM_genes.txt \
# -s 1 -w 1 \
# -o ../../windows \

# # chromatin states within red CS
# cargo run --release -- \
# -m ../../chr_states \
# -g ../../methylome/redCS-SPMRs.bed \
# -s 1 -w 1 \
# -o ../../windows \